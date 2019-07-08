/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#ifndef MC_QUERYING_H_
#define MC_QUERYING_H_

#include <vector>
#include <iostream>
#include <future>
#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "sequence_io.h"
#include "cmdline_utility.h"
#include "query_options.h"
#include "cmdline_utility.h"
#include "sketch_database.h"
#include "candidates.h"

namespace mc {


/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
struct sequence_query
{
    explicit
    sequence_query(query_id qid=0, std::string headerText="",
                   sequence s1= sequence{}, sequence s2 = sequence{}) noexcept
    :
        id{qid}, header(std::move(headerText)),
        seq1(std::move(s1)), seq2(std::move(s2)),
        groundTruth{nullptr}
    {}

    bool empty() const noexcept { return header.empty() || seq1.empty(); }

    query_id id;
    std::string header;
    sequence seq1;
    sequence seq2;  //2nd part of paired-end read
    const taxon* groundTruth;
};



/*************************************************************************//**
 *
 * @param a       : merge sorted ranges in this vector
 * @param offsets : positions where the ranges begin/end
 * @param b       : use this vector as a buffer
 *
 * @details offsets.front() must be 0 and offsets.back() must be a.size()
 *
 *****************************************************************************/
template<class T>
void
merge_sort(std::vector<T>& a, std::vector<size_t>& offsets, std::vector<T>& b) {
    if(offsets.size() < 3) return;
    b.resize(a.size());

    int numChunks = offsets.size()-1;
    for(int s = 1; s < numChunks; s *= 2) {
        for(int i = 0; i < numChunks; i += 2*s) {
            auto begin = offsets[i];
            auto mid = i + s <= numChunks ? offsets[i + s] : offsets[numChunks];
            auto end = i + 2*s <= numChunks ? offsets[i + 2*s] : offsets[numChunks];
            std::merge(a.begin()+begin, a.begin()+mid,
                       a.begin()+mid, a.begin()+end,
                       b.begin()+begin);
        }
        std::swap(a, b);
    }
}



/*************************************************************************//**
 *
 * @brief queries database with batches of reads from one sequence source
 *        produces one match list per sequence
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class LogCallback
>
query_id query_batched(
    const std::string& filename1, const std::string& filename2,
    const database& db, const query_processing_options& opt,
    const query_id& startId,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize,
    LogCallback&& log)
{
    std::mutex finalizeMtx;
    std::atomic<std::int_least64_t> queryLimit{opt.queryLimit};
    std::atomic<query_id> workId{startId};

    // store id and position of most advanced thread
    std::mutex tipMtx;
    query_id qid{startId};
    sequence_pair_reader::stream_positions pos{0,0};

    std::vector<std::future<void>> threads;
    // assign work to threads
    for(int threadId = 0; threadId < opt.numThreads; ++threadId) {
        threads.emplace_back(std::async(std::launch::async, [&] {
            sequence_pair_reader reader{filename1, filename2};

            const auto readSequentially = opt.perThreadSequentialQueries;

            std::vector<sequence_pair_reader::sequence_pair> sequences;
            sequences.reserve(readSequentially);

            match_locations matches;
            match_target_locations matchesBuffer;
            match_target_locations matchesBuffer2;

            while(reader.has_next() && queryLimit > 0) {
                auto batchBuffer = getBuffer();
                bool bufferEmpty = true;

                for(std::size_t i = 0;
                    i < opt.batchSize && queryLimit.fetch_sub(readSequentially) > 0; ++i)
                {
                    query_id myQid;
                    sequence_pair_reader::stream_positions myPos;

                    // get most recent position and query id
                    {
                        std::lock_guard<std::mutex> lock(tipMtx);
                        myQid = qid;
                        myPos = pos;
                    }
                    reader.seek(myPos);
                    if(!reader.has_next()) break;
                    reader.index_offset(myQid);

                    // get work id and skip to this read
                    auto wid = workId.fetch_add(readSequentially);
                    if(myQid != wid) {
                        reader.skip(wid-myQid);
                    }
                    if(!reader.has_next()) break;
                    for(int j = 0; j < readSequentially; ++j) {
                        sequences.emplace_back(reader.next());
                    }

                    // update most recent position and query id
                    myQid = reader.index();
                    myPos = reader.tell();
                    while(myQid > qid) {// reading qid unsafe
                        if(tipMtx.try_lock()) {
                            if(myQid > qid) {// reading qid safe
                                qid = myQid;
                                pos = myPos;
                            }
                            tipMtx.unlock();
                        }
                    }

                    for(auto& seq : sequences) {
                        if(!seq.first.header.empty()) {
                            bufferEmpty = false;
                            matchesBuffer.clear();
                            std::vector<size_t> offsets{0};

                            db.accumulate_matches(seq.first.data, matchesBuffer, offsets);
                            db.accumulate_matches(seq.second.data, matchesBuffer, offsets);

                            merge_sort(matchesBuffer, offsets, matchesBuffer2);

                            matches.clear();
                            for(auto& m : matchesBuffer)
                                matches.emplace_back(db.taxon_of_target(m.tgt), m.win);

                            update(batchBuffer,
                                   sequence_query{seq.first.index,
                                                  std::move(seq.first.header),
                                                  std::move(seq.first.data),
                                                  std::move(seq.second.data)},
                                   matches );
                        }
                    }
                    sequences.clear();
                }
                if(!bufferEmpty) {
                    std::lock_guard<std::mutex> lock(finalizeMtx);
                    finalize(std::move(batchBuffer));
                }
            }
        })); //emplace
    }

    // wait for all threads to finish and catch exceptions
    for(unsigned int threadId = 0; threadId < threads.size(); ++threadId) {
        if(threads[threadId].valid()) {
            try {
                threads[threadId].get();
            }
            catch(file_access_error& e) {
                if(threadId == 0) {
                    log(std::string("FAIL: ") + e.what());
                }
            }
            catch(std::exception& e) {
                log(std::string("FAIL: ") + e.what());
            }
        }
    }

    return qid;
}


/*************************************************************************//**
 *
 * @brief queries database with batches of reads from one sequence source
 *        produces one match list per sequence
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
    template<
            class BufferSource, class BufferClassification, class BufferUpdate, class BufferSink,
            class LogCallback
    >
    query_id query_batched_parallel(
            const std::string& filename1, const std::string& filename2,
            const database& db, const query_processing_options& opt, const classification_options& opt_class,
            const query_id& startId,
            BufferSource&& getBuffer, BufferClassification&& bufclass, BufferUpdate&& update, BufferSink&& finalize,
            LogCallback&& log, int my_id)
    {
        std::mutex results_mapMtx;
        std::mutex addSeqMtx;
        std::mutex writeoutMtx;
        std::map<std::uint_least64_t ,classification_candidates > results_map;
        std::map<std::uint_least64_t ,sequence_query> sequences_map;
        std::mutex finalizeMtx;
        //std::atomic<std::int_least64_t> queryLimit{opt.queryLimit};
        //std::atomic<query_id> workId{startId};

        // store id and position of most advanced thread
        std::mutex tipMtx;
        query_id qid{startId};
        //sequence_pair_reader::stream_positions pos{0,0};

        //unsigned block_size = 50000;
        //std::atomic<unsigned> num_sequences{0};

        std::vector<std::future<void>> threads;

        // Find out number of processes
        int num_procs;
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        //Gather number of items to receive from each rank
        int recvcnts[num_procs];
        int displs[num_procs];
        //int num_items_to_send;

        unsigned queries_added = 0;
        unsigned end_seq = startId + opt.queryLimit;
        // assign work to threads
        int local_thread_id = 0;
        for(int threadId = 0; threadId < opt.numThreads; ++threadId) {
            threads.emplace_back(std::async(std::launch::async, [&] {
                sequence_pair_reader reader{filename1, filename2};

                int my_local_thread_id;

                {
                    std::lock_guard<std::mutex> lock(tipMtx);
                    my_local_thread_id = local_thread_id;
                    //std::cout << "Starting thread" << my_local_thread_id <<std::endl;
                    local_thread_id++;
                }

                //const auto readSequentially = opt.perThreadSequentialQueries;

                std::vector<sequence_pair_reader::sequence_pair> sequences;
                //sequences.reserve(readSequentially);
                reader.skip(startId);

                match_locations matches;
                match_target_locations matchesBuffer;
                match_target_locations matchesBuffer2;

                unsigned current_sequence = startId;

                //while(reader.has_next() && queryLimit > 0 && num_sequences < block_size) {
                while(reader.has_next() && queries_added < opt.queryLimit) {

                    //sequence_pair_reader::sequence_pair current_seq = reader.next();

                    if((current_sequence % opt.numThreads) == (unsigned)my_local_thread_id) {
                        sequences.emplace_back(reader.next());
                        {
                            std::lock_guard<std::mutex> lock(addSeqMtx);
                            queries_added++;
                        }
                    }
                    else {
                        reader.skip(1);
                    }

                    current_sequence++;

                    if(current_sequence >= end_seq) {
                        break;
                    }


                }

                /*{
                    std::lock_guard<std::mutex> lock(writeoutMtx);
                    std::cout << "Number of sequences in thread " << my_local_thread_id << " rank " << my_id << " is "
                              << sequences.size() << std::endl;
                }*/

                for(auto& seq : sequences) {
                    if(!seq.first.header.empty()) {
                        //bufferEmpty = false;
                        matchesBuffer.clear();
                        std::vector<size_t> offsets{0};

                        db.accumulate_matches(seq.first.data, matchesBuffer, offsets);
                        db.accumulate_matches(seq.second.data, matchesBuffer, offsets);

                        // Sort by tgt and win
                        merge_sort(matchesBuffer, offsets, matchesBuffer2);

                        {
                            std::lock_guard<std::mutex> lock(results_mapMtx);

                            sequence_query query = sequence_query{seq.first.index,
                                                                  seq.first.header,
                                                                  seq.first.data,
                                                                  seq.second.data};

                            matches.clear();
                            for(auto& m : matchesBuffer)
                                matches.emplace_back(db.taxon_of_target(m.tgt), m.win);

                            auto batchBuffer = getBuffer();
                            classification_candidates cls = bufclass(batchBuffer, sequence_query{seq.first.index,
                                                                                                 std::move(seq.first.header),
                                                                                                 std::move(seq.first.data),
                                                                                                 std::move(seq.second.data)},
                                                              matches);

                            results_map[seq.first.index] = std::move(cls);//std::move(matches);
                            sequences_map[seq.first.index] = std::move(query);
                        }


                    }

                }

                sequences.clear();
            })); //emplace
        }

        // wait for all threads to finish and catch exceptions
        for(unsigned int threadId = 0; threadId < threads.size(); ++threadId) {
            if(threads[threadId].valid()) {
                try {
                    threads[threadId].get();
                }
                catch(file_access_error& e) {
                    if(threadId == 0) {
                        log(std::string("FAIL: ") + e.what());
                    }
                }
                catch(std::exception& e) {
                    log(std::string("FAIL: ") + e.what());
                }
            }
        }


        qid += results_map.size();

        //std::cout << "Number of items sent from " << my_id << " is: " << results_map.size() << "ID: " << sequences_map.begin()->second.id << ". First: " << sequences_map.begin()->first << ". Content: " << sequences_map.begin()->second.seq2 << std::endl;

        //auto batchBuffer = getBuffer();
        auto current_seq = sequences_map.begin();

        //This way is working fine
        /*
        for (auto item = results_map.begin(); (item != results_map.end()) &&
            current_seq != sequences_map.end(); ++item, ++current_seq) {


            unsigned num_items = item->second.size()*2;

            MPI_Gather(&num_items, 1, MPI_UNSIGNED, recvcnts, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            unsigned total_items = 0;
            //unsigned *items = nullptr;


            unsigned *items_partial = (unsigned *)malloc(num_items* sizeof(unsigned));

            unsigned i = 0;
            for(auto current_item = item->second.begin(); current_item  != item->second.end(); ++current_item) {

                items_partial[i] = current_item->tgt;
                items_partial[i+1] = current_item->win;
                i+=2;


            }

            //int displs[num_procs];
            unsigned *received_locations = nullptr;

            if (my_id == 0) {

                for(int k = 0; k< num_procs; ++k) {
                    total_items += recvcnts[k];

                    if (k == 0) {
                        displs[k] = 0;
                    }
                    else {
                        displs[k] = displs[k-1] + recvcnts[k-1];
                    }
                }

                received_locations = (unsigned *) malloc(total_items*sizeof(unsigned));

            }

            MPI_Gatherv(items_partial, num_items, MPI_UNSIGNED, received_locations, recvcnts, displs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


            if (my_id == 0) {
                match_locations matches;

                for(unsigned j = 0; j< total_items; j+=2) {

                    matches.emplace_back(db.taxon_of_target(received_locations[j]), received_locations[j+1]);
                }

                update(batchBuffer,
                       sequence_query{current_seq->second.id,
                                      std::move(current_seq->second.header),
                                      std::move(current_seq->second.seq1),
                                      std::move(current_seq->second.seq2)},
                       matches );

                //std::cout << "Finalizing in " << my_id << " for seq " << current_seq->second.id <<std::endl;
            }


        }

        if (my_id == 0) {
            //std::cout << "Finalizing in " << my_id << std::endl;
            finalize(std::move(batchBuffer));
        }
        */

        // First gather all data in rank 0
        std::map<std::uint_least64_t, classification_candidates> all_results_map;

        for (auto item = results_map.begin(); (item != results_map.end()) &&
                                              current_seq != sequences_map.end(); ++item, ++current_seq) {


            unsigned num_items = item->second.size()*2;

            MPI_Gather(&num_items, 1, MPI_UNSIGNED, recvcnts, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            unsigned total_items = 0;
            //unsigned *items = nullptr;


            unsigned *items_partial = (unsigned *)malloc(num_items* sizeof(unsigned));

            unsigned i = 0;
            for(auto current_item = item->second.begin(); current_item  != item->second.end(); ++current_item) {

                //items_partial[i] = current_item->tgt;
                //items_partial[i+1] = current_item->win;
                items_partial[i] = current_item->tax->id();
                items_partial[i+1] = current_item->hits;
                i+=2;


            }

            //int displs[num_procs];
            unsigned *received_locations = nullptr;

            if (my_id == 0) {

                for(int k = 0; k< num_procs; ++k) {
                    if (item == results_map.begin()) {
                        std::cout << "Rec items " << recvcnts[k] << std::endl;
                    }

                    total_items += recvcnts[k];

                    if (k == 0) {
                        displs[k] = 0;
                    }
                    else {
                        displs[k] = displs[k-1] + recvcnts[k-1];
                    }
                }

                received_locations = (unsigned *) malloc(total_items*sizeof(unsigned));

            }

            MPI_Gatherv(items_partial, num_items, MPI_UNSIGNED, received_locations, recvcnts, displs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);


            if (my_id == 0) {

                classification_candidates cls;
                candidate_generation_rules rules;

                rules.mergeBelow    = opt_class.lowestRank;
                rules.maxCandidates = opt_class.maxNumCandidatesPerQuery;

                for(unsigned j = 0; j< total_items; j+=2) {
                    match_candidate cand{db.taxon_with_id(received_locations[j]), received_locations[j+1]};
                    cls.insert(cand, db, rules);

                }

                all_results_map[current_seq->first] = std::move(cls);


                /*update(batchBuffer,
                       sequence_query{current_seq->second.id,
                                      std::move(current_seq->second.header),
                                      std::move(current_seq->second.seq1),
                                      std::move(current_seq->second.seq2)},
                       matches );*/

                //std::cout << "Finalizing in " << my_id << " for seq " << current_seq->second.id <<std::endl;
            }


            free(items_partial);
            if (my_id == 0){

                free(received_locations);
            }

        }

        if (my_id == 0) {

            //std::cout << "Number of items in results map: " << all_results_map.size() << std::endl;

            //std::vector<std::future<void>> threads_finalize;
            threads.clear();

            // assign work to threads
            int local_thread_id = 0;
            for(int thread_id = 0; thread_id < opt.numThreads; ++thread_id) {
                threads.emplace_back(std::async(std::launch::async, [&] {
                    bool bufferEmpty = true;
                    //std::cout << "Starting thread" << std::endl;
                    auto batchBuffer = getBuffer();

                    int my_local_thread_id;

                    {
                        std::lock_guard<std::mutex> lock(tipMtx);
                        my_local_thread_id = local_thread_id;
                        //std::cout << "Starting thread" << my_local_thread_id <<std::endl;
                        local_thread_id++;
                    }

                    int current_seq = 0;

                    for(auto & current : all_results_map) {

                        //if(((current_seq % opt.numThreads) == my_local_thread_id) && (!current.second.empty())) {
                        if(((current_seq % opt.numThreads) == my_local_thread_id)) {
                            bufferEmpty = false;

                            update(batchBuffer,
                               sequence_query{sequences_map[current.first].id,//current_seq->second.id,
                                              sequences_map[current.first].header,//std::move(current_seq->second.header),
                                              sequences_map[current.first].seq1, //std::move(current_seq->second.seq1),
                                              sequences_map[current.first].seq2},//std::move(current_seq->second.seq2)},
                               current.second );
                            /*update(batchBuffer, sequences_map[current.first],
                                   current.second );*/
                        }
                        current_seq++;
                    }

                    if(!bufferEmpty) {
                        std::lock_guard<std::mutex> lock(finalizeMtx);
                        finalize(std::move(batchBuffer));
                    }

                }));

            }

            //}


            //std::cout << "Finalizing in " << my_id << std::endl;


            // wait for all threads to finish and catch exceptions
            for(unsigned int threadId = 0; threadId < threads.size(); ++threadId) {
                if(threads[threadId].valid()) {
                    try {
                        threads[threadId].get();
                    }
                    catch(file_access_error& e) {
                        if(threadId == 0) {
                            log(std::string("FAIL: ") + e.what());
                        }
                    }
                    catch(std::exception& e) {
                        log(std::string("FAIL: ") + e.what());
                    }
                }
            }



        }


        all_results_map.clear();
        sequences_map.clear();
        results_map.clear();

        return qid;
    }


/*************************************************************************//**
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class InfoCallback, class ProgressHandler, class LogHandler
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo, ProgressHandler&& showProgress, LogHandler&& log)
{
    const size_t stride = opt.pairing == pairing_mode::files ? 1 : 0;
    const std::string nofile;
    query_id readIdOffset = 0;

    for(size_t i = 0; i < infilenames.size(); i += stride+1) {
        //pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];

        const auto& fname2 = (opt.pairing == pairing_mode::none)
                             ? nofile : infilenames[i+stride];

        if(opt.pairing == pairing_mode::files) {
            showInfo(fname1 + " + " + fname2);
        } else {
            showInfo(fname1);
        }
        showProgress(infilenames.size() > 1 ? i/float(infilenames.size()) : -1);

        readIdOffset = query_batched(fname1, fname2, db, opt, readIdOffset,
                                     std::forward<BufferSource>(bufsrc),
                                     std::forward<BufferUpdate>(bufupdate),
                                     std::forward<BufferSink>(bufsink),
                                     std::forward<LogHandler>(log));
    }
}



/*************************************************************************//**
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink, class InfoCallback
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo)
{
    query_database(infilenames, db, opt,
                   std::forward<BufferSource>(bufsrc),
                   std::forward<BufferUpdate>(bufupdate),
                   std::forward<BufferSink>(bufsink),
                   std::forward<InfoCallback>(showInfo),
                   [] (float p) { show_progress_indicator(std::cerr, p); },
                   [] (const std::string& s) {std::cerr << s << '\n'; }
                   );
}

/*************************************************************************//**
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    receives buffer after batch is finished
 *
 *****************************************************************************/
    template<
            class BufferSource, class BufferClassification, class BufferUpdate, class BufferSink, class InfoCallback, typename identifier
    >
    void query_database_parallel(
            const std::vector<std::string>& infilenames,
            const database& db,
            const query_processing_options& opt, const classification_options& opt_class,
            BufferSource&& bufsrc, BufferClassification&& bufclass, BufferUpdate&& bufupdate, BufferSink&& bufsink,
            InfoCallback&& showInfo, identifier my_id)
    {
        query_database_parallel(infilenames, db, opt, opt_class,
                       std::forward<BufferSource>(bufsrc),
                       std::forward<BufferClassification>(bufclass),
                       std::forward<BufferUpdate>(bufupdate),
                       std::forward<BufferSink>(bufsink),
                       std::forward<InfoCallback>(showInfo),
                       [] (float p) { show_progress_indicator(std::cerr, p); },
                       [] (const std::string& s) {std::cerr << s << '\n'; },
                       my_id
        );
    }

    /*************************************************************************//**
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
    template<
            class BufferSource, class BufferClassification, class BufferUpdate, class BufferSink,
            class InfoCallback, class ProgressHandler, class LogHandler, typename identifier
    >
    void query_database_parallel(
            const std::vector<std::string>& infilenames,
            const database& db,
            const query_processing_options& opt, const classification_options& opt_class,
            BufferSource&& bufsrc, BufferClassification&& bufclass, BufferUpdate&& bufupdate, BufferSink&& bufsink,
            InfoCallback&& showInfo, ProgressHandler&& showProgress, LogHandler&& log, identifier my_id)
    {
        const size_t stride = opt.pairing == pairing_mode::files ? 1 : 0;
        const std::string nofile;
        query_id readIdOffset = 0;
        query_id readIdOffset_two = 0;

        unsigned block_size = opt.queryLimit;

        for(size_t i = 0; i < infilenames.size(); i += stride+1) {
            //pair up reads from two consecutive files in the list
            const auto& fname1 = infilenames[i];

            const auto& fname2 = (opt.pairing == pairing_mode::none)
                                 ? nofile : infilenames[i+stride];

            if(opt.pairing == pairing_mode::files) {
                showInfo(fname1 + " + " + fname2);
            } else {
                showInfo(fname1);
            }
            showProgress(infilenames.size() > 1 ? i/float(infilenames.size()) : -1);

            sequence_pair_reader reader{fname1, fname2};

            unsigned num_seq = 0;

            while(reader.has_next()) {

                reader.next();
                num_seq++;

            }



            if(my_id == 0) {
                std::cout << "Num sequences is: " << num_seq << std::endl;
            }

            for(readIdOffset = 0; readIdOffset < num_seq; readIdOffset += block_size) {

                readIdOffset_two = query_batched_parallel(fname1, fname2, db, opt, opt_class, readIdOffset,
                                                          std::forward<BufferSource>(bufsrc),
                                                          std::forward<BufferClassification>(bufclass),
                                                          std::forward<BufferUpdate>(bufupdate),
                                                          std::forward<BufferSink>(bufsink),
                                                          std::forward<LogHandler>(log),
                                                          my_id);

                if (my_id == 0) {
                    std::cout << "Processed until " << readIdOffset_two << std::endl;
                }

            }
        }
    }


} // namespace mc


#endif
