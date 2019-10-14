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
//#include <omp.h>
#include <unistd.h>
#include <tsl/hopscotch_map.h>

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

    ~sequence_query(){
        header.clear();
        seq1.clear();
        seq2.clear();
        id = -1;
        groundTruth = nullptr;
    }

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
    std::mutex finalizeMtx;

    std::vector<std::vector<classification_candidates>> results_map(opt.numThreads, std::vector<classification_candidates>());
    std::vector<std::vector<sequence_query>> sequences_map(opt.numThreads, std::vector<sequence_query>());

    tsl::hopscotch_map<std::uint_least64_t, classification_candidates> all_results_map(opt.numThreads * opt.queryLimit);

    // store id and position of most advanced thread
    std::mutex tipMtx;
    query_id qid{startId};
    //sequence_pair_reader::stream_positions pos{0,0};

    std::vector<std::future<void>> threads;

    // Find out number of processes
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //Gather number of items to receive from each rank
    int recvcnts[num_procs];
    //int displs[num_procs];
    //int num_items_to_send;


    if (my_id == 0) {

        std::ifstream is {"/proc/self/status"};

        std::string input_line;

        if (is.is_open()) {
            std::string line;
            while (getline(is, line)) {
                if(line.substr(0,6).compare("VmRSS:") == 0 ) {
                    std::cout << "Before start: " << line << std::endl;
                    break;
                }

            }
            is.close();
        }

    }

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

            std::vector<sequence_pair_reader::sequence_pair> sequences;

            unsigned local_startId = startId + (my_local_thread_id*opt.queryLimit);
            unsigned end_seq = local_startId + opt.queryLimit;

            if (my_id == 0) {
                std::cout << "Thread " << my_local_thread_id << " processing from " << local_startId << " to "
                          << end_seq << std::endl;
            }

            reader.skip(local_startId);

            match_locations matches;
            match_target_locations matchesBuffer;
            match_target_locations matchesBuffer2;

            unsigned current_sequence = local_startId;

            sequence_pair_reader::sequence_pair sequence_being_processed;
            //unsigned queries_added = 0;

            //while(reader.has_next() && queryLimit > 0 && num_sequences < block_size) {
            while(reader.has_next() && (current_sequence < end_seq)) {

                sequences.emplace_back(reader.next());

                current_sequence++;

            }

                /*std::cout << "Proc: "<< my_id <<", Thread " << my_local_thread_id << " number of sequences " << sequences.size()
                          << std::endl;*/


            for(auto& seq : sequences) {
                if(!seq.first.header.empty()) {
                    //bufferEmpty = false;
                    matchesBuffer.clear();
                    std::vector<size_t> offsets{0};

                    db.accumulate_matches(seq.first.data, matchesBuffer, offsets);
                    db.accumulate_matches(seq.second.data, matchesBuffer, offsets);

                    // Sort by tgt and win
                    merge_sort(matchesBuffer, offsets, matchesBuffer2);

                    sequence_query query{seq.first.index,
                                                          seq.first.header,
                                                          seq.first.data,
                                                          seq.second.data};

                    matches.clear();
                    for(auto& m : matchesBuffer)
                        matches.emplace_back(db.taxon_of_target(m.tgt), m.win);


                    classification_candidates cls = bufclass(sequence_query{seq.first.index,
                                                                                         seq.first.header,
                                                                                         seq.first.data,
                                                                                         seq.second.data},
                                                      matches);

                    results_map[my_local_thread_id].emplace_back(std::move(cls));
                    sequences_map[my_local_thread_id].emplace_back(std::move(query));

                }

            }

                /*std::cout << "Proc: "<< my_id  << "Thread " << my_local_thread_id << " all sequences processed " << sequences.size()
                          << std::endl;*/

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

    unsigned total_size = 0;

    for(int i = 0; i< opt.numThreads; ++i) {

        total_size += sequences_map[i].size();
    }

    //qid += results_map.size();
    qid += total_size;



    timer time;
    if (my_id == 0) {
        time.start();
        //std::cout << "Getting data: " << std::endl;
    }

    // First gather all data in rank 0
    std::vector<std::uint32_t> items_partial_vector;
    std::uint32_t *items_partial = nullptr;


    if (my_id != 0) {
        for(int i = 0; i< opt.numThreads; ++i) {

            for(unsigned j = 0; j< sequences_map[i].size(); ++j) {

                for(auto &current_class : results_map[i][j]){
                    items_partial_vector.push_back(sequences_map[i][j].id);
                    items_partial_vector.push_back(current_class.tax->id());
                    items_partial_vector.push_back(current_class.hits);

                }

                //items_partial_vector.push_back(current_item.hits);
            }

        }

        items_partial = (std::uint32_t *) malloc(items_partial_vector.size() * sizeof(std::uint32_t));

        std::copy(items_partial_vector.begin(), items_partial_vector.end(), items_partial);

        //items_partial_vector.clear();
    }




    std::uint32_t *received_locations = nullptr;//(unsigned *) calloc(254*num_procs, sizeof(unsigned));

    unsigned num_items = items_partial_vector.size();

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(&num_items, 1, MPI_UINT32_T, recvcnts, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);



    candidate_generation_rules rules;

    rules.mergeBelow    = opt_class.lowestRank;
    rules.maxCandidates = opt_class.maxNumCandidatesPerQuery;


    if (my_id == 0) {

        std::uint_least64_t current_seq_id;

        for(int i = 0; i< opt.numThreads; ++i) {

            for(unsigned j = 0; j< sequences_map[i].size(); ++j) {

                for(auto &current_class : results_map[i][j]){

                    current_seq_id = (std::uint_least64_t)sequences_map[i][j].id;

                    match_candidate cand{db.taxon_with_id(current_class.tax->id()), current_class.hits};

                    auto current_item = all_results_map.find(current_seq_id);
                    if (current_item != all_results_map.end()) {
                        //all_results_map[current_seq_id] = classification_candidates();
                        current_item.value().insert(cand, db, rules);
                    }
                    else {
                        classification_candidates  candidates;
                        candidates.insert(cand, db, rules);
                        //all_results_map[current_seq_id] = std::move(candidates);
                        all_results_map.insert(std::make_pair(current_seq_id, std::move(candidates)));
                    }

                }

                //items_partial_vector.push_back(current_item.hits);
            }

        }

        results_map.clear();

    }

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;


    for(int j = 1; j< num_procs; ++j) {

        if (my_id == 0) {
            std::uint_least64_t current_seq_id;

            received_locations = (std::uint32_t *) malloc(recvcnts[j] * sizeof(std::uint32_t));

            MPI_Recv(received_locations, (int)recvcnts[j], MPI_UINT32_T, j, j+num_procs, MPI_COMM_WORLD, &status);

            //std::cout << "Rank " << my_id << " Recv from rank " << j << ". Items:: " << recvcnts[j] << std::endl;

            for(unsigned k = 0; k< (unsigned)recvcnts[j]; k+=3) {

                current_seq_id = (std::uint_least64_t)received_locations[k];

                match_candidate cand{db.taxon_with_id(received_locations[k + 1]), received_locations[k + 2]};

                auto current_item = all_results_map.find(current_seq_id);
                if (current_item != all_results_map.end()) {
                    //all_results_map[current_seq_id] = classification_candidates();
                    current_item.value().insert(cand, db, rules);
                }
                else {
                    classification_candidates  candidates;
                    candidates.insert(cand, db, rules);
                    //all_results_map[current_seq_id] = std::move(candidates);
                    all_results_map.insert(std::make_pair(current_seq_id, std::move(candidates)));
                }
            }

            //std::cout << "Rank " << my_id << " Size of map " << all_results_map.size() << std::endl;

            free(received_locations);


        }
        else if (my_id == j) {
            //std::cout << "Sending from Rank " << my_id << " to rank 0. Items:: " << num_items << std::endl;
            MPI_Send(items_partial, num_items, MPI_UINT32_T, 0, j+num_procs, MPI_COMM_WORLD);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_id != 0) {
        free(items_partial);
        items_partial_vector.clear();
        results_map.clear();
        //std::cout << "All cleared in rank " << my_id  << std::endl;
    } else {
        time.stop();

        std::cout << "Time involved in receive and process data: " << time.seconds() << "s." << std::endl;

        // From here, multi-thread
        threads.clear();

        // assign work to threads
        local_thread_id = 0;



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

                //int current_seq = 0;

                for(auto &current_sequence : sequences_map[my_local_thread_id]) {

                    //if(((current_seq % opt.numThreads) == my_local_thread_id) && (!current.second.empty())) {
                    //if(((current_seq % opt.numThreads) == my_local_thread_id)) {
                        classification_candidates cls;

                        auto current = all_results_map.find(current_sequence.id);

                        if (current != all_results_map.end()) {
                            bufferEmpty = false;
                            cls = current->second;
                        }

                        update(batchBuffer,
                               sequence_query{current_sequence.id,//current_seq->second.id,
                                              current_sequence.header,
                                              current_sequence.seq1,
                                              current_sequence.seq2},
                               cls);
                    //}
                    //current_seq++;
                }

                if(!bufferEmpty) {
                    std::lock_guard<std::mutex> lock(finalizeMtx);
                    finalize(std::move(batchBuffer));
                }

            }));

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


        // Common code
        //free(received_locations);

    }


    items_partial_vector.clear();
    sequences_map.clear();
    /*for(int i = 0; i< opt.numThreads; ++i) {
        sequences_map[i].clear();
        results_map[i].clear();
    }*/

    all_results_map.clear();


    if (my_id == 0) {
        std::ifstream is{"/proc/self/status"};

        if (is.is_open()) {
            std::string line;
            while (getline(is, line)) {
                if (line.substr(0, 6).compare("VmRSS:") == 0) {
                    std::cout << "Final memory: " << line << std::endl;
                    break;
                }

            }
            is.close();
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
template<class BufferSource, class BufferClassification, class BufferUpdate, class BufferSink, class LogCallback>
query_id query_batched_parallel2(
            const std::string& filename1, const std::string& filename2,
            const database& db, const query_processing_options& opt, const classification_options& opt_class,
            const query_id& startId,
            BufferSource&& getBuffer, BufferClassification&& bufclass, BufferUpdate&& update, BufferSink&& finalize,
            LogCallback&& log, int my_id)
{
    std::mutex finalizeMtx;

    std::vector<std::vector<classification_candidates>> results_map(opt.numThreads, std::vector<classification_candidates>());
    std::vector<std::vector<sequence_query>> sequences_map(opt.numThreads, std::vector<sequence_query>());

    tsl::hopscotch_map<std::uint_least64_t, classification_candidates> all_results_map(opt.numThreads * opt.queryLimit);

    // store id and position of most advanced thread
    std::mutex tipMtx;
    query_id qid{startId};
    //sequence_pair_reader::stream_positions pos{0,0};

    std::vector<std::future<void>> threads;

    // Find out number of processes
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

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

            std::vector<sequence_pair_reader::sequence_pair> sequences;

            unsigned local_startId = startId + (my_local_thread_id*opt.queryLimit);
            unsigned end_seq = local_startId + opt.queryLimit;

            if (my_id == 0) {
                std::cout << "Thread " << my_local_thread_id << " processing from " << local_startId << " to "
                          << end_seq << std::endl;
            }

            reader.skip(local_startId);

            match_locations matches;
            match_target_locations matchesBuffer;
            match_target_locations matchesBuffer2;

            unsigned current_sequence = local_startId;

            sequence_pair_reader::sequence_pair sequence_being_processed;
            //unsigned queries_added = 0;

            //while(reader.has_next() && queryLimit > 0 && num_sequences < block_size) {
            while(reader.has_next() && (current_sequence < end_seq)) {

                sequences.emplace_back(reader.next());

                current_sequence++;

            }

            for(auto& seq : sequences) {
                if(!seq.first.header.empty()) {
                    //bufferEmpty = false;
                    matchesBuffer.clear();
                    std::vector<size_t> offsets{0};

                    db.accumulate_matches(seq.first.data, matchesBuffer, offsets);
                    db.accumulate_matches(seq.second.data, matchesBuffer, offsets);

                    // Sort by tgt and win
                    merge_sort(matchesBuffer, offsets, matchesBuffer2);

                    sequence_query query{seq.first.index,
                                         seq.first.header,
                                         seq.first.data,
                                         seq.second.data};

                    matches.clear();
                    for(auto& m : matchesBuffer)
                        matches.emplace_back(db.taxon_of_target(m.tgt), m.win);


                    classification_candidates cls = bufclass(sequence_query{seq.first.index,
                                                                            seq.first.header,
                                                                            seq.first.data,
                                                                            seq.second.data},
                                                             matches);

                    results_map[my_local_thread_id].emplace_back(std::move(cls));
                    sequences_map[my_local_thread_id].emplace_back(std::move(query));

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

    unsigned total_size = 0;

    for(int i = 0; i< opt.numThreads; ++i) {

        total_size += sequences_map[i].size();
    }

    //qid += results_map.size();
    qid += total_size;



    timer time;
    if (my_id == 0) {
        time.start();
        //std::cout << "Getting data: " << std::endl;
    }

    // Senders and receivers
    std::set<int> senders;
    std::set<int> receivers;

    for(int i = 0; i < num_procs; ++i) {
        if( i % 2 == 0) {
            receivers.insert(i);
        }
        else {
            senders.insert(i);
        }
    }

    // Rules
    candidate_generation_rules rules;

    rules.mergeBelow    = opt_class.lowestRank;
    rules.maxCandidates = opt_class.maxNumCandidatesPerQuery;

    // Where to store data to send/rec
    std::vector<std::uint32_t> items_partial_vector;
    std::uint32_t *items_partial = nullptr;

    MPI_Status status;

    for(int k = num_procs; k > 1; k = k/2) {

        std::vector<int> senders_to_process;
        std::vector<int> receivers_to_process;

        // Iterate over senders and receivers
        for(auto current_sender = senders.begin(), current_receiver = receivers.begin(); (current_sender != senders.end()) && (current_receiver != receivers.end()); current_sender++, current_receiver++) {

            int sender = *current_sender;
            int receiver = *current_receiver;

            std::uint32_t num_items_to_send;

            if (my_id == receiver) {

                std::uint_least64_t current_seq_id;

                // Only in first iteration
                if (all_results_map.empty()) {
                    for (int i = 0; i < opt.numThreads; ++i) {

                        for (unsigned j = 0; j < sequences_map[i].size(); ++j) {

                            for (auto &current_class : results_map[i][j]) {

                                current_seq_id = (std::uint_least64_t) sequences_map[i][j].id;

                                match_candidate cand{db.taxon_with_id(current_class.tax->id()), current_class.hits};

                                auto current_item = all_results_map.find(current_seq_id);
                                if (current_item != all_results_map.end()) {
                                    //all_results_map[current_seq_id] = classification_candidates();
                                    current_item.value().insert(cand, db, rules);
                                } else {
                                    classification_candidates candidates;
                                    candidates.insert(cand, db, rules);
                                    //all_results_map[current_seq_id] = std::move(candidates);
                                    all_results_map.insert(std::make_pair(current_seq_id, std::move(candidates)));
                                }

                            }

                            //items_partial_vector.push_back(current_item.hits);
                        }

                    }

                    results_map.clear();
                }


                // Rec number of items to obtain from sender
                MPI_Recv(&num_items_to_send, 1, MPI_UINT32_T, sender, receiver+num_procs, MPI_COMM_WORLD, &status);

                // Reserve memory
                auto *received_locations = (std::uint32_t *) malloc(num_items_to_send * sizeof(std::uint32_t));

                // Receive items
                MPI_Recv(received_locations, (int)num_items_to_send, MPI_UINT32_T, sender, sender+num_procs, MPI_COMM_WORLD, &status);
                //std::cout << "Rank " << my_id << " Recv from rank " << sender << ". Items:: " << num_items_to_send << std::endl;

                // Process received data
                for(unsigned w = 0; w < (unsigned)num_items_to_send; w+=3) {

                    current_seq_id = (std::uint_least64_t)received_locations[w];

                    match_candidate cand{db.taxon_with_id(received_locations[w + 1]), received_locations[w + 2]};

                    auto current_item = all_results_map.find(current_seq_id);
                    if (current_item != all_results_map.end()) {
                        //all_results_map[current_seq_id] = classification_candidates();
                        current_item.value().insert(cand, db, rules);
                    }
                    else {
                        classification_candidates  candidates;
                        candidates.insert(cand, db, rules);
                        //all_results_map[current_seq_id] = std::move(candidates);
                        all_results_map.insert(std::make_pair(current_seq_id, std::move(candidates)));
                    }
                }
            }

            else if (my_id == sender) {

                // In first iteration add items stored locally
                if (all_results_map.empty()) {
                    for(int i = 0; i< opt.numThreads; ++i) {

                        for(unsigned j = 0; j< sequences_map[i].size(); ++j) {

                            for(auto &current_class : results_map[i][j]){
                                items_partial_vector.push_back(sequences_map[i][j].id);
                                items_partial_vector.push_back(current_class.tax->id());
                                items_partial_vector.push_back(current_class.hits);

                            }

                            //items_partial_vector.push_back(current_item.hits);
                        }

                    }
                }
                else { // In the rest of iterations, store items that are already in the hashmap

                    for (auto current_entry = all_results_map.begin(); current_entry != all_results_map.end(); ++current_entry) {

                        std::uint32_t current_seq_id = current_entry.key();
                        for(auto current_entry_in_this_seq = current_entry.value().begin(); current_entry_in_this_seq != current_entry.value().end(); ++current_entry_in_this_seq) {
                            items_partial_vector.push_back(current_seq_id);
                            items_partial_vector.push_back(current_entry_in_this_seq->tax->id());
                            items_partial_vector.push_back(current_entry_in_this_seq->hits);
                        }


                    }

                    all_results_map.clear();
                }

                items_partial = (std::uint32_t *) malloc(items_partial_vector.size() * sizeof(std::uint32_t));

                std::copy(items_partial_vector.begin(), items_partial_vector.end(), items_partial);

                num_items_to_send = items_partial_vector.size();

                MPI_Send(&num_items_to_send, 1, MPI_UINT32_T, receiver, receiver+num_procs, MPI_COMM_WORLD);

                MPI_Send(items_partial, num_items_to_send, MPI_UINT32_T, receiver, sender+num_procs, MPI_COMM_WORLD);

                // Free memory
                free(items_partial);
                items_partial_vector.clear();

            }

            senders_to_process.emplace_back(sender);
            receivers_to_process.emplace_back(receiver);
            //}

        }

        MPI_Barrier(MPI_COMM_WORLD);

        for (int current_item : senders_to_process) {

            auto sender_to_delete = senders.find(current_item);

            if (sender_to_delete != senders.end()) {
                senders.erase(sender_to_delete);
            }

        }

        bool do_delete = false;

        for(int current_item: receivers_to_process) {

            if (do_delete) {

                auto receiver_to_delete = receivers.find(current_item);

                if (receiver_to_delete != receivers.end()) {

                    receivers.erase(receiver_to_delete);
                    senders.insert(current_item);
                    do_delete = false;

                }
            }
            else {
                do_delete = true;
            }


        }

        senders_to_process.clear();
        receivers_to_process.clear();

    }

    MPI_Barrier(MPI_COMM_WORLD);



    if (my_id == 0) {
        time.stop();

        std::cout << "Time involved in receive and process data: " << time.seconds() << "s." << std::endl;

        // From here, multi-thread
        threads.clear();

        // assign work to threads
        local_thread_id = 0;

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

                //int current_seq = 0;

                for(auto &current_sequence : sequences_map[my_local_thread_id]) {

                    //if(((current_seq % opt.numThreads) == my_local_thread_id) && (!current.second.empty())) {
                    //if(((current_seq % opt.numThreads) == my_local_thread_id)) {
                    classification_candidates cls;

                    auto current = all_results_map.find(current_sequence.id);

                    if (current != all_results_map.end()) {
                        bufferEmpty = false;
                        cls = current->second;
                    }

                    update(batchBuffer,
                           sequence_query{current_sequence.id,//current_seq->second.id,
                                          current_sequence.header,
                                          current_sequence.seq1,
                                          current_sequence.seq2},
                           cls);
                    //}
                    //current_seq++;
                }

                if(!bufferEmpty) {
                    std::lock_guard<std::mutex> lock(finalizeMtx);
                    finalize(std::move(batchBuffer));
                }

            }));

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


        // Common code
        //free(received_locations);

    }


    items_partial_vector.clear();
    sequences_map.clear();
    /*for(int i = 0; i< opt.numThreads; ++i) {
        sequences_map[i].clear();
        results_map[i].clear();
    }*/

    all_results_map.clear();


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

        /*std::map<std::uint_least64_t ,classification_candidates> *results_map = new std::map<std::uint_least64_t ,classification_candidates>();
        std::map<std::uint_least64_t ,sequence_query> *sequences_map = new std::map<std::uint_least64_t ,sequence_query>();
        std::map<std::uint_least64_t, classification_candidates> *all_results_map = new std::map<std::uint_least64_t, classification_candidates>();*/

        std::cout << "Starting to process files in rank: " << my_id << std::endl;

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
            //showProgress(infilenames.size() > 1 ? i/float(infilenames.size()) : -1);

            sequence_pair_reader reader{fname1, fname2};

            unsigned num_seq = 0;

            while(reader.has_next()) {

                reader.next();
                num_seq++;

            }

            reader.close();


            if(my_id == 0) {
                std::cout << "Num sequences is: " << num_seq << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for(readIdOffset = 0; readIdOffset < num_seq; readIdOffset += (block_size * opt.numThreads)) {

                /*readIdOffset_two = query_batched_parallel_basic(fname1, fname2, db, opt, opt_class, readIdOffset,
                                                          std::forward<BufferSource>(bufsrc),
                                                          std::forward<BufferClassification>(bufclass),
                                                          std::forward<BufferUpdate>(bufupdate),
                                                          std::forward<BufferSink>(bufsink),
                                                          std::forward<LogHandler>(log),
                                                          my_id);//, results_map, sequences_map, all_results_map);*/

                readIdOffset_two = query_batched_parallel2(fname1, fname2, db, opt, opt_class, readIdOffset,
                                                                std::forward<BufferSource>(bufsrc),
                                                                std::forward<BufferClassification>(bufclass),
                                                                std::forward<BufferUpdate>(bufupdate),
                                                                std::forward<BufferSink>(bufsink),
                                                                std::forward<LogHandler>(log),
                                                                my_id);//, results_map, sequences_map, all_results_map);



                if (my_id == 0) {
                    std::cout << "Processed until " << readIdOffset_two << std::endl;
                }

            }
        }
    }


} // namespace mc


#endif
