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

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>

#include "timer.h"
#include "args_handling.h"
#include "args_parser.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "io_error.h"
#include "io_options.h"
#include "config.h"
#include "sequence_io.h"
#include "taxonomy_io.h"


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;



/*************************************************************************//**
 *
 * @brief database creation parameters
 *
 *****************************************************************************/
struct build_options
{
    int kmerlen = 16;
    int sketchlen = 16;
    int winlen = 128;
    int winstride = 113;

    float maxLoadFactor = -1;           //< 0 : use database default
    int maxLocationsPerFeature = database::max_supported_locations_per_feature();
    bool removeOverpopulatedFeatures = true;

    taxon_rank removeAmbigFeaturesOnRank = taxon_rank::none;
    int maxTaxaPerFeature = 1;

    taxonomy_options taxonomy;
    bool resetParents = false;

    info_level infoLevel = info_level::moderate;

    string dbfile;
    std::vector<string> infiles;
};



/*************************************************************************//**
 *
 * @brief command line args -> database creation parameters
 *
 *****************************************************************************/
build_options
get_build_options(const args_parser& args, const build_options& defaults = {})
{
    build_options opt;

    opt.dbfile = database_name(args);

    opt.infiles = sequence_filenames(args);

    if(args.contains("silent")) {
        opt.infoLevel = info_level::silent;
    } else if(args.contains("verbose")) {
        opt.infoLevel = info_level::verbose;
    }

    opt.kmerlen   = args.get<int>({"kmerlen"}, defaults.kmerlen);
    opt.sketchlen = args.get<int>({"sketchlen"}, defaults.sketchlen);
    opt.winlen    = args.get<int>({"winlen"}, defaults.winlen);
    opt.winstride = args.get<int>({"winstride"}, opt.winlen - opt.kmerlen + 1);

    opt.maxLoadFactor = args.get<float>({"max-load-fac", "max_load_fac",
                                         "maxloadfac"},
                                        defaults.maxLoadFactor);

    opt.maxLocationsPerFeature = args.get<int>({"max-locations-per-feature",
                                                "max_locations_per_feature" },
                                                 defaults.maxLocationsPerFeature);

    opt.removeOverpopulatedFeatures = args.contains({"remove-overpopulated-features",
                                                     "remove_overpopulated_features" });

    opt.removeAmbigFeaturesOnRank = taxonomy::rank_from_name(
        args.get<string>({"remove-ambig-features",
                               "remove_ambig_features"},
                    taxonomy::rank_name(defaults.removeAmbigFeaturesOnRank)));

    opt.maxTaxaPerFeature = args.get<int>({"max-ambig-per-feature",
                                           "max_ambig_per_feature"},
                                            defaults.maxTaxaPerFeature);

    opt.resetParents = args.contains({"reset-parents", "reset_parents"});

    opt.taxonomy = get_taxonomy_options(args);

    return opt;
}



/*************************************************************************//**
 *
 * @brief extract db creation parameters
 *
 *****************************************************************************/
build_options
get_build_options_from_db(const database& db)
{
    build_options opt;

    opt.kmerlen   = db.target_sketcher().kmer_size();
    opt.sketchlen = db.target_sketcher().sketch_size();
    opt.winlen    = db.target_window_size();
    opt.winstride = db.target_window_stride();

    opt.maxLoadFactor = db.max_load_factor();

    opt.maxLocationsPerFeature = db.max_locations_per_feature();

    return opt;
}



/*************************************************************************//**
 *
 * @brief Alternative way of providing taxonomic mappings.
 *        The input file must be a text file with each line in the format:
 *        accession  accession.version taxid gi
 *        (like the NCBI's *.accession2version files)
 *
 *****************************************************************************/
void rank_targets_with_mapping_file(database& db,
                                    std::set<const taxon*>& targetTaxa,
                                    const string& mappingFile)
{
    if(targetTaxa.empty()) return;

    std::ifstream is {mappingFile};
    if(!is.good()) return;

    const auto fsize = file_size(mappingFile);
    auto nextStatStep = fsize / 1000;
    auto nextStat = nextStatStep;

    // Find out process rank
    int my_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    bool showProgress = (fsize > 100000000) && (my_id == 0);

    cout << "Try to map sequences to taxa in rank " << my_id << " using '" << mappingFile
         << "' (" << std::max(std::streamoff(1),
                              fsize/(1024*1024)) << " MB)" << endl;

    if(showProgress) show_progress_indicator(cout, 0);

    string acc;
    string accver;
    std::uint64_t taxid;
    string gi;

    //skip header
    getline(is, acc);
    acc.clear();

    while(is >> acc >> accver >> taxid >> gi) {
        //target in database?
        //accession.version is the default
        const taxon* tax = db.taxon_with_name(accver);

        if(!tax) {
            tax = db.taxon_with_similar_name(acc);
            if(!tax) tax = db.taxon_with_name(gi);
        }

        //if in database then set parent
        if(tax) {
            auto i = targetTaxa.find(tax);
            if(i != targetTaxa.end()) {
                db.reset_parent(*tax, taxid);
                targetTaxa.erase(i);
                if(targetTaxa.empty()) break;
            }
        }

        if(showProgress) {
            auto pos = is.tellg();
            if(pos >= nextStat) {
                show_progress_indicator(cout, pos / float(fsize));
                nextStat = pos + nextStatStep;
            }
        }
    }

    if(showProgress) clear_current_line(cout);
}



/*************************************************************************//**
 *
 * @return target taxa that have no parent taxon assigned
 *
 *****************************************************************************/
std::set<const taxon*>
unranked_targets(const database& db)
{
    auto res = std::set<const taxon*>{};

    for(const auto& tax : db.target_taxa()) {
        if(!tax.has_parent()) res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @return all target taxa
 *
 *****************************************************************************/
std::set<const taxon*>
all_targets(const database& db)
{
    auto res = std::set<const taxon*>{};

    for(const auto& tax : db.target_taxa()) {
        res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @brief try to assign parent taxa to target taxa using mapping files
 *
 *****************************************************************************/
void try_to_rank_unranked_targets(database& db, const build_options& opt, int my_id = 0)
{
    std::set<const taxon*> unranked;
    if(opt.resetParents)
        unranked = all_targets(db);
    else
        unranked = unranked_targets(db);

    if(!unranked.empty()) {
        if((opt.infoLevel != info_level::silent) && (my_id == 0)) {
            cout << unranked.size()
                 << " targets are unranked." << endl;
        }

        for(const auto& file : opt.taxonomy.mappingPostFiles) {
            rank_targets_with_mapping_file(db, unranked, file);
            if(unranked.empty()) break;
        }
    }

    unranked = unranked_targets(db);
    if(opt.infoLevel != info_level::silent) {
        if(unranked.empty()) {
            cout << "All targets are ranked." << endl;
        }
        else {
            cout << unranked.size()
                 << " targets remain unranked." << (*unranked.begin())->name() << endl;
        }
    }
}



/*************************************************************************//**
 *
 * @brief look up taxon id based on an identifier (accession number etc.)
 *
 *****************************************************************************/
taxon_id find_taxon_id(
    const std::map<string,taxon_id>& name2tax, 
    const string& name)
{
    if(name2tax.empty()) return taxonomy::none_id();
    if(name.empty()) return taxonomy::none_id();

    //try to find exact match
    auto i = name2tax.find(name);
    if(i != name2tax.end()) return i->second;

    //find nearest match
    i = name2tax.upper_bound(name);
    if(i == name2tax.end()) return taxonomy::none_id();

    //if nearest match contains 'name' as prefix -> good enough
    //e.g. accession vs. accession.version
    if(i->first.compare(0,name.size(),name) != 0) return taxonomy::none_id();
    return i->second;
}



/*************************************************************************//**
 *
 * @brief adds reference sequences from *several* files to database
 *
 *****************************************************************************/
void add_targets_to_database(database& db,
    const std::vector<string>& infiles,
    const std::map<string,taxon_id>& sequ2taxid,
    info_level infoLvl = info_level::moderate,
    int my_id = 0,
    int num_procs = 1)
{
    int n = infiles.size();
    int i = 0;
    //unsigned long current_sequence = 0;

    //read sequences
    for(const auto& filename : infiles) {
        if(infoLvl == info_level::verbose) {
            cout << "  " << filename << " in rank " << my_id <<" ... " << flush;
        } else if((infoLvl != info_level::silent) && (my_id == 0)) {
            show_progress_indicator(cout, i/float(n));
        }

        try {
            auto reader = make_sequence_reader(filename);
            while(reader->has_next()) {

                auto sequ = reader->next();
                if(!sequ.data.empty()) {
                    auto seqId  = extract_accession_string(sequ.header);
                    auto fileId = extract_accession_string(filename);

                    //make sure sequence id is not empty,
                    //use entire header if necessary
                    if(seqId.empty()) seqId = sequ.header;

                    taxon_id parentTaxId = find_taxon_id(sequ2taxid, seqId);

                    if(parentTaxId == taxonomy::none_id())
                        parentTaxId = find_taxon_id(sequ2taxid, fileId);

                    if(parentTaxId == taxonomy::none_id())
                        parentTaxId = extract_taxon_id(sequ.header);

                    if(infoLvl == info_level::verbose) {
                        cout << "[" << seqId;
                        if(parentTaxId > 0) cout << ":" << parentTaxId;
                        cout << "] ";
                    }

                    //try to add to database if sequence belong to this rank
                    //if ((current_sequence % num_procs) == (unsigned long)my_id) {
                        bool added = db.add_target_distributed(
                                sequ.data, seqId, parentTaxId,
                                database::file_source{filename, reader->index()},
                                num_procs, my_id);

                        if (infoLvl == info_level::verbose && !added) {
                            cout << seqId << " not added to database" << endl;
                        }
                    //}
                    //current_sequence++;
                }
            }
            if((infoLvl == info_level::verbose) && (my_id == 0)) {
                cout << "done." << endl;
            }
        }
        catch(database::target_limit_exceeded_error&) {
            cout << endl;
            cerr << "Reached maximum number of targets per database ("
                 << db.max_target_count() << ").\n"
                 << "See 'README.md' on how to compile MetaCache with "
                 << "support for databases with more reference targets.\n"
                 << endl;
            break;
        }
        catch(std::exception& e) {
            if(infoLvl == info_level::verbose) {
                cout << "FAIL: " << e.what() << endl;
            }
        }

        ++i;
    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build
 *
 *****************************************************************************/
void prepare_database(database& db, const build_options& opt)
{
    if(opt.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
    }

    if(opt.maxLoadFactor > 0.4 && opt.maxLoadFactor < 0.99) {
        db.max_load_factor(opt.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << opt.maxLoadFactor << endl;
    }

    if(!opt.taxonomy.path.empty()) {
        db.reset_taxa_above_sequence_level(
            make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel) );

        if(opt.infoLevel != info_level::silent) {
            cout << "Taxonomy applied to database." << endl;
        }
    }

    if(db.non_target_taxon_count() < 1 && opt.infoLevel != info_level::silent) {
        cout << "The datbase doesn't contain a taxonomic hierarchy yet.\n"
                  << "You can add one or update later via:\n"
                  << "   ./metacache add <database> -taxonomy <directory>"
                  << endl;
    }

    if(opt.removeAmbigFeaturesOnRank != taxon_rank::none &&
       opt.infoLevel != info_level::silent)
    {
        if(db.non_target_taxon_count() > 1) {
            cout << "Ambiguous features on rank "
                      << taxonomy::rank_name(opt.removeAmbigFeaturesOnRank)
                      << " will be removed afterwards.\n";
        } else {
            cout << "Could not determine amiguous features "
                      << "due to missing taxonomic information.\n";
        }
    }
}



/*************************************************************************//**
 *
 * @brief database features post-processing
 *
 *****************************************************************************/
void post_process_features(database& db, const build_options& opt)
{
    const bool notSilent = opt.infoLevel != info_level::silent;

    if(opt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();
        auto maxlpf = db.max_locations_per_feature() - 1;
        if(maxlpf > 0) { //always keep buckets with size 1
            if(notSilent) {
                cout << "\nRemoving features with more than "
                     << maxlpf << " locations... " << flush;
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if(notSilent) {
                cout << rem << " of " << old << " removed." << endl;
                if(rem != old) print_content_properties(db);
            }
        }
    }

    if(opt.removeAmbigFeaturesOnRank != taxon_rank::none &&
        db.non_target_taxon_count() > 1)
    {
        if(notSilent) {
            cout << "\nRemoving ambiguous features on rank "
                 << taxonomy::rank_name(opt.removeAmbigFeaturesOnRank)
                 << "... " << flush;
        }

        auto old = db.feature_count();
        auto rem = db.remove_ambiguous_features(opt.removeAmbigFeaturesOnRank,
                                                opt.maxTaxaPerFeature);

        if(notSilent) {
            cout << rem << " of " << old << "." << endl;
            if(rem != old) print_content_properties(db);
        }

    }
}

/*************************************************************************//**
 *
 * @brief database features post-processing
 *
 *****************************************************************************/
    void post_process_features_distributed(database& db, const build_options& opt, std::uint32_t overpopulated_keys[], std::uint32_t num_overpopulated_keys)
    {
        const bool notSilent = opt.infoLevel != info_level::silent;

        if(opt.removeOverpopulatedFeatures) {

            auto rem = db.remove_features_with_more_locations_than_distributed(overpopulated_keys, num_overpopulated_keys);

            if(notSilent) {
                cout << rem << " features removed." << endl;
            }

        }
/*
        if(opt.removeAmbigFeaturesOnRank != taxon_rank::none &&
           db.non_target_taxon_count() > 1)
        {
            if(notSilent) {
                cout << "\nRemoving ambiguous features on rank "
                     << taxonomy::rank_name(opt.removeAmbigFeaturesOnRank)
                     << "... " << flush;
            }

            auto old = db.feature_count();
            auto rem = db.remove_ambiguous_features(opt.removeAmbigFeaturesOnRank,
                                                    opt.maxTaxaPerFeature);

            if(notSilent) {
                cout << rem << " of " << old << "." << endl;
                if(rem != old) print_content_properties(db);
            }

        }*/
    }

/*************************************************************************//**
 *
 * @brief prepares datbase for build, adds targets and writes database to disk
 *
 *****************************************************************************/
void add_to_database(database& db, const build_options& opt, int my_id, int num_procs)
{
    prepare_database(db, opt);

    const bool notSilent = (opt.infoLevel != info_level::silent) && (my_id == 0);
    if(notSilent) print_static_properties(db);

    timer time;
    time.start();

    if(!opt.infiles.empty()) {
        if(notSilent) cout << "Processing reference sequences." << endl;

        const auto initNumTargets = db.target_count();

        auto taxonMap = make_sequence_to_taxon_id_map(
                            opt.taxonomy.mappingPreFiles,
                            opt.infiles, opt.infoLevel);

        add_targets_to_database(db, opt.infiles, taxonMap, opt.infoLevel, my_id, num_procs);

        if(notSilent) {
            clear_current_line(cout);
            cout << "Added "
                 << (db.target_count() - initNumTargets) << " reference sequences "
                 << "in " << time.seconds() << " s" << endl;

            print_content_properties(db);
        }

    }

    try_to_rank_unranked_targets(db, opt, my_id);

    cout << "[JMABUIN] Waiting at proc " << my_id  << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if(opt.removeOverpopulatedFeatures) {
        cout << "[JMABUIN] At proc " << my_id << " :: Before Creating array :: " << endl;

        // Obtain features and number of items from all ranks
        std::uint32_t *items_size = (std::uint32_t *) malloc(db.feature_count() * 2 * sizeof(std::uint32_t));
        cout << "[JMABUIN] At proc " << my_id << " :: Before get_keys_num_items :: " << endl;
        db.get_keys_num_items(items_size);

        //Gather number of items to receive from each rank
        std::uint32_t recvcnts[num_procs];
        std::uint32_t num_items_to_send;

        num_items_to_send = db.feature_count() * 2;

        cout << "[JMABUIN] At proc " << my_id << " :: Before Gather :: " << num_items_to_send << endl;

        MPI_Gather(&num_items_to_send, 1, MPI_UINT32_T, recvcnts, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

        cout << "[JMABUIN] At proc " << my_id << " :: After Gather :: Sent: " <<num_items_to_send << endl;

        if (my_id == 0) {

            for(int i = 0; i< num_procs; ++i) {
                cout << "[JMABUIN] At proc " << my_id << " :: received: " <<recvcnts[i] << endl;

            }

        }

        std::uint32_t total_items_rec = 0;


        std::uint32_t displs[num_procs];
        std::uint32_t *items_size_total = nullptr;

        if (my_id == 0) {

            for(int i = 0; i< num_procs; ++i) {
                total_items_rec += recvcnts[i];

                if (i == 0) {
                    displs[i] = 0;
                }
                else {
                    displs[i] = displs[i-1] + recvcnts[i-1];
                }
            }

            items_size_total = (std::uint32_t *) malloc(total_items_rec * sizeof(std::uint32_t));

        }
        cout << "[JMABUIN] At proc " << my_id << " :: Before Gatherv :: " << endl;

        //https://blogs.cisco.com/performance/can-i-mpi_send-and-mpi_recv-with-a-count-larger-than-2-billion
        //MPI_Gatherv(items_size, num_items_to_send, MPI_UNSIGNED, items_size_total, recvcnts, displs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        MPI_Status status;

        if (my_id == 0) {

            memcpy(items_size_total, items_size, num_items_to_send * sizeof(std::uint32_t));

            cout << "[JMABUIN] At proc " << my_id << " :: After memcpy :: " << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        std::uint32_t *offset = items_size_total;

        for(int j = 1; j< num_procs; ++j) {

            if (my_id == 0) {

                /*
                 * MPI_Recv(
        void* data,
        int count,
        MPI_Datatype datatype,
        int source,
        int tag,
        MPI_Comm communicator,
        MPI_Status* status)

                 *
                 */
                offset += recvcnts[j-1];
                MPI_Recv(offset, (int)recvcnts[j], MPI_UINT32_T, j, j+num_procs, MPI_COMM_WORLD, &status);
                cout << "[JMABUIN] At proc " << my_id << " :: After Recv from rank :: " << j << ". Offset: " << offset << endl;
            }
            else if (my_id == j) {
                /*
                 * MPI_Send(
        void* data,
        int count,
        MPI_Datatype datatype,
        int destination,
        int tag,
        MPI_Comm communicator)

                 */
                MPI_Send(items_size, num_items_to_send, MPI_UINT32_T, 0, j+num_procs, MPI_COMM_WORLD);
                cout << "[JMABUIN] At proc " << my_id << " :: After send from rank :: " << j << endl;
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        cout << "[JMABUIN] At proc " << my_id << " :: After Gatherv :: " << endl;

        std::unordered_map<std::uint32_t, std::uint32_t> items_map;
        std::uint32_t *items_to_delete = nullptr;
        auto maxlpf = db.max_locations_per_feature() - 1;
        std::uint32_t total_to_delete;

        if (my_id == 0) {
            total_to_delete = 0;
            for(std::uint32_t i = 0; i< total_items_rec; i+=2) {

                if(items_map.find(items_size_total[i]) != items_map.end()) {
                    items_map[items_size_total[i]] = items_map[items_size_total[i]] + items_size_total[i+1];

                }
                else {
                    items_map.insert(std::pair<std::uint32_t, std::uint32_t>(items_size_total[i], items_size_total[i+1]));
                }
            }


            free(items_size_total);


            for(auto current_item = items_map.begin(); current_item != items_map.end(); ++current_item) {
                if (current_item->second > (std::uint32_t)maxlpf) {
                    total_to_delete++;
                }

            }

        }

        free(items_size);

        MPI_Bcast(&total_to_delete, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        cout << "[JMABUIN] At proc " << my_id << " :: After Bcast1 :: Received: "<< total_to_delete << endl;

        items_to_delete = (std::uint32_t *) malloc(total_to_delete * sizeof(std::uint32_t));

        if (my_id == 0) {
            unsigned k = 0;
            for(auto current_item = items_map.begin(); current_item != items_map.end(); ++current_item) {
                if (current_item->second > (std::uint32_t)maxlpf) {
                    items_to_delete[k] = current_item->first;
                    ++k;
                }

            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        cout << "[JMABUIN] At proc " << my_id << " :: Before Bcast2 :: " << endl;
        MPI_Bcast(items_to_delete, total_to_delete, MPI_UINT32_T, 0, MPI_COMM_WORLD);

        post_process_features_distributed(db, opt, items_to_delete, total_to_delete);

        MPI_Barrier(MPI_COMM_WORLD);
        cout << "[JMABUIN] Post processed features at proc " << my_id << " :: Now writing database :: " << endl;
    }


    string mpi_file_name;
    std::stringstream mpi_file_name_stream;
    mpi_file_name_stream <<  opt.dbfile << "_" << my_id;
    mpi_file_name = mpi_file_name_stream.str();

    if(notSilent) {
        cout << "Writing database to file '" << mpi_file_name << "' ... " << flush;
    }
    try {
        db.write(mpi_file_name);
        if(notSilent) cout << "done." << endl;
    }
    catch(const file_access_error&) {
        if(notSilent) cout << "FAIL" << endl;
        cerr << "Could not write database file!" << endl;
    }

    cout << "[JMABUIN] Database wrote at proc " << my_id << endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    time.stop();

    if(notSilent) {
        cout << "Total build time: " << time.seconds() << " s" << endl;
    }

    //prevents slow deallocation
    db.clear_without_deallocation();
}



/*************************************************************************//**
 *
 * @brief adds reference sequences to an existing database
 *
 *****************************************************************************/
void main_mode_build_modify(const args_parser& args, int my_id, int num_procs)
{
    auto dbfile = database_name(args);

    if(dbfile.empty()) {
        throw std::invalid_argument{"No database filename provided."};
    }

    cout << "Modify database " << dbfile << endl;

    auto db = make_database<database>(dbfile);

    auto opt = get_build_options(args, get_build_options_from_db(db));

    if(opt.infoLevel != info_level::silent && !opt.infiles.empty()) {
        cout << "Adding reference sequences to database..." << endl;
    }

    add_to_database(db, opt, my_id, num_procs);
}



/*************************************************************************//**
 *
 * @brief builds a database from reference input sequences
 *
 *****************************************************************************/
void main_mode_build(const args_parser& args, int my_id, int num_procs)
{
    auto opt = get_build_options(args);

    if((opt.infoLevel != info_level::silent) && (my_id == 0)) {
        cout << "Building new database from reference sequences." << endl;
    }

    if(opt.dbfile.empty()) {
        throw std::invalid_argument{"No database filename provided."};
    }

    if(opt.infiles.empty()) {
        throw std::invalid_argument{
            "Nothing to do - no reference sequences provided."};
    }

    //configure sketching scheme
    auto sketcher = database::sketcher{};
    sketcher.kmer_size(opt.kmerlen);
    sketcher.sketch_size(opt.sketchlen);

    auto db = database{sketcher};
    db.target_window_size(opt.winlen);
    db.target_window_stride(opt.winstride);

    add_to_database(db, opt, my_id, num_procs);
}


} // namespace mc
