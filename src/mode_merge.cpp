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
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "cmdline_utility.h"
#include "filesys_utility.h"
#include "args_handling.h"
#include "config.h"
#include "query_options.h"
#include "taxonomy_io.h"
#include "classification.h"
#include "io_error.h"
#include "candidates.h"
#include "printing.h"


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


/*************************************************************************//**
 *
 * @brief merge parameters
 *
 *****************************************************************************/
struct merge_options
{
    taxonomy_options taxonomy;

    info_level infoLevel = info_level::moderate;

	query_options query;
};



/*************************************************************************//**
 *
 * @brief command line args -> merge options
 *
 *****************************************************************************/
merge_options
get_merge_options(const args_parser& args)
{
    merge_options opt;

    if(args.contains("silent")) {
        opt.infoLevel = info_level::silent;
    } else if(args.contains("verbose")) {
        opt.infoLevel = info_level::verbose;
    }

    opt.taxonomy = get_taxonomy_options(args);

    opt.query = get_query_options(args, {});

    //TODO: get hitsmin from file
    if(opt.query.classify.hitsMin == 0)
        opt.query.classify.hitsMin = 5;
    if(opt.query.classify.lowestRank < taxon_rank::Species)
        opt.query.classify.lowestRank = taxon_rank::Species;
    if(opt.query.output.lowestRank < taxon_rank::Species)
        opt.query.output.lowestRank = taxon_rank::Species;
    //TODO? multi threading
    if(opt.query.process.numThreads > 1)
        opt.query.process.numThreads = 1;

    return opt;
}



/*************************************************************************//**
 *
 * @brief forward stream past the next occurence of a specific character
 *
 *****************************************************************************/
inline void
forward(std::istream& is, char c)
{
    is.ignore(std::numeric_limits<std::streamsize>::max(), c);
}



/*************************************************************************//**
 *
 * @brief information about results file
 *
 *****************************************************************************/
struct results_source {
    std::string filename;
    std::streampos resultsBegin = 0;
    std::size_t numQueries = 0;
    int tophitsColumn = 0;
};



/*************************************************************************//**
 *
 * @brief parse classification result file
 *
 *****************************************************************************/
results_source
get_results_file_properties(const string& filename)
{
    results_source res;
    res.filename = filename;

    std::ifstream ifs(filename);
    if(!ifs.good()) throw io_error("could not open file " + filename);

    string line;

    //check classification rank
    while(ifs.good()) {
        getline(ifs, line);
        if(line[0] != '#') {
            throw io_format_error("classificaion ranks not found in file " + filename);
        }
        if(line.compare(0,16,"# Classification") == 0) {
            auto pos = line.find("sequence");
            if(pos != string::npos)
                throw io_format_error("cannot merge results on sequence level");
            break;
        }
    }

    //get layout
    while(ifs.good()) {
        getline(ifs, line);
        if(line[0] != '#') {
            throw io_format_error("TABLE_LAYOUT not found in file " + filename);
        }
        if(line.compare(0,15,"# TABLE_LAYOUT:") == 0) {
            std::stringstream lineStream(line.substr(15));
            string column;
            lineStream >> column;
            if(column != "query_id") {
                throw io_format_error("no query_id in file " + filename);
            }
            int col = 0;
            while(lineStream.good()) {
                forward(lineStream, '|');
                lineStream >> column;
                ++col;
                if(column == "top_hits") {
                    res.tophitsColumn = col;
                    break;
                }
            }
            break;
        }
    }
    if(res.tophitsColumn < 1)
        throw io_format_error("no top_hits in file " + filename);

    char lineBegin = ifs.peek();
    //skip comments
    while(ifs.good() && lineBegin == '#') {
        forward(ifs, '\n');
        lineBegin = ifs.peek();
    }
    //count query results
    res.resultsBegin = ifs.tellg();
    while(ifs.good()) {
        if(lineBegin != '#') ++res.numQueries;
        forward(ifs, '\n');
        lineBegin = ifs.peek();
    }

    return res;
}



/*************************************************************************//**
 *
 * @brief extract query headers and candidates from results file
 *
 *****************************************************************************/
void read_results(const results_source& res,
                  const database& db,
                  const candidate_generation_rules& rules,
                  vector<string>& queryHeaders,
                  vector<classification_candidates>& queryCandidates)
{
    std::ifstream ifs(res.filename);
    if(!ifs.good()) throw io_error("could not re-open file " + res.filename);
    ifs.seekg(res.resultsBegin);
    if(!ifs.good()) throw io_format_error("could not process file " + res.filename);

    queryCandidates.resize(res.numQueries);
    queryHeaders.resize(res.numQueries);

    char lineBegin = ifs.peek();
    while(ifs.good()) {
        if(lineBegin != '#') {
            size_t queryId;
            ifs >> queryId;
            if(queryId > 0) queryId--;

            forward(ifs, '|');
            // if we don't have a query header yet -> read from file
            if(queryHeaders[queryId].empty()) {
                string header;
                ifs >> header;
                if(!header.empty()) queryHeaders[queryId] = std::move(header);
            }

            // skip to tophits
            for(int i = 1; i < res.tophitsColumn; ++i) forward(ifs, '|');
            forward(ifs, '\t');

            // get tophits
            char nextChar = ifs.peek();
            while(nextChar != '\t') {
                taxon_id taxid;
                ifs >> taxid;
                forward(ifs, ':');
                match_candidate::count_type hits;
                ifs >> hits;

                const taxon* tax = db.taxon_with_id(taxid);
                if(tax) {
                    queryCandidates[queryId].insert(match_candidate{tax, hits}, db, rules);
                } else {
                    cerr << "taxid not found" << endl;
                }
                nextChar = ifs.get();
            }
            // should be at end of tophits
        }
        forward(ifs, '\n');
        lineBegin = ifs.peek();
    }
}



/*************************************************************************//**
 *
 * @brief merge classification result files
 *
 *****************************************************************************/
void merge_result_files(const vector<string>& infiles,
                        const database& db,
						const query_options& opt,
                        classification_results& results)
{
    vector<string> queryHeaders;
    vector<classification_candidates> queryCandidates;

    candidate_generation_rules rules;

    rules.mergeBelow    = opt.classify.lowestRank;
    if(opt.classify.maxNumCandidatesPerQuery > 0)
        rules.maxCandidates = opt.classify.maxNumCandidatesPerQuery;
    //else default to 2

    const auto& comment = opt.output.format.comment;

    cerr << " max canddidates: " << rules.maxCandidates << endl;
    cerr << " number of files: " << infiles.size() << endl;

    results.perReadOut << comment << "Merging " << infiles.size() << " files:\n";
    for(const auto& filename : infiles) {
        results.perReadOut << comment << filename << '\n';
    }

    for(size_t i = 0; i < infiles.size(); ++i) {
        show_progress_indicator(cerr, infiles.size() > 1 ? i/float(infiles.size()) : -1);

        read_results(get_results_file_properties(infiles[i]),
                     db, rules, queryHeaders, queryCandidates);
    }
    clear_current_line(cerr);

    map_candidates_to_targets(queryHeaders, queryCandidates, db, opt, results);
}



/*************************************************************************//**
 *
 * @brief runs merge on input files; sets output target streams
 *
 *****************************************************************************/
void process_result_files(const vector<string>& infiles,
                          const database& db, const query_options& opt,
                          const string& queryMappingsFilename,
                          const string& abundanceFilename)
{
    std::ostream* perReadOut   = &cout;
    std::ostream* perTargetOut = &cout;
    std::ostream* perTaxonOut  = &cout;
    std::ostream* status       = &cerr;

    std::ofstream mapFile;
    if(!queryMappingsFilename.empty()) {
        mapFile.open(queryMappingsFilename, std::ios::out);

        if(mapFile.good()) {
            cout << "Per-Read mappings will be written to file: " << queryMappingsFilename << endl;
            perReadOut = &mapFile;
            //default: auxiliary output same as mappings output
            perTargetOut = perReadOut;
            perTaxonOut  = perReadOut;
        }
        else {
            throw file_write_error{"Could not write to file " + queryMappingsFilename};
        }
    }

    std::ofstream abundanceFile;
    if(!abundanceFilename.empty()) {
        abundanceFile.open(abundanceFilename, std::ios::out);

        if(abundanceFile.good()) {
            cout << "Per-Taxon mappings will be written to file: " << abundanceFilename << endl;
            perTaxonOut = &abundanceFile;
        }
        else {
            throw file_write_error{"Could not write to file " + abundanceFilename};
        }
    }

    classification_results results {*perReadOut,*perTargetOut,*perTaxonOut,*status};

    if(opt.output.showQueryParams) {
        show_query_parameters(results.perReadOut, opt);
    }

    results.flush_all_streams();

    results.time.start();
    merge_result_files(infiles, db, opt, results);
    results.time.stop();

    clear_current_line(results.status);
    results.status.flush();

    if(opt.output.showSummary) show_summary(opt, results);

    results.flush_all_streams();
}



/*************************************************************************//**
 *
 * @brief runs merge on input files;
 *        handles output file split
 *
 *****************************************************************************/
void process_result_files(const vector<string>& infiles,
                         const database& db,
                         const query_options& opt)
{
    if(infiles.empty()) {
        cerr << "No input filenames provided!" << endl;
        return;
    }
    else {
        bool noneReadable = std::none_of(infiles.begin(), infiles.end(),
                           [](const auto& f) { return file_readable(f); });
        if(noneReadable) {
            throw std::runtime_error{
                        "None of the query sequence files could be opened"};
        }
    }

    //process all input files at once
    process_result_files(infiles, db, opt,
                         opt.output.queryMappingsFile,
                         opt.output.abundanceFile);
}



/*************************************************************************//**
 *
 * @brief merge classification result files
 *
 *****************************************************************************/
void main_mode_merge(const args_parser& args)
{
    const auto nargs = args.non_prefixed_count();
    if(nargs < 2) {
        cerr << "No result filenames provided. Abort." << endl;
    	return;
    }
    auto infiles = result_filenames(args);
    std::sort(infiles.begin(), infiles.end());

    merge_options opt = get_merge_options(args);

    database db;

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
    else {
        cout << "No taxonomy specified. Unable to perform merge." << endl;
        return;
    }

    //TODO parallelize?
    // ranks cache would need to be updated at species level before
    // running the result processing concurrently

    if(infiles.size() >= 2) {
        cerr << "Merging result files." << endl;

        process_result_files(infiles, db, opt.query);
    }
    else {
        cout << "Please provide at least two files to be merged!" << endl;
    }
}

} // namespace mc
