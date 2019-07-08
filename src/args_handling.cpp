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

#include <stdexcept>
#include <algorithm>

#include "filesys_utility.h"
#include "args_handling.h"


namespace mc {


//-------------------------------------------------------------------
std::string
database_name(const args_parser& args)
{
    //first 2 non-prefixed args are mode and database name
    if(args.non_prefixed_count() > 1) {
        auto filename = args.non_prefixed(1);
        if(filename.find(".db") == std::string::npos) {
            filename += ".db";
        }
        return filename;
    }
    return "";
}



//-------------------------------------------------------------------
std::vector<std::string>
input_filenames(const args_parser& args, const std::size_t skip = 2)
{
    //first <skip> non-prefixed args are ignored
    //all other non-prefixed args that are not
    //preceded by options ("-xyz") should be file or folder names
    auto n = args.non_prefixed_count();

    auto files = std::vector<std::string>{};
    files.reserve(n-1);
    for(std::size_t i = skip; i < n; ++i) {
        if(!args.is_preceded_by_prefixed_arg(i) && //no option
           !args.is_preceded_by_prefixed_arg(i-1)) //no option value
        {
            auto name = args.non_prefixed(i);

            auto fnames = files_in_directory(name);

            if(!fnames.empty()) {
                files.insert(files.end(), std::make_move_iterator(fnames.begin()),
                                          std::make_move_iterator(fnames.end()));
            } else {
                files.push_back(name);
            }
        }
    }

    return files;
}



//-------------------------------------------------------------------
std::vector<std::string>
sequence_filenames(const args_parser& args)
{
    //first 2 non-prefixed args are mode and database name
    //all other non-prefixed args that are not
    //preceded by options ("-xyz") should be sequence file or folder names
    return input_filenames(args, 2);
}



//-------------------------------------------------------------------
std::vector<std::string>
result_filenames(const args_parser& args)
{
    //first non-prefixed arg is mode
    //all other non-prefixed args that are not
    //preceded by options ("-xyz") should be result file or folder names
    return input_filenames(args, 1);
}



//-------------------------------------------------------------------
taxonomy_options
get_taxonomy_options(const args_parser& args)
{
    taxonomy_options opt;

    opt.path = args.get<std::string>("taxonomy", std::string(""));
    if(!opt.path.empty() &&
        opt.path.back() != '/') opt.path += '/';

    opt.nodesFile = opt.path + "nodes.dmp";
    opt.namesFile = opt.path + "names.dmp";
    opt.mergeFile = opt.path + "merged.dmp";

    opt.mappingPreFiles.push_back("assembly_summary.txt");

    //manually added accession to taxon map file names
    auto postm = args.get<std::string>({"taxpostmap", "taxonomy-postmap"}, "");

    if(!postm.empty()) {
        opt.mappingPostFiles.push_back(postm);
    }

    //default NCBI accession to taxon map file names
    opt.mappingPostFiles.push_back(opt.path + "nucl_gb.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_wgs.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_est.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_gss.accession2taxid");

    //find additional maps by file extension ".accession2taxid"
    for(const auto f : files_in_directory(opt.path)) {
        if(f.find(".accession2taxid") != std::string::npos) {
            if(std::find(opt.mappingPostFiles.begin(),
                         opt.mappingPostFiles.end(), f)
               == opt.mappingPostFiles.end())
            {
                opt.mappingPostFiles.push_back(f);
            }
        }
    }

    return opt;
}



} // namespace mc
