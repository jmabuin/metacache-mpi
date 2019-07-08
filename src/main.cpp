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
#include <stdexcept>
#include <string>

//extern "C"
//{
#include <mpi.h>
//}


#include "args_parser.h"
#include "modes.h"


/*************************************************************************//**
 *
 * @brief  select mode or display quick help
 *
 *****************************************************************************/
int main(int argc, char* argv[])
{
    using namespace mc;

    int         my_id;
    int         num_procs;

    MPI_Init( &argc, &argv );
    // Find out number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Find out process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);


    //std::cerr << "\nMy_id: " << my_id << "\nNum_procs: " << num_procs << std::endl;
    args_parser args{argc, argv};
    
    bool nomode = false;
    auto modestr = std::string{};
    if(args.non_prefixed_count() > 0) {
        modestr = args.non_prefixed(0);
    }

    try {
        if(modestr == "help") {
            main_mode_help(args);
        }
        else if(modestr == "build") {
            main_mode_build(args, my_id, num_procs);
        }
        else if(modestr == "modify") {
            main_mode_build_modify(args, my_id, num_procs);
        }
        else if(modestr == "query") {
            main_mode_query(args);
        }
        else if(modestr == "info") {
            main_mode_info(args);
        }
        else if(modestr == "annotate") {
            main_mode_annotate(args);
        }
        else if(modestr == "merge") {
            main_mode_merge(args);
        }
        else {
            nomode = true;
        }
    }
    catch(std::runtime_error& e) {
        std::cerr << "\nABORT: " << e.what() << "!" << std::endl;
    }
    catch(std::invalid_argument& e) {
        std::cerr << "\nABORT: " << e.what() << "!" << std::endl;
        nomode = true;
    }

    if(nomode) {
       main_mode_help(args);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;

}
