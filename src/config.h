/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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

/*************************************************************************//**
 *
 * @file contains all important type aliases
 *       controlling database space consumption and basic sketching scheme
 *
 *****************************************************************************/

#ifndef MC_CONFIG_H_
#define MC_CONFIG_H_


#include "sketch_database.h"
#include "hash_dna.h"


namespace mc {


/********************************************************************
 * @brief limits max. k (= number of letters) per k-mer
 */
#ifdef MC_KMER_TYPE
    using kmer_type = MC_KMER_TYPE ;
#else
    using kmer_type = std::uint32_t;
#endif


/********************************************************************
 * @brief limits number of targets (reference sequences) per database
 */
#ifdef MC_TARGET_ID_TYPE
    using target_id = MC_TARGET_ID_TYPE ;
#else
    using target_id = std::uint16_t;
#endif


/********************************************************************
 * @brief limits max. number of windows per target (reference sequence)
 */
#ifdef MC_WINDOW_ID_TYPE
    using window_id = MC_WINDOW_ID_TYPE ;
#else
    using window_id = std::uint32_t;
#endif


/********************************************************************
 * @brief limits max. number of locations per feature
 */
#ifdef MC_LOCATION_LIST_SIZE_TYPE
    using loclist_size_t = MC_LOCATION_LIST_SIZE_TYPE ;
#else
    using loclist_size_t = std::uint8_t;
#endif


/**************************************************************************
 * @brief nucleotide sequence storage type
 */
using sequence = std::string;
using query_id = std::uint_least64_t;


/**************************************************************************
 * @brief controls how nucleotide sequences are transformed into 'features'
 *        (called "h_1" in the paper)
 */
using sketching_hash = same_size_hash<kmer_type>;

//using sketcher = single_function_min_hasher<kmer_type,sketching_hash>;
using sketcher = single_function_unique_min_hasher<kmer_type,sketching_hash>;


/**************************************************************************
 * @brief hash function for database hash multi-map (called "h_2" in the paper)
 *        note: std::hash<SomeUnsignedIntegerType>
 *              is mostly implemented as the identity function
 */
//using feature_hash = std::hash<typename sketcher::feature_type>;
using feature_hash = same_size_hash<typename sketcher::feature_type>;


/**************************************************************************
 * @brief stores a multimap that maps features to locations
 *        within targets (= reference sequences)
 */
using database = sketch_database<sequence,sketcher,feature_hash,
                                 target_id,window_id,loclist_size_t>;

/** @brief pull some types from database into global namespace */
using taxon                  = database::taxon;
using taxon_rank             = database::taxon_rank;
using taxon_id               = database::taxon_id;
using ranked_lineage         = database::ranked_lineage;
using match_locations        = database::match_locations;
using match_target_locations = database::match_target_locations;


/**************************************************************************
 * @brief controls how a classification is derived from a location hit list;
 *        default is a top 2 voting scheme;
 *        forward declarations (breaks cycle "config.h" <-> "candidates.h")
 */

class distinct_matches_in_contiguous_window_ranges;
class best_distinct_matches_in_contiguous_window_ranges;


using classification_candidates =
    best_distinct_matches_in_contiguous_window_ranges;
//    distinct_matches_in_contiguous_window_ranges;


} // namespace mc


#endif
