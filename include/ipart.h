/*
 * Copyright (C) 2025 Commissariat à l'énergie atomique et aux énergies
 * alternatives (CEA) and Institut national de recherche en sciences et
 * technologies du numérique (INRIA)
 * Contributor(s): Julien Rodriguez <julien.ro34@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
*/


/************************************************************/
/**                                                        **/
/**   NAME       : ipart.h                                 **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are initial partitioning    **/
/**                functions declarations for red-black    **/
/**                hypergraph.                             **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#ifndef IPART_H

#define IPART_H

#include "commons.h"
#include "rbh.h"
#include "crbh.h"
#include "a.h"
#include "matrix.h"
#include "pqueue.h"
#include "math.h"
#include "uf.h"

int
derived_breadth_first_search_multilevel(Hypergraph * h, 
                                       List      ** neighbors, 
                                       List      ** in_neighbors,
                                       PART * partition, 
                                       INT          k, 
                                       INT epsilon);

int 
derived_breadth_first_search(Hypergraph * h, 
                             List      ** neighbors, 
                             List      ** in_neighbors, 
                             PART  * partition, 
                             INT          k, 
                             INT          epsilon);

int 
derived_depth_first_search(Hypergraph * h, 
                           List      ** neighbors, 
                           List      ** in_neighbors, 
                           PART  * partition, 
                           INT          k, 
                           INT          epsilon);

int 
derived_depth_first_searchMultilevel(Hypergraph * h, 
                                     List      ** neighbors, 
                                     List      ** in_neighbors, 
                                     PART       * partition, 
                                     INT          k, 
                                     INT          epsilon);

int 
critical_connected_component_partitioning(Hypergraph * h,
                                          List      ** neighbors, 
                                          List      ** in_neighbors, 
                                          INT        * umap, 
                                          INT          k, 
                                          INT          epsilon, 
                                          INT          bound);

int 
critical_connected_component_partitioning_multilevel(Hypergraph * h, 
                                                     List ** neighbors, 
                                                     List ** in_neighbors, 
                                                     INT  * umap,
                                                     INT k, 
                                                     INT epsilon, 
                                                     INT bound);

#endif

