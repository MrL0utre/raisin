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
/**   NAME       : crbh.h                                  **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are clustering functions    **/
/**                declarations for red-black hypergraph.  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 13 oct 2023     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#ifndef CRBH_H
#define CRBH_H

#define SPEED_DELAY_CLUSTERING 1

#include "commons.h"
#include "rbh.h"
#include "a.h"
#include "matrix.h"
#include "pqueue.h"

/**
 * @brief Compute a reduced hypergraph (h2) according to clustering of 
 * an input red-black hypergraph (h1).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
best_phi_clustering(Hypergraph *h1, 
                    Hypergraph *h2,
                    Arch        *a,
                    List       **neighbors, 
                    List       **in_neighbors,
                    INT         *sort, 
                    INT          umap[],
                    INT          cluster_size, 
                    INT          epsilon, 
                    INT          relax);

int 
coarsening(Hypergraph   *h1, 
           Hypergraph   *h2,
           List        **neighbors, 
           List        **in_neighbors,
           unsigned int *map, 
           unsigned int  k);

/**
 * @brief Compute a matching (map array) of a red-black hypergraph (h1) vertices according 
 * to criticality weighting and construct a new contracted 
 * red-black hypergraph (h2).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 * @param D                Delay factor.
 *
 * @return Return a flag.
 */        
int 
heavy_edge_matching(Hypergraph *h1, 
                    Hypergraph *h2,
                    List      **neighbors, 
                    List      **in_neighbors,
                    INT        *map, 
                    INT         k, 
                    INT         epsilon);

/**
 * @brief Sort function (quicksort) of an array of integers.
 *
 * @param array       Array of integers.
 * @param keys        Arrays of integers(keys)
 * @param low         Lowerbound.
 * @param high        Upperbound.
 *
 * @return Partition position.
 */
void 
quick_sort(int array[], 
           int keys[], 
           int low, 
           int high);

/**
 * @brief Function to find the partition position for quicksort.
 *
 * @param array       Array of integers.
 * @param keys        Arrays of integers(keys)
 * @param low         Lowerbound.
 * @param high        Upperbound.
 *
 * @return Partition position.
 */
int 
partition(int array[], 
          int keys[], 
          int low, 
          int high);

/**
 * @brief Function to swap elements.
 *
 * @param a            Pointer.
 * @param b            Pointer.
 *
 * @return void.
 */
void 
swap(int *a, 
     int *b);

/**
 * @brief Compute a reduced hypergraph (h2) according to clustering of 
 * an input red-black hypergraph (h1).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
compute_hypergraph_kclustering(Hypergraph *h1, 
                               Hypergraph *h2,
                               List **neighbors, 
                               List **in_neighbors,
                               INT *sort, 
                               INT *cluster, 
                               INT nvp);

/**
 * @brief Compute a reduced hypergraph (h2) according to matching of 
 * an input red-black hypergraph (h1).
 * 
 * @attention This algorithm only work for matching (cluster of size 2 or 1)  
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
compute_hypergraph_clustering(Hypergraph *h1, 
                              Hypergraph *h2,
                              INT *cluster, 
                              INT nvp);

/**
 * @brief Compute the critical path length of a clustering partition
 *        in a red-black hypergraph.
 *
 * @param h                Input hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 * @param D                Delay factor.
 *
 * @return The critical path length (maximum accumulated delay).
 */
int 
compute_clustering_criticality(Hypergraph *h,
                               List      **neighbors, 
                               List      **in_neighbors,
                               INT        *sort, 
                               INT        *partition, 
                               INT         D);

int 
compute_subhypergraph(Hypergraph * h1, 
                      Hypergraph * h2, 
                      INT        * cluster, 
                      INT          nvp);


#endif /* CRBH_H */
