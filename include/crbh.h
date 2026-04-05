/*
 * Copyright (C) 2025 Commissariat à l'énergie atomique et aux énergies
 * alternatives (CEA) and Institut national de recherche en sciences et
 * technologies du numérique (INRIA)
 * Contributor(s): Julien Rodriguez <julien.ro34@gmail.com>
 *
 * Licensed under the GPL, Version 3 (the "License");
 * You may obtain a copy of the License at:
 * https://www.gnu.org/licenses/gpl-3.0.txt
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/


/** 
 * @file crbh.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  This file is part of the project.    
 *        It contains declarations of clustering functions
 *        for red-black hypergraphs partitioning software (raisin).
 *      
 */
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
