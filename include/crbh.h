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
#define __SPEED_DELAY_CLUSTERING__ 1

#include "commons.h"
#include "rbh.h"
#include "a.h"
#include "matrix.h"
#include"pqueue.h"

int bestPhiClustering(Hypergraph * h1, Hypergraph * h2, Arch *a, List ** neighbors, List ** in_neighbors, INT * sort, INT umap[], INT cluster_size, INT epsilon, INT relax);
int coarsening(Hypergraph * h1, Hypergraph * h2, List ** neighbors, List ** in_neighbors, unsigned int * map, unsigned int k);
int heavyEdgeMatching(Hypergraph * h1, Hypergraph * h2, List ** neighbors, List ** in_neighbors, INT * map, INT k, INT epsilon);
void quickSort(int array[], int keys[], int low, int high);
int partition(int array[], int keys[], int low, int high);
void swap(int *a, int *b);
int compute_hypergraph_kclustering(Hypergraph * h1, Hypergraph * h2, List ** neighbors, List ** in_neighbors, INT * sort, INT * cluster, INT nvp);
int compute_hypergraph_clustering(Hypergraph * h1, Hypergraph * h2, INT * cluster, INT nvp);
int compute_clustering_criticality(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * partition, INT D);


#endif
