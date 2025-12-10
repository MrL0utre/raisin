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
/**   NAME       : ipart.h                                 **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are initial partitioning    **/
/**                functions declarations for red-black    **/
/**                hypergraph.                             **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 13 oct 2023     **/
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

derivedBreadthFirstSearchMultilevel(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon);
int derivedBreadthFirstSearch(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon);
int derivedDepthFirstSearch(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon);
int derivedDepthFirstSearchMultilevel(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon);
int criticalConnectedComponentPartitioning(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon, INT bound);
criticalConnectedComponentPartitioningMultilevel(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon, INT bound);

#endif
