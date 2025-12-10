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
/**   NAME       : kfm.h                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are derived local search    **/
/**                functions declaration based on FM       **/
/**                (Fiduccia & Mattheyses).                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : xx xxx xxxx     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#ifndef KFM_H

#define KFM_H
#include <math.h>
#include "commons.h"
#include "rbh.h"
#include "a.h"
#include "dlist.h"
#include "pqueue.h"

INT compute_partition_criticality(
Hypergraph * h, 
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition);

/* compute the local critical path through the vertex (vertex) according to a partition and a target topology (partition) */
INT compute_local_partition_criticality(
Hypergraph * h, 
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition,
INT vertex);

INT dkfmFast(
Hypergraph * h,  
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance,
INT k);

INT dkfm(
Hypergraph * h,  
Arch *, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance, 
INT k);

INT kfm(
Hypergraph * h,  
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance,
INT k);

INT compute_partition_cut(
Hypergraph * h, 
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition,
INT * lambda,
INT k);
#endif
