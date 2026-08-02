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
/**   NAME       : kfm.h                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are derived local search    **/
/**                functions declaration based on FM       **/
/**                (Fiduccia & Mattheyses).                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

/** 
 * @file fm.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  FM based refinement algorithms and functions  
 *        declarations for red-black hypergraph.
 *      
 */

#ifndef KFM_H

#define KFM_H
#include <math.h>
#include "commons.h"
#include "rbh.h"
#include "a.h"
#include "dlist.h"
#include "pqueue.h"

INT 
compute_partition_criticality(Hypergraph * h, 
                              Arch * a, 
                              List ** neighbors, 
                              List ** in_neighbors, 
                              INT * sort, 
                              PART * partition);

/* compute the local critical path through the vertex (vertex) according to a partition and a target topology (partition) */
INT 
compute_local_partition_criticality(Hypergraph * h, 
                                    Arch * a, 
                                    List ** neighbors, 
                                    List ** in_neighbors, 
                                    INT * sort, 
                                    PART * partition,
                                    INT vertex);

INT 
dkfm_fast(Hypergraph * h,  
          Arch * a, 
          List ** neighbors, 
          List ** in_neighbors, 
          INT * sort, 
          PART * partition, 
          INT perform, 
          INT tolerance,
          INT k);

INT 
dkfm(Hypergraph * h,  
     Arch *, 
     List ** neighbors, 
     List ** in_neighbors, 
     INT * sort, 
     PART * partition, 
     INT perform, 
     INT tolerance, 
     INT k);

INT 
kfm(Hypergraph * h,  
    Arch * a, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * sort, 
    PART * partition, 
    INT perform, 
    INT tolerance,
    INT k);

INT 
compute_partition_cut(Hypergraph * h, 
                      Arch * a, 
                      List ** neighbors, 
                      List ** in_neighbors, 
                      INT * sort, 
                      PART * partition,
                      INT * lambda,
                      INT k);

INT
compute_partition_criticality_with_buf(Hypergraph *h,
                                       Arch       *a,
                                       List      **out_neighbors,
                                       List      **in_neighbors,
                                       INT        *sort,
                                       PART * partition,
                                       INT        *delays,
                                       bool       *is_red,
                                       bool       *flag);
#endif
