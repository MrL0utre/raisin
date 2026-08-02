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
/**   NAME       : m.c                                     **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are definitions for the     **/
/**                multi-level schemes.                    **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 01 jan 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/************************************************************/




/** 
 * @file m.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 05 apr 2025
 * @brief  These lines are definitions for the  
 *        multi-level schemes.
 *      
 */
#ifndef M_H

#define M_H

/*
**  The defines and includes.
*/
#include <math.h>
#include "commons.h"
#include "rbh.h"
#include "crbh.h"
#include "ipart.h"
#include "a.h"
#include "dlist.h"
#include "pqueue.h"
#include "fm.h"

INT 
multilevel(Hypergraph  * h,
           Arch        * a,
           List       ** neighbors, 
           List       ** in_neighbors, 
           INT         * sort, 
           PART * partition, 
           INT           perform, 
           INT           tolerance,
           INT           k,
           INT           epsilon,
           INT           algo_cluster,
           int (*intialPart)(Hypergraph * h, List ** neighbors, List ** in_neighbors, PART * partition, INT k, INT epsilon),
           int (*refinement)(Hypergraph * h, Arch * a, List ** neighbors, List ** in_neighbors, INT * sort, PART * partition, INT perform, INT tolerance, INT k));

INT 
multilevel_cut(Hypergraph  * h,
               Arch        * a,
               List       ** neighbors, 
               List       ** in_neighbors, 
               INT         * sort, 
               PART * partition, 
               INT           perform, 
               INT           tolerance,
               INT           k,
               INT           epsilon,
               INT           algo_cluster);

INT 
multilevel_pmax(Hypergraph  * h,
                Arch        * a,
                List       ** neighbors, 
                List       ** in_neighbors, 
                INT         * sort, 
                PART * partition, 
                INT           perform, 
                INT           tolerance,
                INT           k,
                INT           epsilon,
                INT           algo_cluster);

#endif

