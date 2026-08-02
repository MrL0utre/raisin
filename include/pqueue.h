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
/**   NAME       : vector.h                                **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are priority queue          **/
/**                functions declarations.                 **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : xx xxx xxxx     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#ifndef PQUEUE_H
#include <stdlib.h>
#include "commons.h"

#define PQUEUE_H

void 
pqueue_sift_up(INT * queue, 
               INT * queue_w, 
               INT   i, 
               INT   n);

void 
pqueue_sift_down(INT * queue, 
                 INT * queue_w, 
                 INT   i, 
                 INT   n);

void 
pqueue_heapify(INT * queue, 
               INT * queue_w, 
               INT   n);

INT 
pqueue_dequeue(INT * queue, 
               INT * queue_w, 
               INT   n);

void 
pqueue_add_element(INT * queue, 
                   INT * queue_w, 
                   INT   i, 
                   INT   n);

#endif

