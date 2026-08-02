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
/**   NAME       : io.h                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are read and write function **/
/**                declarations.                           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 04 feb 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


/** 
 * @file io.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 05 apr 2025
 * @brief  These lines are read and write function  
 *        declarations.
 *      
 */
#ifndef IO_H

#define IO_H
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"

int 
read_line(FILE * file, 
          char * buffer, 
          int size);

int 
load_partition(INT i_vertices, 
               PART * partition, 
               const char * s_path);

int 
write_partition(INT i_vertices,
                PART * partition, 
                const char * file_path);

#endif
