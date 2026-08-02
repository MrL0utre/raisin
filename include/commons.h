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


/** 
 * @file commons.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 05 apr 2025
 * @brief  This file is part of the project.    
 *        It contains declarations of commons macros and functions
 *        for red-black hypergraphs partitioning software (raisin).
 *      
 */
#ifndef COMMONS_H

#define COMMONS_H
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

#define INT int


#ifndef RAISIN_PART_INT
typedef uint8_t PART;
#define RAISIN_PART_MAX ((INT)UINT8_MAX)
#else
typedef int     PART;
#define RAISIN_PART_MAX INT_MAX
#endif

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MAX3(X, Y, Z) (MAX(X, MAX(Y, Z)))
#define RAISIN_UNUSED(x) ((void)(x))
#define MEM_ERROR(x) do { if (!(x)) { fprintf(stderr, __FILE__ ":%d:unable to allocate buffer `" #x "'\n", __LINE__); exit(EXIT_FAILURE); } } while (0)
#define FATAL(cond, s) do { if (cond) { fprintf(stderr, __FILE__ ":%d:" s "\n", __LINE__); exit(EXIT_FAILURE); } } while (0)
#define BUFSIZE 16777216

#endif
