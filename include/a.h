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
/**   NAME       : a.h                                     **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are architecture functions  **/
/**                declaration.                            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : xx xxx xxxx     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#ifndef A_H

#define A_H
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include<stdbool.h>
#include "commons.h"
#include "io.h"

/*
**  The type and structure definitions.
*/


/*+ The Architecture class type. +*/


typedef struct arch{
  char *            s_arch_name;          /*+ Architecture name.              +*/
  INT  *            ti_capacity;          /*+ Capacities for each element.    +*/
  INT  *            ti_delay;             /*+ Time crossing between elements. +*/
  INT               i_m;                  /*+ Dimension of architecture.      +*/
  INT               i_n;                  /*+ Dimension of architecture.      +*/
} Arch;

/*
**  The function prototypes.
*/

int               arch_init   (Arch * this, INT m, INT n);
int               arch_free   (Arch * this);
int               arch_load   (Arch * this, const char * s_path, bool verbose);

/*
**  The macro definitions.
*/

#endif
