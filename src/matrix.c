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
/**   NAME       : matrix.c                                **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are matrix                  **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#include "matrix.h"


void new_matrix(Matrix * mat, INT m, INT n) {
  
  mat->v = (Vector*)malloc(sizeof(Vector)*m);
  MEM_ERROR(mat->v);
  for (INT i=0; i < m; ++i) {
    mat->v[i].v = (INT*)malloc(sizeof(INT)*n);
    MEM_ERROR(mat->v[i].v);
    mat->v[i].size = n;
  }

  mat->m = m;
  mat->n = n;
  return;
}

void delete_matrix(Matrix * mat) {

  for (INT i=0; i < mat->m; ++i) {

     free(mat->v[i].v);

  }
  free(mat->v);
  free(mat);
  return;
}
