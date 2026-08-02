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
/**   NAME       : vth.h                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   :  Vertex-to-Hyperedge CSR index (VtH).   **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/
#ifndef VTH_H
#define VTH_H

#include "commons.h"
#include "rbh.h"

typedef struct vth_s {
    INT *data;   /* hyperedge indices, packed                     */
    INT *idx;    /* idx[v]..idx[v+1] = slice of data for vertex v */
    INT  nv;
    INT  total;  /* == i_pins                                     */
} VtH;

VtH *vth_build(const Hypergraph *h);
void vth_free  (VtH *vth);

static inline INT        vth_deg  (const VtH *vth, INT v) { return vth->idx[v+1] - vth->idx[v]; }
static inline const INT *vth_begin(const VtH *vth, INT v) { return vth->data + vth->idx[v]; }
static inline const INT *vth_end  (const VtH *vth, INT v) { return vth->data + vth->idx[v+1]; }

#endif /* VTH_H */
