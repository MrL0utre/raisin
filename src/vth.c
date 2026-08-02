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
/**   NAME       : vth.c                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are declarations for the    **/
/**                vertex-to-hyperedge structure.         **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 05 apr 2025     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/************************************************************/

#include "vth.h"
#include <string.h>

VtH *
vth_build(const Hypergraph *h)
{
    INT nv = h->i_vertices;
    INT ne = h->i_hyperedges;

    VtH *vth = (VtH *)malloc(sizeof(VtH));
    MEM_ERROR(vth);
    vth->nv = nv;

    INT *deg = (INT *)calloc((size_t)nv, sizeof(INT));
    MEM_ERROR(deg);

    /* Pass 1 – count how many hyperedges touch each vertex */
    for (INT j = 0; j < ne; j++) {
        INT base = h->ti_idx_hyperedges[j];
        INT sz   = h->ti_hyperedges[base + 1];
        for (INT i = 0; i < sz; i++)
            deg[h->ti_hyperedges[base + 2 + i]]++;
    }

    /* Prefix sum → idx */
    vth->idx = (INT *)malloc((nv + 1) * sizeof(INT));
    MEM_ERROR(vth->idx);
    vth->idx[0] = 0;
    for (INT v = 0; v < nv; v++)
        vth->idx[v + 1] = vth->idx[v] + deg[v];

    vth->total = vth->idx[nv];
    vth->data  = (INT *)malloc(vth->total * sizeof(INT));
    MEM_ERROR(vth->data);

    /* Pass 2 – fill */
    memset(deg, 0, nv * sizeof(INT));
    for (INT j = 0; j < ne; j++) {
        INT base = h->ti_idx_hyperedges[j];
        INT sz   = h->ti_hyperedges[base + 1];
        for (INT i = 0; i < sz; i++) {
            INT v = h->ti_hyperedges[base + 2 + i];
            vth->data[vth->idx[v] + deg[v]++] = j;
        }
    }

    free(deg);
    return vth;
}

void
vth_free(VtH *vth)
{
    if (!vth) return;
    free(vth->data);
    free(vth->idx);
    free(vth);
}
