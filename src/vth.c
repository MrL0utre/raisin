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

    /* Persistent buffer to avoid use-after-free when vth_build is called in loops */
    static INT *deg_p = NULL;
    static INT  deg_cap = 0;
    if (nv > deg_cap) {
        deg_p = (INT *)realloc(deg_p, nv * sizeof(INT));
        deg_cap = nv;
    }
    memset(deg_p, 0, nv * sizeof(INT));
    INT *deg = deg_p;

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
