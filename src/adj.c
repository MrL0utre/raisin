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
/**   NAME       : adj.c                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are adjacency neighbour     **/
/**                storing functions declarations for      **/
/**                red-black hypergraph.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "adj.h"
#include <string.h>

/* ------------------------------------------------------------------ */
/* Outgoing adjacency                                                  */
/* ------------------------------------------------------------------ */

AdjCSR *
adj_csr_build_out(Hypergraph *h)
{
    INT nv      = h->i_vertices;
    INT ne      = h->i_hyperedges;
    INT *e      = h->ti_hyperedges;
    INT *idx_e  = h->ti_idx_hyperedges;

    AdjCSR *adj  = (AdjCSR *)malloc(sizeof(AdjCSR));
    MEM_ERROR(adj);
    adj->nv  = nv;

    INT *deg  = (INT  *)calloc(nv + 1, sizeof(INT));
    MEM_ERROR(deg);

    /* seen[v] = id of current source being processed, or -1 */
    INT *seen = (INT *)malloc(sizeof(INT) * nv);
    MEM_ERROR(seen);
    memset(seen, -1, sizeof(INT) * nv);

    /* Temporary stack to reset seen[] entries per source vertex */
    INT *touched     = (INT *)malloc(sizeof(INT) * nv);
    MEM_ERROR(touched);

    /* ---------- Pass 1: count unique out-neighbours per source ---------- */
    for (INT j = 0; j < ne; j++) {
        INT idx_j  = idx_e[j];
        INT size_j = e[idx_j + 1];
        INT u      = e[idx_j + 2];          /* source vertex */

        INT n_touched = 0;

        for (INT i = 0; i < size_j - 1; i++) {
            INT v = e[idx_j + 3 + i];       /* sink vertex */
            if (seen[v] != u) {
                seen[v]          = u;
                touched[n_touched++] = v;
                deg[u]++;
            }
        }

        /* reset only the entries we touched */
        for (INT t = 0; t < n_touched; t++)
            seen[touched[t]] = -1;
    }

    /* ---------- Prefix sum → idx[] ---------- */
    adj->idx = (INT *)malloc(sizeof(INT) * (nv + 1));
    MEM_ERROR(adj->idx);

    adj->idx[0] = 0;
    for (INT v = 0; v < nv; v++)
        adj->idx[v + 1] = adj->idx[v] + deg[v];

    INT total   = adj->idx[nv];
    adj->total  = total;
    adj->data   = (INT *)malloc(sizeof(INT) * total);
    MEM_ERROR(adj->data);

    /* Reuse deg[] as a fill pointer (offset from idx[v]) */
    memset(deg, 0, sizeof(INT) * nv);
    memset(seen, -1, sizeof(INT) * nv);

    /* ---------- Pass 2: fill data[] ---------- */
    for (INT j = 0; j < ne; j++) {
        INT idx_j  = idx_e[j];
        INT size_j = e[idx_j + 1];
        INT u      = e[idx_j + 2];

        INT n_touched = 0;

        for (INT i = 0; i < size_j - 1; i++) {
            INT v = e[idx_j + 3 + i];
            if (seen[v] != u) {
                seen[v]              = u;
                touched[n_touched++] = v;
                adj->data[adj->idx[u] + deg[u]++] = v;
            }
        }

        for (INT t = 0; t < n_touched; t++)
            seen[touched[t]] = -1;
    }

    free(deg);
    free(seen);
    free(touched);
    return adj;
}

/* ------------------------------------------------------------------ */
/* Incoming adjacency (built from outgoing)                            */
/* ------------------------------------------------------------------ */

AdjCSR *
adj_csr_build_in(const AdjCSR *out, INT nv)
{
    AdjCSR *in = (AdjCSR *)malloc(sizeof(AdjCSR));
    MEM_ERROR(in);
    in->nv    = nv;
    in->total = out->total;   /* same number of edges */

    INT *deg = (INT *)calloc(nv + 1, sizeof(INT));
    MEM_ERROR(deg);

    /* Count in-degree of each vertex */
    for (INT u = 0; u < nv; u++)
        for (INT i = out->idx[u]; i < out->idx[u + 1]; i++)
            deg[out->data[i]]++;

    in->idx = (INT *)malloc(sizeof(INT) * (nv + 1));
    MEM_ERROR(in->idx);

    in->idx[0] = 0;
    for (INT v = 0; v < nv; v++)
        in->idx[v + 1] = in->idx[v] + deg[v];

    in->data = (INT *)malloc(sizeof(INT) * in->total);
    MEM_ERROR(in->data);

    memset(deg, 0, sizeof(INT) * nv);

    for (INT u = 0; u < nv; u++) {
        for (INT i = out->idx[u]; i < out->idx[u + 1]; i++) {
            INT v = out->data[i];
            in->data[in->idx[v] + deg[v]++] = u;
        }
    }

    free(deg);
    return in;
}

/* ------------------------------------------------------------------ */
/* Free                                                                */
/* ------------------------------------------------------------------ */

void
adj_csr_free(AdjCSR *adj)
{
    if (!adj) return;
    free(adj->data);
    free(adj->idx);
    free(adj);
}
