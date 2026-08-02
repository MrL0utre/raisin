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
/**   NAME       : adj.h                                   **/
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

/**
 * @file adj.h
 * @brief Flat CSR adjacency structure replacing the linked-list neighbour
 *        representation used throughout rbh.c and fm.c.
 *
 * Compatibility
 * -------------
 * Functions that still take List** are kept unchanged; new code should
 * prefer AdjCSR.  A thin adapter adj_to_list() is provided for the
 * transition period.
 */

#ifndef ADJ_H
#define ADJ_H

#include "commons.h"
#include "rbh.h"

/**
 * Compressed Sparse Row adjacency.
 *
 *   Neighbours of vertex v : data[ idx[v] .. idx[v+1] )
 */
typedef struct {
    INT *data;   /* concatenated neighbour lists              */
    INT *idx;    /* idx[v] = start of v's list in data[]     */
    INT  nv;     /* number of vertices                        */
    INT  total;  /* total number of (src,dst) pairs = |data| */
} AdjCSR;

/**
 * Build outgoing CSR from a hypergraph.
 * Multi-edges are collapsed (each (u,v) pair appears at most once).
 * O(pins) time, O(nv + pins) space.
 */
AdjCSR *adj_csr_build_out(Hypergraph *h);

/**
 * Build incoming CSR from an already-built outgoing CSR.
 * O(nv + pins) time.
 */
AdjCSR *adj_csr_build_in(const AdjCSR *out, INT nv);

/** Free an AdjCSR (also frees the struct itself). */
void    adj_csr_free(AdjCSR *adj);

/**
 * Degree of vertex v.  Inline for use in hot loops.
 */
static inline INT
adj_deg(const AdjCSR *adj, INT v)
{
    return adj->idx[v + 1] - adj->idx[v];
}

/**
 * Pointer to first neighbour of v (for range-based iteration).
 *   for (INT *nb = adj_begin(adj,v); nb != adj_end(adj,v); nb++)
 *       ... use *nb ...
 */
static inline const INT *
adj_begin(const AdjCSR *adj, INT v)
{
    return adj->data + adj->idx[v];
}

static inline const INT *
adj_end(const AdjCSR *adj, INT v)
{
    return adj->data + adj->idx[v + 1];
}

#endif /* ADJ_H */
