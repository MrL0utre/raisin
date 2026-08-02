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
