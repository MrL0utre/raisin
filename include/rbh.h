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
/**   NAME       : rbh.h                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are declarations for the    **/
/**                red-black hypergraph structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 04 feb 2022     **/
/**                                 to   : 13 oct 2023     **/
/**                                                        **/
/************************************************************/


#ifndef RBH_H
#define RBH_H
#include "commons.h"
#include "io.h"
#include "matrix.h"
#include "list.h"

/*
**  The type and structure definitions.
*/


/* Forward declaration — vth.h includes rbh.h so we break the cycle */
typedef struct vth_s VtH;


/*+ The red-black hypergraph class type. +*/

typedef struct hypergraph{
  char *            s_rbh_name;                /*+ Hypergraph name.   +*/
  INT               i_vertices;                /*+ Vertex number.     +*/
  INT               i_hyperedges;              /*+ Hyperedges number. +*/
  INT               i_reds;                    /*+ Red vertex number. +*/
  INT               i_pins;                    /*+ Pin number.        +*/
  INT               i_weights;                 /*+ Weight dimension.  +*/

  INT *             ti_hyperedges;             /*+ Hyperedges.        +*/
  INT *             ti_idx_hyperedges;         /*+ Index hyperedges.  +*/
  INT *             ti_delays;                 /*+ Delays.            +*/
  INT *             ti_criticalities_right;    /*+ Criticalities.     +*/
  INT *             ti_criticalities_left;     /*+ Criticalities.     +*/
  INT *             ti_reds;                   /*+ Red vertex.        +*/
  INT *             ti_weights;                /*+ Vertex weight.     +*/

  bool *            is_red;                    /*+ Static red/black mask.   +*/
  VtH  *            vth;                       /*+ Vertex→Hyperedge index.  +*/

  INT               (* rbh_init)  (struct hypergraph *); /*+ Hypergraph init function. +*/
  INT               (* rbh_free)  (struct hypergraph *); /*+ Hypergraph free function. +*/
  INT               (* rbh_load)  (struct hypergraph *, const char *, INT, bool); /*+ Hypergraph loading function. +*/
  INT               (* rbh_save)  (struct hypergraph *, const char *, INT, bool); /*+ Hypergraph saving function. +*/

} Hypergraph;

/*
**  The structure for function parameters.
*/
typedef struct {
    Hypergraph *      h;                  /*+ Hypergraph. +*/
    char *            s_path;             /*+ File path. +*/
    INT               i_baseval;          /*+ Vertex indexation. +*/
    bool              b_verbose;          /*+ Verbose. +*/
} rbhLoad_args;

int               var_rbh_load  (rbhLoad_args);

/*
**  The function prototypes.
*/
int               rbh_init           (Hypergraph * this);
int               rbh_free           (Hypergraph * this);
int               rbh_build_precomputed(Hypergraph * h);
int               rbh_validate       (const Hypergraph * h);
int               rbh_load_base       (Hypergraph * this, const char * const s_path, INT i_baseval, bool b_verbose);
int               rbh_save           (Hypergraph * this, const char * const s_path, INT i_baseval, bool b_verbose);
int               compute_neighbors  (Hypergraph * this, Matrix * neighbors);
int               compute_in_neighbors(Hypergraph * this, Matrix * neighbors, Matrix * in_neighbors);
int               compute_list_neighbors  (Hypergraph * this, List ** neighbors);
int               compute_list_in_neighbors(Hypergraph * this, List ** neighbors, List ** in_neighbors);
int               topological_sort   (Hypergraph * this, List ** neighbors, List ** in_neighbors, INT * sort);
int               compute_criticality (Hypergraph * h, List ** neighbors, List **  in_neighbors, INT * sort);
int               compute_subpmax     (Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, bool * vertices);
int compute_list_neighbors_unalloc(Hypergraph * this, List ** neighbors_list);
int compute_path_length(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmax, INT * lmax, INT * depth, float * avg, float * stdw);
int compute_pmax(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelpmax, INT * lpmax, float * avg, float * stdw);
int compute_maxdeg(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmaxdeg, INT * lmaxdeg, float * avg, float * stdw);
int compute_maxcon(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmaxdeg, INT * lmaxdeg, float * avg, float * stdw);
int compute_list_in_neighbors_unalloc(Hypergraph * this, List ** neighbors_list, List ** in_neighbors_list);
int compute_neighbors_unalloc(Hypergraph * h, Matrix * neighbors);

int compute_in_neighbors_unalloc(
Hypergraph * h, 
Matrix     * neighbors, 
Matrix     * in_neighbors);

/*
**  The macro definitions.
*/

#define rbh_load(...) var_rbh_load((rbhLoad_args){__VA_ARGS__})



#endif
