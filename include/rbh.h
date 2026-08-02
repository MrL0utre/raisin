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

  INT               (* rbh_init)  ();           /*+ Hypergraph init function.    +*/
  INT               (* rbh_free)  ();           /*+ Hypergraph free function.    +*/
  INT               (* rbh_load)  ();            /*+ Hypergraph loading function. +*/
  INT               (* rbh_save)  ();           /*+ Hypergraph saving function.  +*/

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
