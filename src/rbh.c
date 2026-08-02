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
/**   NAME       : rbh.c                                   **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are declarations for the    **/
/**                red-black hypergraph structure.         **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 01 jan 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/************************************************************/


/*
**  The defines and includes.
*/
#include "rbh.h"
#include "vth.h"

/*
**  The static definitions.
*/

/**
 * @brief Function initializing a red-black 
 * hypergraph.    
 *
 * @param this        Red-black hypergraph.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion.
 */
int 
rbh_init(Hypergraph * this)
{
  this->ti_hyperedges = (INT*)calloc((2 * this->i_hyperedges + this->i_pins), sizeof(INT));
  
  MEM_ERROR(this->ti_hyperedges);

  this->ti_idx_hyperedges = (INT*)calloc(this->i_hyperedges, sizeof(INT));
    
  MEM_ERROR(this->ti_idx_hyperedges);

  this->ti_delays = (INT*)calloc(this->i_vertices, sizeof(INT));
    
  MEM_ERROR(this->ti_delays);

  this->ti_criticalities_right = (INT*)calloc(this->i_vertices, sizeof(INT));
  
  MEM_ERROR(this->ti_criticalities_right);

  this->ti_criticalities_left = (INT*)calloc(this->i_vertices, sizeof(INT));
    
  MEM_ERROR(this->ti_criticalities_left);

  this->ti_reds = (INT*)calloc(this->i_reds, sizeof(INT));

  MEM_ERROR(this->ti_reds);

  this->ti_weights = (INT*)calloc(this->i_vertices * this->i_weights, sizeof(INT));
    
  MEM_ERROR(this->ti_weights);

  /* is_red and vth are built after data is loaded (in rbh_load_base) */
  this->is_red = NULL;
  this->vth    = NULL;

  return(0);
}

/**
 * @brief Function deleting a red-black 
 * hypergraph.    
 *
 * @param this        Red-black hypergraph.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
rbh_free(Hypergraph * this)
{
    /* Opt: free precomputed structures if present */
    if (this->is_red) { free(this->is_red); this->is_red = NULL; }
    if (this->vth)    { vth_free(this->vth); this->vth = NULL; }

    free(this->ti_weights);

    free(this->ti_reds);
    
    free(this->ti_criticalities_right);
    
    free(this->ti_criticalities_left);
    
    free(this->ti_delays);
    
    free(this->ti_idx_hyperedges);
    
    free(this->ti_hyperedges);
    
    if(this->s_rbh_name != NULL) 
      {
        free(this->s_rbh_name);
      }

    return(0);
}

 /**
 * @brief Function for building precomputed structures for a red-black hypergraph.
 * Safe to call multiple times; frees any existing data first.
 * Call after filling a Hypergraph's data arrays (e.g. after HEM contraction).
 *       
 *
 * @param h        A red-black hypergraph.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion.
 */
int
rbh_build_precomputed(Hypergraph *h)
{
    INT i;
    if (!h) return 1;
    if (h->is_red) { free(h->is_red); h->is_red = NULL; }
    if (h->vth)    { vth_free(h->vth); h->vth = NULL; }

    h->is_red = (bool *)calloc(h->i_vertices, sizeof(bool));
    if (!h->is_red) return 1;
    for (i = 0; i < h->i_reds; i++)
        h->is_red[h->ti_reds[i]] = true;

    h->vth = vth_build(h);
    if (!h->vth) return 1;

    return 0;
}


/**
 * @brief Function for loading a red-black 
 * hypergraph written in a hygr file format. In 
 * this format, the weight vector is use for 
 * saving for each node the following parameters: 
 * the criticality, the color, the delay and then,
 * the weights.     
 *
 * @param in        A structure containing the following 
 *                  parameters: the hypergraph, file path, 
 *                  vertex indexation, verbose.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
var_rbh_load(rbhLoad_args in)
{
    Hypergraph * h_out             = in.h ? in.h                 : NULL;
    
    char *     s_path_out          = in.s_path ? in.s_path       : NULL;
    
    INT        i_baseval_out       = in.i_baseval ? in.i_baseval : 0;
    
    bool       b_verbose_out       = in.b_verbose ? in.b_verbose : false;

    return rbh_load_base(h_out, s_path_out, i_baseval_out, b_verbose_out);
}

/**
 * @brief Function for validating a red-black hypergraph.
 *
 * Rules checked:
 *  1. Every black (combinational) vertex must have at least one incoming arc.
 *     A black vertex with indeg=0 is undriven — likely a malformed netlist.
 *  2. No purely combinational cycles (cycles with only black vertices).
 *     Red (register) vertices break cycles in synchronous circuits.
 *     Cycles r->b->...->r are valid; b->b->... cycles without red are not.
 *   
 *
 * @param h         A red-black hypergraph. 
 * 
 *  
 * @return number of violations (0 = structurally valid). Messages go to stderr.
 */
int
rbh_validate(const Hypergraph *h)
{
    if (!h || !h->is_red) return -1;

    INT   nv    = h->i_vertices;
    INT   ne    = h->i_hyperedges;
    INT  *e     = h->ti_hyperedges;
    INT  *idx_e = h->ti_idx_hyperedges;
    bool *is_red = h->is_red;
    int   issues = 0;

    /* Compute in-degree of every vertex */
    INT *indeg = (INT *)calloc(nv, sizeof(INT));
    if (!indeg) return -1;
    for (INT j = 0; j < ne; j++) {
        INT base = idx_e[j];
        INT size = e[base + 1];
        for (INT i = 1; i < size; i++) {   /* dsts only (skip src at i=0) */
            INT dst = e[base + 2 + i];
            if (dst >= 0 && dst < nv) indeg[dst]++;
        }
    }

    /* Rule 1: black sources (indeg=0) */
    INT black_sources = 0;
    for (INT u = 0; u < nv; u++)
        if (!is_red[u] && indeg[u] == 0) black_sources++;

    if (black_sources > 0) {
        fprintf(stderr,
            "[rbh_validate] INFO: %d combinational (black) vertices have "
            "no incoming arcs (indeg=0).\n"
            "  These are likely primary inputs (PI) of the circuit — "
            "this is normal.\n",
            black_sources);
        /* Not counted as issues: primary inputs are valid in a partial netlist */
    }

    /* Rule 2: purely combinational cycles — Kahn on black-only subgraph */
    INT  *indeg_b = (INT *)calloc(nv, sizeof(INT));
    bool *vis     = (bool *)calloc(nv, sizeof(bool));
    INT  *queue   = (INT *)malloc((nv + 1) * sizeof(INT));
    if (!indeg_b || !vis || !queue) {
        free(indeg); free(indeg_b); free(vis); free(queue);
        return -1;
    }

    /* Build in-degree restricted to black→black arcs */
    for (INT j = 0; j < ne; j++) {
        INT base = idx_e[j];
        INT size = e[base + 1];
        INT src  = e[base + 2];
        if (src < 0 || src >= nv || is_red[src]) continue;
        for (INT i = 1; i < size; i++) {
            INT dst = e[base + 2 + i];
            if (dst >= 0 && dst < nv && !is_red[dst]) indeg_b[dst]++;
        }
    }

    /* Seed with black vertices that have no black predecessors */
    INT qh = 0, qt = 0;
    for (INT u = 0; u < nv; u++)
        if (!is_red[u] && indeg_b[u] == 0) { queue[qt++] = u; vis[u] = true; }

    /* Process: reduce indeg_b of black successors via VtH index */
    VtH *vth_v = h->vth;
    while (qh < qt) {
        INT u = queue[qh++];
        if (!vth_v) continue;
        /* Iterate over all hyperedges incident to u */
        for (INT k = vth_v->idx[u]; k < vth_v->idx[u + 1]; k++) {
            INT j    = vth_v->data[k];        /* hyperedge index */
            INT base = idx_e[j];
            if (e[base + 2] != u) continue;  /* u must be the source */
            INT size = e[base + 1];
            for (INT i = 1; i < size; i++) {
                INT dst = e[base + 2 + i];
                if (dst < 0 || dst >= nv || is_red[dst] || vis[dst]) continue;
                if (--indeg_b[dst] == 0) { queue[qt++] = dst; vis[dst] = true; }
            }
        }
    }

    INT comb_cycles = 0;
    for (INT u = 0; u < nv; u++)
        if (!is_red[u] && !vis[u]) comb_cycles++;

    if (comb_cycles > 0) {
        fprintf(stderr,
            "[rbh_validate] WARNING: %d combinational (black) vertices are in "
            "purely combinational cycles.\n"
            "  Combinational feedback loops are impossible in synchronous circuits.\n"
            "  The netlist is malformed.\n",
            comb_cycles);
        issues += comb_cycles;
    }

    free(queue); free(vis); free(indeg_b); free(indeg);

    if (issues == 0)
        fprintf(stderr,
            "[rbh_validate] OK: hypergraph structurally valid "
            "(%d vertices, %d red, %d black).\n",
            nv, h->i_reds, nv - h->i_reds);

    return issues;
}


int 
rbh_load_base(Hypergraph *       this, 
            const char * const s_path, 
            INT                i_baseval, 
            bool               b_verbose)
{

  /*
  ** Errors.
  */
  if(this == NULL) 
    {
      return(1);
    }
    
  if(s_path == NULL) 
    {
      return(2);
    }
    
  if(i_baseval < 0) 
    {
      return(3);
    }
  /*
  ** End errors.
  */
  char * buffer = (char*)malloc(BUFSIZE);
  
  MEM_ERROR(buffer);

  FILE * in = NULL;
  in        = fopen(s_path, "r");

  if(in == NULL) 
    {
      return(1);
    }

  if(b_verbose) 
    {
      printf("File : %s\n", s_path);
    }

  INT i_pins, i_vertices, i_hyperedges, i_reds, i_weights;
  
  read_line(in, buffer, BUFSIZE);

  int n = sscanf(buffer, "%d %d %d %d %d", &i_pins, &i_vertices, &i_hyperedges, &i_reds, &i_weights);
    
  if(n < 0) 
    {
      return(10);
    }

  FATAL(n != 5, "Invalid graph format!");
  
  assert(n == 5);

  if(b_verbose) 
    {
      printf("#pins : %d\n#vertex : %d\n#edges : %d\n#red : %d\n", i_pins, i_vertices, i_hyperedges, i_reds);
    }

    this->i_vertices       = i_vertices;
    this->i_hyperedges     = i_hyperedges;
    this->i_reds           = i_reds;
    this->i_pins           = i_pins;
    this->i_weights        = i_weights;

    rbh_init(this);

    char   token[2];
    
    strcpy(token, " ");
    
    char * raw;
    
    char * buffer2 = (char*)malloc(sizeof(char) * BUFSIZE);
    
    MEM_ERROR(buffer2);

    INT  idx = 0;
    INT  len;
    bool is_in;
    INT  cur_len;
    INT  idx_e;

    for(INT j = 0; j < this->i_hyperedges; j++)
      {
        this->ti_idx_hyperedges[j] = idx;
        
        idx_e                      = idx;

        read_line(in, buffer, BUFSIZE);

        if(b_verbose)
          {
            printf("Hyperedges : j(%d)\n%s\n", j, buffer);
          }

        strcpy(buffer2, buffer);

        raw = strtok(buffer, token);
        
        len = 0;

        while(raw != NULL)
          {
            /* len count */
            len++;
            
            raw = strtok(NULL, token);
          }

        raw = strtok(buffer2, token);
        
        /*
        **  weight source sinks
        */
        this->ti_hyperedges[idx++] = atoi(raw);                  /* weight */
        
        raw = strtok(NULL, token);
        
        this->ti_hyperedges[idx++] = len - 1;                    /* hyperedge size */
        
        cur_len                    = 0;
        
        for(INT i = 0; i < len - 1; i++) 
          {
            is_in                      = false;
            this->ti_hyperedges[idx++] = atoi(raw) - i_baseval;  /* vertex */
            
            for(INT ii = 0; ii < cur_len; ii++) 
              {
                if(this->ti_hyperedges[idx_e + 2 + ii] == this->ti_hyperedges[idx - 1]) 
                  {
                    is_in = true;
                  }
              }
            cur_len++;

            if(is_in) 
              {
                cur_len--;
                idx--;
              }
            raw = strtok(NULL, token);
          }
        this->ti_hyperedges[idx_e + 1] = cur_len;
      }
    idx = 0;

    for(INT i = 0; i < this->i_vertices; i++) 
      {
        /*
        ** red (0/1) delay criticality weights
        */
        read_line(in, buffer, BUFSIZE);
        if(b_verbose) 
          {
            printf("vertex : %d\t %s\n", i, buffer);
          }
        raw = strtok(buffer, token);
        
        if(atoi(raw) == 1) 
          {
            this->ti_reds[idx++] = i;
          }
        raw = strtok(NULL, token);
        
        this->ti_delays[i] = atoi(raw);
        
        raw = strtok(NULL, token);
        
        this->ti_criticalities_right[i] = atoi(raw);

        for(INT iw = 0; iw < this->i_weights; iw ++) 
          {
            raw = strtok(NULL, token);
            this->ti_weights[i * i_weights + iw] = atoi(raw);
          }
      }

    free(buffer2);
    free(buffer);
    fclose(in);

    this->is_red = (bool *)calloc(this->i_vertices, sizeof(bool));
    MEM_ERROR(this->is_red);
    for (INT _i = 0; _i < this->i_reds; _i++)
        this->is_red[this->ti_reds[_i]] = true;

    this->vth = vth_build(this);
    MEM_ERROR(this->vth);

    return(0);
}

/**
 * @brief Function for saving a red-black 
 * hypergraph written in a hygr file format. In 
 * this format, the weight vector is use for 
 * saving for each node the following parameters: 
 * the criticality, the color, the delay and then,
 * the weights.     
 *
 * @param this         A red-black hypergraph. 
 * @param s_path       String (path of the loaded file).
 * @param i_baseval    Integer (base value for node index).
 * @param b_verbose    Boolean (if true, print more informations).
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
rbh_save(Hypergraph *       this, 
        const char * const s_path, 
        INT                i_baseval, 
        bool               b_verbose)
{
    /*
    ** Errors.
    */
    if(this == NULL) 
      {
        return(1);
      }

    if(s_path == NULL) 
      {
        return(2);
      }

    if(i_baseval < 0) 
      {
        return(3);
      }
    /*
    ** End errors.
    */

    FILE * out = NULL;
    out        = fopen(s_path, "w+");

    if(out == NULL) 
      {
        return(1);
      }

    if(b_verbose) 
      {
        printf("File : %s\n", s_path);
      }

    int n = fprintf(out, "%d %d %d %d %d\n", this->i_pins, this->i_vertices, this->i_hyperedges, this->i_reds, this->i_weights);

    if(n < 0) 
      {
        return(10);
      }

    if(b_verbose) 
      {
        printf("#pins : %d\n#vertex : %d\n#edges : %d\n#red : %d\n", this->i_pins, this->i_vertices, this->i_hyperedges, this->i_reds);
      }

    
    bool *is_red = this->is_red;

    INT idx, len;

    for(INT j = 0; j < this->i_hyperedges; j++) 
      {
        idx = this->ti_idx_hyperedges[j];

        len = this->ti_hyperedges[idx+1];

        fprintf(out, "%d", this->ti_hyperedges[idx++]);
        idx++;
        
        for(INT i = 0; i < len; i++) 
          {
            fprintf(out, " %d", this->ti_hyperedges[idx++]);
          }
        fprintf(out, "\n");
      }

    idx = 0;

    for(INT i = 0; i < this->i_vertices; i++) 
      {
        /*
        ** red (0/1) delay criticality weights
        */
        if(is_red[i]) 
          {
            fprintf(out, "1");
          }
        else 
          {
            fprintf(out, "0");
          }

        fprintf(out, " %d", this->ti_delays[i]);
        
        fprintf(out, " %d", this->ti_criticalities_right[i]);
        
        for(INT iw = 0; iw < this->i_weights; iw ++) 
          {
            fprintf(out, " %d", this->ti_weights[i * this->i_weights + iw]);
          }
        fprintf(out, "\n");
      }

    
    fclose(out);

    return(0);
}

/**
 * @brief Function computing the list of 
 * outneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a list (neighbors_list).
 *
 * @param this            A red-black hypergraph. 
 * @param neighbors_list  List.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_list_neighbors(Hypergraph * this, 
                     List      ** neighbors_list) 
{
    INT nv    = this->i_vertices;
    INT ne    = this->i_hyperedges;
    INT *e    = this->ti_hyperedges;
    INT *idx_e = this->ti_idx_hyperedges;

    INT *seen    = (INT *)malloc(sizeof(INT) * nv); MEM_ERROR(seen);
    INT *touched = (INT *)malloc(sizeof(INT) * nv); MEM_ERROR(touched);
    for(INT i = 0; i < nv; i++) seen[i] = -1;

    for(INT j = 0; j < ne; j++) {
        INT idx_j  = idx_e[j];
        INT size_j = e[idx_j + 1];
        INT u      = e[idx_j + 2];
        INT n_touched = 0;
        for(INT i = 0; i < size_j - 1; i++) {
            INT v = e[idx_j + 3 + i];
            if(seen[v] != u) {
                seen[v]              = u;
                touched[n_touched++] = v;
                list_add_element(neighbors_list[u], v);
            }
        }
        for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
    }
    return(0);
}

/**
 * @brief Function computing the list of 
 * inneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a list (in_neighbors_list).
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors_list    List.
 * @param in_neighbors_list List.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_list_in_neighbors(Hypergraph * this, 
                       List      ** neighbors_list, 
                       List      ** in_neighbors_list) 
{
    INT nv = this->i_vertices;

    INT *seen    = (INT *)malloc((size_t)nv * sizeof(INT));
    INT *touched = (INT *)malloc((size_t)nv * sizeof(INT));
    MEM_ERROR(seen);
    MEM_ERROR(touched);
    for(INT i = 0; i < nv; i++) seen[i] = -1;

    for(INT v = 0; v < nv; v++) {
        List *cell = neighbors_list[v];
        INT n_touched = 0;
        while(cell != NULL && cell->size > 0) {
            INT u = cell->i;
            if(seen[u] != v) {
                seen[u]              = v;
                touched[n_touched++] = u;
                list_add_element(in_neighbors_list[u], v);
            }
            cell = cell->next;
        }
        for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
    }
    
    free(touched);
    free(seen);
    return(0);
}

/**
 * @brief Function computing the list of 
 * inneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a matrix (neighbors).
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors         Matrix.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_neighbors(Hypergraph * this, 
                 Matrix     * neighbors) 
{
  INT nv       = this->i_vertices;           /* number of vertices.                        */
  INT ne       = this->i_hyperedges;         /* number of hyperarcs.                       */
  INT * e      = this->ti_hyperedges;      /* array containing the hyperegdes.           */
  INT * idx_e  = this->ti_idx_hyperedges;  /* array containing the index of hyperedges.  */

  new_matrix(neighbors, nv, nv);

  for(INT i = 0; i < nv; i++)
      neighbors->v[i].size = 0;

  INT *seen    = (INT *)malloc((size_t)nv * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)nv * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < nv; i++) seen[i] = -1;

  for(INT j = 0; j < ne; j++) {
      INT idx_j  = idx_e[j];
      INT size_j = e[idx_j + 1];
      INT u      = e[idx_j + 2];
      INT n_touched = 0;

      for(INT i = 0; i < size_j - 1; ++i) {
          INT v = e[idx_j + 3 + i];
          if(seen[v] != u) {
              seen[v]              = u;
              touched[n_touched++] = v;
              /* grow row if needed */
              if(neighbors->v[u].size == neighbors->n) {
                  INT new_n = neighbors->n * 2;
                  for(INT x = 0; x < neighbors->m; x++) {
                      void *tmp = realloc(neighbors->v[x].v, sizeof(INT) * new_n);
                      if(!tmp) { free(touched); free(seen); return(2); }
                      neighbors->v[x].v = (INT *)tmp;
                  }
                  neighbors->n = new_n;
              }
              neighbors->v[u].v[neighbors->v[u].size++] = v;
          }
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }

  free(touched);
  free(seen);
  return(0);
}

/**
 * @brief Function computing the list of 
 * inneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a matrix (neighbors).
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors         Matrix.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_in_neighbors(Hypergraph * this, 
                   Matrix     * neighbors, 
                   Matrix     * in_neighbors) 
{
  INT   nv     = this->i_vertices;
  
  new_matrix(in_neighbors, nv, nv);

  for(INT v = 0; v < nv; v++)
      in_neighbors->v[v].size = 0;

  /* Fix: O(1) seen[] marker replaces O(in_deg) linear scan */
  INT *seen    = (INT *)malloc((size_t)nv * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)nv * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < nv; i++) seen[i] = -1;

  for(INT v = 0; v < nv; v++) {
      INT n_touched = 0;
      for(INT i = 0; i < neighbors->v[v].size; i++) {
          INT u = neighbors->v[v].v[i];
          if(seen[u] != v) {
              seen[u]              = v;
              touched[n_touched++] = u;
              in_neighbors->v[u].v[in_neighbors->v[u].size++] = v;
          }
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }
  free(touched);
  free(seen);
  return(0);
}

/**
 * @brief Function computing a topological 
 * sort of the vertices of the red-black 
 * hypergraph (this) and saving it in an
 * array of integers. 
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors         List.
 * @param in_neighbors      List.
 * @param sort              Array of integers.
 * 
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
topological_sort(Hypergraph * this, 
                List      ** neighbors, 
                List      ** in_neighbors, 
                INT        * sort) 
{
  
  INT nv     = this->i_vertices;    /* number of vertices                        */
  INT ne     = this->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = this->i_reds;        /* number of red vertices                    */
  INT np     = this->i_pins;        /* number of pins                            */
  INT nwv    = this->i_weights;     /* number of vertices weight                 */
  INT * reds = this->ti_reds;       /* arrays of red vertices                    */

  INT i,u,v;

  INT  *queue      = (INT *)malloc((size_t)nv * sizeof(INT));
  INT  *main_queue = (INT *)malloc((size_t)nv * sizeof(INT));
  INT  *indeg      = (INT *)malloc((size_t)nv * sizeof(INT));
  bool *flag       = (bool *)malloc((size_t)nv * sizeof(bool));
  MEM_ERROR(queue);
  MEM_ERROR(main_queue);
  MEM_ERROR(indeg);
  MEM_ERROR(flag);

  bool *is_red = this->is_red;

  for(i = 0; i < nv; i++) 
    {
      indeg[i]   = in_neighbors[i]->size;
      flag[i]    = false;
      queue[i]   = -1;
      main_queue[i] = -1;
    }

    INT main_queue_start = 0;
    
    INT main_queue_end   = 0;
    
    for(u = 0; u < nv; u++) 
      {
        if(is_red[u]) 
          {
            main_queue[main_queue_end++] = u;
            flag[u]                      = true;
          }
      }
    
    INT index = 0, queue_start = 0, queue_end = 0;

    while(main_queue_start != main_queue_end) 
      {
        queue_start = 0;
        queue_end = 0;

        for(i = 0; i < main_queue_end; i++) 
          {
            // all red vertices visited are explored next
            queue[queue_end++] = main_queue[i];
          }

        main_queue_start = 0;
        
        main_queue_end   = 0;
        
        while(queue_start != queue_end) 
          {
            u             = queue[queue_start++];
            
            sort[index++] = u;
            
            flag[u]       = true;
            
            List * cell   = neighbors[u];
            
            for(INT x = 0; x < neighbors[u]->size; x++) 
              {
                v = cell->i;
                
                indeg[v]--;
                
                if(indeg[v] <= 0 && !is_red[v] && !flag[v]) 
                  {
                    flag[v]            = true;
                    queue[queue_end++] = v;
                  }
                
                if(is_red[v] && !flag[v]) 
                  {
                    main_queue[main_queue_end++] = v;
                    flag[v]                      = true;
                  }
                cell = cell->next;
              }
          }
      }

    if(index < nv) 
      {
        // if reds vertex u,v st.
        //  u --> v --> u
        main_queue_start = 0;
        
        main_queue_end   = 0;
        
        for(u = 0; u < nv; u++) 
          {
            if(is_red[u] && !flag[u]) 
              {
                main_queue[main_queue_end++] = u;
                flag[u]                      = true;
              }
            
            if(indeg[u] == 0 && !flag[u]) 
              {
                main_queue[main_queue_end++] = u;
                flag[u]                      = true;
              }
          }
        /* if exists interdependance between u --> v --> u with v black */
        /* !flag[v] = true but indeg[v]>0  */
        for(u = 0; u < nv; u++) 
          {
            if(flag[u] && is_red[u]) 
              {
                List * cell   = neighbors[u];

                for(INT x = 0; x < neighbors[u]->size; x++) 
                  {
                    v = cell->i;
                    
                    if(!flag[v])
                      {
                        main_queue[main_queue_end++] = v;
                        flag[v]                      = true;
                      }
                    cell = cell->next;
                  }
              }
          }

        while(main_queue_start != main_queue_end) 
          {
            queue_start = 0;
            
            queue_end = 0;
            
            for(i = 0; i < main_queue_end; i++) 
              {
                // all red vertices visited are explored next
                queue[queue_end++] = main_queue[i];
              }

            main_queue_start = 0;
            
            main_queue_end   = 0;
            
            while(queue_start != queue_end) 
              {
                u             = queue[queue_start++];
                
                sort[index++] = u;
                
                flag[u]       = true;
                
                List * cell   = neighbors[u];
                
                for(INT x = 0; x < neighbors[u]->size; x++) 
                  {
                    v = cell->i;
                    
                    if(!is_red[v] && !flag[v]) 
                      {
                        flag[v]            = true;
                        queue[queue_end++] = v;
                      }
                    
                    if(is_red[v] && !flag[v]) 
                      {
                        main_queue[main_queue_end++] = v;
                        flag[v]                      = true;
                      }
                    cell = cell->next;
                  }
              }
          }
      }
    /* In the case of black cycle, we push the vertices in the end
       of topological sort to not taking account these vertices */
    for(INT i = 0; i < nv; i++) 
      {
        if(!flag[i]) 
          {
            sort[index++] = i;
          }
      }

    assert(index == nv);

    free(flag);
    free(indeg);
    free(main_queue);
    free(queue);
    return(0);
}

/**
 * @brief Function computing the list of 
 * outneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a list 
 * (neighbors_list). The list is not allocated.
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors_list    List.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_list_neighbors_unalloc(Hypergraph * this, 
                            List      ** neighbors_list) 
{
  INT nv      = this->i_vertices;
  INT ne      = this->i_hyperedges;
  INT *e      = this->ti_hyperedges;
  INT *idx_e  = this->ti_idx_hyperedges;

  for(INT i = 0; i < nv; i++) {
      
      List *_cur = neighbors_list[i]->next;
      while(_cur != NULL) {
          List *_nxt = _cur->next;
          free(_cur);
          _cur = _nxt;
      }
      neighbors_list[i]->next = NULL;
      neighbors_list[i]->size = 0;
      neighbors_list[i]->i    = -1;
  }

  INT seen_size = nv;
  for(INT _j = 0; _j < ne; _j++) {
      INT _base = idx_e[_j], _sz = e[_base + 1];
      for(INT _i = 0; _i < _sz; _i++) {
          INT _vid = e[_base + 2 + _i];
          if(_vid >= seen_size) seen_size = _vid + 1;
      }
  }
  
  INT *seen    = (INT *)malloc((size_t)seen_size * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)seen_size * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < seen_size; i++) seen[i] = -1;

  INT idx_j, size_j, u, v;

  for(INT j = 0; j < ne; j++) {
      idx_j  = idx_e[j];
      size_j = e[idx_j + 1];
      u      = e[idx_j + 2];

      INT n_touched = 0;
      for(INT i = 0; i < size_j - 1; i++) {
          v = e[idx_j + 3 + i];
          if(seen[v] != u) {
              seen[v]              = u;
              touched[n_touched++] = v;
              list_add_element(neighbors_list[u], v);
          }
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }

  free(touched);
  free(seen);
  return(0);
}

/**
 * @brief Function computing the list of 
 * inneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a list 
 * (in_neighbors_list). The list was already allocated.
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors_list    List.
 * @param in_neighbors_list List.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_list_in_neighbors_unalloc(Hypergraph * this, 
                              List ** neighbors_list, 
                              List ** in_neighbors_list) 
{
  INT nv = this->i_vertices;

  for(INT i = 0; i < nv; i++) {
      List *_ci = in_neighbors_list[i]->next;
      while(_ci != NULL) { List *_ni = _ci->next; free(_ci); _ci = _ni; }
      in_neighbors_list[i]->next = NULL;
      in_neighbors_list[i]->size = 0;
  }

  INT seen_size = nv;
  for(INT v = 0; v < nv; v++) {
      List *_tmp = neighbors_list[v];
      for(INT _x = 0; _x < neighbors_list[v]->size; _x++) {
          if(_tmp->i >= seen_size) seen_size = _tmp->i + 1;
          _tmp = _tmp->next;
      }
  }
  
  INT *seen    = (INT *)malloc((size_t)seen_size * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)seen_size * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < seen_size; i++) seen[i] = -1;

  for(INT v = 0; v < nv; v++) {
      List *cell     = neighbors_list[v];
      INT   deg_v    = neighbors_list[v]->size;

      INT n_touched = 0;
      for(INT i = 0; i < deg_v; i++) {
          INT u = cell->i;
          if(seen[u] != v) {
              seen[u]              = v;
              touched[n_touched++] = u;
              list_add_element(in_neighbors_list[u], v);
          }
          cell = cell->next;
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }

  free(touched);
  free(seen);
  return(0);
}

/**
 * @brief Function computing the list of 
 * inneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a matrix 
 * (in_neighbors). The list was already allocated.
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors         Matrix.
 * @param in_neighbors      Matrix.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_in_neighbors_unalloc(Hypergraph * h, 
                          Matrix     * neighbors, 
                          Matrix     * in_neighbors) 
{   
  INT   nv = h->i_vertices;
  INT   ne = h->i_hyperedges;
  
  in_neighbors->m = nv;
  in_neighbors->n = nv;

  for(INT v = 0; v < nv; v++)
      in_neighbors->v[v].size = 0;

  /* Fix: O(1) seen[] marker replaces O(in_deg) linear scan */
  INT *seen    = (INT *)malloc((size_t)nv * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)nv * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < nv; i++) seen[i] = -1;

  for(INT v = 0; v < nv; v++) {
      INT n_touched = 0;
      for(INT i = 0; i < neighbors->v[v].size; i++) {
          INT u = neighbors->v[v].v[i];
          if(seen[u] != v) {
              seen[u]              = v;
              touched[n_touched++] = u;
              in_neighbors->v[u].v[in_neighbors->v[u].size++] = v;
          }
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }
  free(touched);
  free(seen);
  return 0;
}

/**
 * @brief Function computing the list of 
 * outneighbors of the red-black hypergraph in
 * paramter (this) and saving it in a matrix 
 * (neighbors). The matrix was already allocated.
 *
 * @param this              A red-black hypergraph. 
 * @param neighbors         Matrix.
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_neighbors_unalloc(Hypergraph * h, 
                        Matrix * neighbors) 
{
  INT   nv = h->i_vertices;
  INT   ne = h->i_hyperedges;
  
  INT * idx_hyperedges = h->ti_idx_hyperedges;
  INT * hyperedges     = h->ti_hyperedges;

  neighbors->m = nv;
  neighbors->n = nv;

  for(INT i = 0; i < nv; i++) 
    {
      neighbors->v[i].size = 0;
    }

  INT u, v, idx_j, weight, size_j;

  /* Fix: O(1) seen[] marker replaces O(deg) linear scan */
  INT *seen    = (INT *)malloc((size_t)nv * sizeof(INT));
  INT *touched = (INT *)malloc((size_t)nv * sizeof(INT));
  MEM_ERROR(seen);
  MEM_ERROR(touched);
  for(INT i = 0; i < nv; i++) seen[i] = -1;

  for(INT j = 0; j < ne; j++) {
      idx_j  = idx_hyperedges[j];
      weight = hyperedges[idx_j];
      size_j = hyperedges[idx_j + 1];
      u      = hyperedges[idx_j + 2];
      INT n_touched = 0;

      for(INT i = 1; i < size_j; ++i) {
          v = hyperedges[idx_j + 2 + i];
          if(seen[v] != u) {
              seen[v]              = u;
              touched[n_touched++] = v;
              if(neighbors->v[u].size == neighbors->n) {
                  INT new_n = neighbors->n * 2;
                  for(int x = 0; x < neighbors->m; x++) {
                      void *tmp = realloc(neighbors->v[x].v, sizeof(INT) * new_n);
                       if(!tmp) { free(touched); free(seen); return(2); }
                      neighbors->v[x].v = (INT *)tmp;
                  }
                  neighbors->n = new_n;
              }
              neighbors->v[u].v[neighbors->v[u].size++] = v;
          }
      }
      for(INT t = 0; t < n_touched; t++) seen[touched[t]] = -1;
  }
  free(touched);
  free(seen);
  return 0;
}

/**
 * @brief Function computing the criticality
 * for each vertices in the red-black hypergraph. 
 * The criticality is a weight proportional to the
 * maximum length of paths passing through the 
 * vertices.
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 *
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_criticality(Hypergraph * h, 
                    List      ** neighbors, 
                    List      ** in_neighbors, 
                    INT        * sort) 
{
  bool verbose = false;

  INT nvp    = 0;
  INT nep    = 0;
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  bool *is_red = h->is_red;

  bool *flag       = (bool *)malloc((size_t)nv * sizeof(bool));
  INT  *delays_fwd = (INT *)malloc((size_t)nv * sizeof(INT));
  INT  *delays_bwd = (INT *)malloc((size_t)nv * sizeof(INT));
  MEM_ERROR(flag);
  MEM_ERROR(delays_fwd);
  MEM_ERROR(delays_bwd);

  INT crit_max = 0;
    
  for(INT i = 0; i < nv; i++) 
    {
      delays_fwd[i] = 0;
      
      delays_bwd[i] = 0;
      
      flag[i]       = false;
    }

    INT u,v;
    
    for(INT i = 0; i < nv; i++) {
        
        v       = sort[i];
        
        flag[v] = true;
        
        if(is_red[v]) 
          {
            List * cell = neighbors[v];
            
            while(cell != NULL && neighbors[v]->size) 
              {
                u = cell->i;
                
                if(!is_red[u]) 
                  {
                    delays_fwd[u] = MAX(delays_fwd[u], h->ti_delays[v] + h->ti_delays[u]);
                  }
                cell = cell->next;
              }
          }
        else 
          {
            List * cell = neighbors[v];
            while(cell != NULL && neighbors[v]->size) 
              {
                u = cell->i;
                
                if(!is_red[u]) 
                  {
                    delays_fwd[u] = MAX(delays_fwd[u], delays_fwd[v] + h->ti_delays[u]);
                  }
                cell = cell->next;
              }
          }
      }

    INT dhat;

    // retro propagation
    for(INT i = nv - 1; i > 0; i--) 
      {
        v       = sort[i];
       
        dhat    = 0;
        
        if(is_red[v] || delays_bwd[v] == 0)
          {
            delays_bwd[v] = delays_fwd[v];
          }

        List * cell = in_neighbors[v];
        
        while(cell != NULL && in_neighbors[v]->size) 
          { 
            u = cell->i;
            
            if(dhat < delays_fwd[u]) 
              {
                dhat = delays_fwd[u];
              }
            cell = cell->next;
          }
        cell = in_neighbors[v];
        
        while(cell != NULL && in_neighbors[v]->size) 
          {
            u = cell->i;
            
            if(!is_red[u]) 
              {
                delays_bwd[u] = MAX(delays_bwd[u], delays_bwd[v] - (dhat - delays_fwd[u]));
              }
            cell = cell->next;
          }
      }

    for(v = 0; v < nv; v++) 
      {
        h->ti_criticalities_right[v] = MAX(delays_bwd[v], 1);
        
        h->ti_criticalities_left[v]  = MAX(delays_fwd[v], 1);
      }

    for(INT j = 0; j < ne; j++) 
      {
        INT idx_e  = h->ti_idx_hyperedges[j];
        
        INT source = h->ti_hyperedges[idx_e + 2];
        
        h->ti_hyperedges[idx_e] = h->ti_criticalities_right[source];
      }

      
  
    free(delays_bwd);
    free(delays_fwd);
    free(flag);
    return 0;
}


/**
 * @brief Function computing the maximum 
 * connectinivity of hyperedges in the 
 * red-black hypergraph. 
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 * @param sizelmaxcon    Integer (maximum hyperedge size).
 * @param lmaxcon        Array (hyperedge with maximum size).
 * @param avg            Integer (average of hyperedges size).
 * @param stdw           Float (standard deviation of the hyperedge size).
 * 
 * @return Integer, 0 if no probelm occurs during the 
 * execution of this funtion..
 */
int 
compute_maxcon(Hypergraph * h, 
               List      ** neighbors, 
               List      ** in_neighbors, 
               INT        * sort, 
               INT        * sizelmaxcon, 
               INT        * lmaxcon, 
               float      * avg, 
               float      * stdw) 
{  
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  INT max_con=0;
  
  (*avg) = 0;
  
  float var = 0;
  
  INT * e      = h->ti_hyperedges;
  
  INT * idx_e  = h->ti_idx_hyperedges;

  for(INT j = 0; j < ne; j++) 
    {
     INT idx_j  = idx_e[j];
     
     INT size_j = e[idx_j + 1];
     
     (*avg) += size_j;
     
     var += size_j * size_j;
     
     if(max_con < size_j) 
       {
         max_con = size_j;
       }  
    }
  
  (*avg) /= (float)nv;

  (*stdw) = sqrt(var/(float)nv - (*avg) * (*avg)); 
   
  for(INT j = 0; j < ne; j++) 
    {
      
      INT idx_j  = idx_e[j];
      
      INT size_j = e[idx_j + 1];
     
      if(max_con == size_j) 
        {
          lmaxcon[(*sizelmaxcon)++] = j;
        }
    }
  return max_con;
}

/**
 * @brief Function computing the maximum 
 * degree of vertices in the 
 * red-black hypergraph. 
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 * @param sizelmaxcon    Integer (maximum hyperedge size).
 * @param lmaxcon        Array (hyperedge with maximum size).
 * @param avg            Integer (average of hyperedges size).
 * @param stdw           Float (standard deviation of the hyperedge size).
 * 
 * @return Integer, the maximum degree of vertices. 
 */
int 
compute_maxdeg(Hypergraph * h, 
               List      ** neighbors, 
               List      ** in_neighbors, 
               INT        * sort, 
               INT        * sizelmaxdeg, 
               INT        * lmaxdeg, 
               float      * avg, 
               float      * stdw) 
{
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  INT *degrees = (INT*)malloc(nv * sizeof(INT)); MEM_ERROR(degrees);

  INT max_deg = 0;
  
  (*avg) = 0;
  
  float var = 0;

  for(INT i = 0; i<nv; i++)
    {
      degrees[i] = neighbors[i]->size;
      degrees[i] += in_neighbors[i]->size;
    
      (*avg) += degrees[i];
      
      var += degrees[i] * degrees[i];
    
      if(degrees[i] > max_deg) 
        {
          max_deg = degrees[i];
        }
    }
  
  (*avg) /= (float)nv;
  
  (*stdw) = sqrt(var/(float)nv - (*avg) * (*avg));
  for(INT v = 0; v<nv; v++)
    {
	  if(degrees[v] == max_deg)
        {
		  lmaxdeg[(*sizelmaxdeg)++] = v;
		}
    }
  free(degrees);

  return max_deg;
}

/**
 * @brief Function computing the maximum 
 * path length by propagating the delay of vertices 
 * in the red-black hypergraph. 
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 * @param sizelmaxcon    Integer (maximum hyperedge size).
 * @param lmaxcon        Array (hyperedge with maximum size).
 * @param avg            Integer (average of hyperedges size).
 * @param stdw           Float (standard deviation of the hyperedge size).
 * 
 * @return Integer, the maximum path delay. 
 */
int 
compute_pmax(Hypergraph * h, 
             List ** neighbors, 
             List ** in_neighbors, 
             INT * sort, 
             INT * sizelpmax, 
             INT * lpmax,  
             float * avg, 
             float * stdw)  
{
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  INT  *delays = (INT *)malloc((size_t)nv * sizeof(INT));
  bool *is_red = (bool *)malloc((size_t)nv * sizeof(bool));
  bool *flag   = (bool *)malloc((size_t)nv * sizeof(bool));
  MEM_ERROR(delays);
  MEM_ERROR(is_red);
  MEM_ERROR(flag);

  INT   u, v;

  INT ncuts = 0;
  INT idx_pmax = 0;

  (*avg) = 0;

  for(INT i = 0; i< nv; i++)
    {
      delays[i] =  h->ti_delays[i];
      
      flag[i]   = false;
      
      is_red[i]    = false;
    }

  for(INT i = 0; i < nr; i++) 
    {
      is_red[reds[i]] = true;
    }

  for(INT i = 0; i<nv; i++)
    {
      v       = sort[i];
    
      flag[v] = true;

      if(is_red[v])
        {
          List * cell_neighbors = neighbors[v];
            
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;
                       
              delays[u] = MAX(delays[u], h->ti_delays[v] + h->ti_delays[u]);

              cell_neighbors = cell_neighbors->next;
            }
        }
      else
        {
          List * cell_neighbors = neighbors[v];
          
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;

              delays[u] = MAX(delays[u], delays[v] + h->ti_delays[u]);
                  
              cell_neighbors = cell_neighbors->next;
            }
        }
    }

  int max_v = -1;
  
  int max_delays = 0;
  
  float var = 0;
  
  for(INT v = 0; v<nv; v++)
    {
      (*avg) += delays[v];
      
      var        += delays[v] * delays[v];
		
      if(delays[v] > max_delays)
        {
		  max_delays = delays[v];
		  
          max_v = v;
		}
    }
  
  (*avg) /= (float)nv;
  
  (*stdw) = sqrt(var / (float)nv - (*avg) * (*avg));
  
  for(INT v = 0; v < nv; v++)
    {
	  if(delays[v] == max_delays)
        {
		  lpmax[idx_pmax++] = v;
		}
    }
  
  (*sizelpmax) = idx_pmax;


  free(flag);
  free(is_red);
  free(delays);
  return max_delays;
}

/**
 * @brief Function computing the maximum 
 * path length by counting the number of vertices  
 * along each path in each DAH of the in the 
 * red-black hypergraph. 
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 * @param sizelmax       Integer (Adress).
 * @param lmax           Integer (Adress).
 * @param depth          Integer (Adress).
 * @param avg            Integer (average of hyperedges size).
 * @param stdw           Float (standard deviation of the hyperedge size).
 * 
 * @return Integer, the value of the longest path. 
 */
int 
compute_path_length(Hypergraph * h, 
                    List      ** neighbors, 
                    List      ** in_neighbors, 
                    INT        * sort, 
                    INT        * sizelmax, 
                    INT        * lmax, 
                    INT        * depth, 
                    float      * avg, 
                    float      * stdw) 
{
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  bool *is_red = (bool *)malloc((size_t)nv * sizeof(bool));
  MEM_ERROR(is_red);

  INT   u, v;

  INT ncuts = 0;
  INT idx_lmax = 0;
  
  (*avg) = 0;

  for(INT i = 0; i < nv; i++)
    {
      depth[i] = 0;
      is_red[i]    = false;
    }

  for(INT i = 0; i < nr; i++) 
    {
      is_red[reds[i]] = true;
    }
  
  /* Init depth[] - caller may not zero-init */
  for(INT i = 0; i < nv; i++) depth[i] = 0;

  for(INT i = 0; i < nv; i++)
    {
      v = sort[i];
      
      if(is_red[v])
        {
          List * cell_neighbors = neighbors[v];
        
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;
                       
              depth[u] = MAX(depth[u], 2);

              cell_neighbors = cell_neighbors->next;
            }
        }
      else
        {
          List * cell_neighbors = neighbors[v];
        
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;

              depth[u] = MAX(depth[u], depth[v] + 1);
                  
              cell_neighbors = cell_neighbors->next;
            }
        }
    }  
      
  INT max_length = 0;
  float var = 0;
  
  for(INT v = 0; v < nv; v++)
    {
       (*avg) += depth[v];
       
       var += (float)depth[v] * (float)depth[v];
		
       if(depth[v] > max_length)
         {
			max_length = depth[v] ;	
		 }
    }
  
  (*avg) /= (float)nv;
  
  (*stdw) = sqrt(var / (float)nv - (*avg) * (*avg));
  
  for(INT v = 0; v < nv; v++)
    {
	  if(depth[v] == max_length)
        {
		  lmax[idx_lmax++] = v;	
		}
    }

  (*sizelmax) = idx_lmax;
  
   
  free(is_red);
  return max_length;
}

/**
 * @brief Function computing the maximum 
 * path delay length on vetices in a subhypergraph. 
 * The vertices are identified by an array of boolean, 
 * indicating if a vertex i belong to the subhypergraph. 
 * The delay length is computed by propagating the delay along
 * the topological sort of the vertices.  
 *
 * @param h              A red-black hypergraph. 
 * @param neighbors      List.
 * @param in_neighbors   List.
 * @param sort           Array of integer (topologial sort of vertices).
 * @param vertices       Array of boolean (subhypergraph).
 * 
 * @return Integer, the value of the longest path. 
 */
int 
compute_subpmax(Hypergraph * h, 
                List      ** neighbors, 
                List      ** in_neighbors, 
                INT        * sort, 
                bool       * vertices) 
{
  
  INT nv     = h->i_vertices;    /* number of vertices                        */
  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h->i_reds;        /* number of red vertices                    */
  INT np     = h->i_pins;        /* number of pins                            */
  INT nwv    = h->i_weights;     /* number of vertices weight                 */
  INT * reds = h->ti_reds;       /* arrays of red vertices                    */

  INT  *delays = (INT *)malloc((size_t)nv * sizeof(INT));
  bool *is_red = (bool *)malloc((size_t)nv * sizeof(bool));
  bool *flag   = (bool *)malloc((size_t)nv * sizeof(bool));
  MEM_ERROR(delays);
  MEM_ERROR(is_red);
  MEM_ERROR(flag);

  INT u, v;

  INT ncuts = 0;

  for(INT i = 0; i < nv; i++)
    {
      delays[i] = 0;
    
      if(vertices[i])
        {
          delays[i] = h->ti_delays[i];
        }

      flag[i]   = false;
      
      is_red[i]    = false;
    }

  for(INT i = 0; i < nr; i++) 
    {
      is_red[reds[i]] = true;
    }

  for(INT i = 0; i < nv; i++)
    {
      v       = sort[i];
    
      flag[v] = true;

      if(is_red[v])
        {
          List * cell_neighbors = neighbors[v];
            
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;
                        
              if(vertices[u] && vertices[v])
                {
                  delays[u] = MAX(delays[u], h->ti_delays[v] + h->ti_delays[u]);
                }
              
              cell_neighbors = cell_neighbors->next;
            }
        }
      else
        {
          List * cell_neighbors = neighbors[v];
      
          while(cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;
                  
              if(vertices[u] && vertices[v])
                {
                  delays[u] = MAX(delays[u], delays[v] + h->ti_delays[u]);
                }
              cell_neighbors = cell_neighbors->next;
            }
        }
    }

  int max_v = -1;
  
  int max_delays = 0;
  
  for(INT v = 0; v < nv; v++)
    {
	  if(delays[v] > max_delays)
        {
		  max_delays = delays[v];
			
          max_v = v;
		}
    }


  free(flag);
  free(is_red);
  free(delays);
  return max_delays;
}
