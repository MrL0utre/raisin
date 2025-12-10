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
/**                                 to   : 13 oct 2023     **/
/**                                                        **/
/************************************************************/


/*
**  The defines and includes.
*/
#include "rbh.h"

/*
**  The static definitions.
*/

int rbhInit(Hypergraph * this){

    this->ti_hyperedges           = (INT*)calloc((2*this->i_hyperedges+this->i_pins), sizeof(INT));
    MEM_ERROR(this->ti_hyperedges);

    this->ti_idx_hyperedges       = (INT*)calloc(this->i_hyperedges, sizeof(INT));
    MEM_ERROR(this->ti_idx_hyperedges);

    this->ti_delays               = (INT*)calloc(this->i_vertices, sizeof(INT));
    MEM_ERROR(this->ti_delays);

    this->ti_criticalities_right  = (INT*)calloc(this->i_vertices, sizeof(INT));
    MEM_ERROR(this->ti_criticalities_right);

    this->ti_criticalities_left   = (INT*)calloc(this->i_vertices, sizeof(INT));
    MEM_ERROR(this->ti_criticalities_left);

    this->ti_reds                 = (INT*)calloc(this->i_reds, sizeof(INT));
    MEM_ERROR(this->ti_reds);

    this->ti_weights              = (INT*)calloc(this->i_vertices*this->i_weights, sizeof(INT));
    MEM_ERROR(this->ti_weights);

    return (0);
}

int rbhFree(Hypergraph *  this){

    free(this->ti_weights);
    free(this->ti_reds);
    free(this->ti_criticalities_right);
    free(this->ti_criticalities_left);
    free(this->ti_delays);
    free(this->ti_idx_hyperedges);
    free(this->ti_hyperedges);
    if (this->s_rbh_name!=NULL) {
        free(this->s_rbh_name);
    }


    return (0);
}

int varRbhLoad(rbhLoad_args in){

    Hypergraph * h_out             = in.h ? in.h                 : NULL;
    char *     s_path_out          = in.s_path ? in.s_path       : NULL;
    INT        i_baseval_out       = in.i_baseval ? in.i_baseval : 0;
    bool       b_verbose_out       = in.b_verbose ? in.b_verbose : false;

    return rbhLoadBase(h_out, s_path_out, i_baseval_out, b_verbose_out);
}

int rbhLoadBase(Hypergraph *  this, const char * const s_path, INT i_baseval, bool b_verbose){

    /*
    ** Errors.
    */
    if (this == NULL) {
        return (1);
    }
    if (s_path == NULL) {
        return (2);
    }
    if (i_baseval < 0) {
        return (3);
    }
    /*
    ** End errors.
    */

    char * buffer = (char*)malloc(BUFSIZE);
    MEM_ERROR(buffer);

    FILE * in = NULL;
    in        = fopen(s_path, "r");

    if (in == NULL) {
        fprintf(stderr, "Cannot open file %s\n", s_path);
        exit(0);
    }

    if (b_verbose) {
        printf("File : %s\n", s_path);
    }

    INT i_pins, i_vertices, i_hyperedges, i_reds, i_weights;
    readLine(in, buffer, BUFSIZE);

    int n = sscanf(buffer, "%d %d %d %d %d", &i_pins, &i_vertices, &i_hyperedges, &i_reds, &i_weights);
    if (n<0) {
        return (10);
    }
    FATAL(n!=5, "Invalid graph format!");
    assert(n==5);

    if (b_verbose) {
        printf("#pins : %d\n#vertex : %d\n#edges : %d\n#red : %d\n", i_pins, i_vertices, i_hyperedges, i_reds);
    }

    this->i_vertices       = i_vertices;
    this->i_hyperedges     = i_hyperedges;
    this->i_reds           = i_reds;
    this->i_pins           = i_pins;
    this->i_weights        = i_weights;

    rbhInit(this);

    char   token[2];
    strcpy(token, " ");
    char * raw;
    char * buffer2      = (char*)malloc(sizeof(char)*BUFSIZE);
    MEM_ERROR(buffer2);

    INT  idx            = 0;
    INT  len;
    bool is_in;
    INT  cur_len;
    INT  idx_e;
    for(INT j = 0; j < this->i_hyperedges; j++) {

        this->ti_idx_hyperedges[j] = idx;
        idx_e                      = idx;

        readLine(in, buffer, BUFSIZE);

        if (b_verbose) {
            printf("Hyperedges : j(%d)\n%s\n", j, buffer);
        }

        strcpy(buffer2, buffer);

        raw = strtok(buffer, token);
        len = 0;

        while (raw != NULL)
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
        for (INT i = 0; i < len - 1; i++) {
            is_in                      = false;
            this->ti_hyperedges[idx++] = atoi(raw) - i_baseval;  /* vertex */
            for (INT ii = 0; ii < cur_len; ii++) {
                if (this->ti_hyperedges[idx_e+2+ii] == this->ti_hyperedges[idx-1]) {
                    is_in = true;
                }
            }
            cur_len++;
            if (is_in) {

                cur_len--;
                idx--;
            }
            raw = strtok(NULL, token);
        }
        this->ti_hyperedges[idx_e+1] = cur_len;
    }

    idx = 0;

    for(INT i = 0; i < this->i_vertices; i++) {

        /*
        ** red (0/1) delay criticality weights
        */

        readLine(in, buffer, BUFSIZE);
        if(b_verbose) {
            printf("vertex : %d\t %s\n", i, buffer);
        }
        raw = strtok(buffer, token);
        if(atoi(raw)==1) {
            this->ti_reds[idx++] = i;
        }
        raw = strtok(NULL, token);
        this->ti_delays[i]              = atoi(raw);
        raw = strtok(NULL, token);
        this->ti_criticalities_right[i] = atoi(raw);

        for(INT iw = 0; iw < this->i_weights; iw ++) {
            raw = strtok(NULL, token);
            this->ti_weights[i*i_weights+iw]      = atoi(raw);
        }

    }


    free  (buffer2);
    free  (buffer);
    fclose(in);

    return (0);

}

int rbhSave(Hypergraph *  this, const char * const s_path, INT i_baseval, bool b_verbose){
    /*
    ** Errors.
    */
    if (this == NULL) {
        return (1);
    }
    if (s_path == NULL) {
        return (2);
    }
    if (i_baseval < 0) {
        return (3);
    }
    /*
    ** End errors.
    */

    FILE * out = NULL;
    out        = fopen(s_path, "w+");

    if (out == NULL) {
        fprintf(stderr, "Cannot open file %s\n", s_path);
        exit(0);
    }

    if (b_verbose) {
        printf("File : %s\n", s_path);
    }

    int n = fprintf(out, "%d %d %d %d %d\n", this->i_pins, this->i_vertices, this->i_hyperedges, this->i_reds, this->i_weights);


    if (n<0) {
        return (10);
    }


    if (b_verbose) {
        printf("#pins : %d\n#vertex : %d\n#edges : %d\n#red : %d\n", this->i_pins, this->i_vertices, this->i_hyperedges, this->i_reds);
    }

    bool * is_red = (bool*)malloc(sizeof(bool)*this->i_vertices);
    MEM_ERROR(is_red);

    for (INT i = 0; i < this->i_vertices; i++) {
        is_red[i] = false;
    }
    for (INT i = 0; i < this->i_reds; i++) {
        is_red[this->ti_reds[i]] = true;
    }

    INT idx, len;

    for(INT j = 0; j < this->i_hyperedges; j++) {

        idx = this->ti_idx_hyperedges[j];

        len = this->ti_hyperedges[idx+1];

        fprintf(out, "%d", this->ti_hyperedges[idx++]);
        idx++;
        for (INT i = 0; i < len; i++) {
            fprintf(out, " %d", this->ti_hyperedges[idx++]);
        }
        fprintf(out, "\n");
    }

    idx = 0;

    for(INT i = 0; i < this->i_vertices; i++) {

        /*
        ** red (0/1) delay criticality weights
        */
        if (is_red[i]) {
            fprintf(out, "1");
        }else {
            fprintf(out, "0");
        }

        fprintf(out, " %d", this->ti_delays[i]);
        fprintf(out, " %d", this->ti_criticalities_right[i]);
        for(INT iw = 0; iw < this->i_weights; iw ++) {
            fprintf(out, " %d", this->ti_weights[i*this->i_weights+iw]);
        }
        fprintf(out, "\n");
    }


    free(is_red);
    fclose(out);

    return (0);

}

int computeListNeighbors(Hypergraph * this, List ** neighbors_list) {

    bool is_in;
    INT   nv     = this->i_vertices;
    INT   ne     = this->i_hyperedges;
    INT * e      = this->ti_hyperedges;
    INT * idx_e  = this->ti_idx_hyperedges;

    for (INT j = 0; j < ne; j++) {
        INT idx_j  = idx_e[j];
        INT size_j = e[idx_j+1];
        INT u      = e[idx_j+2];

        for(INT i = 0; i < size_j-1; i++) {
            is_in = false;
            INT v     = e[idx_j+3+i];
            List * cell = neighbors_list[u];
            while(cell != NULL) {
                if (cell->i == v) {
                    is_in = true;
                }
                cell = cell->next;
            }
            if (!is_in) {
                list_add_element(neighbors_list[u], v);
            }
        }
    }

    return (0);
}

int computeListInNeighbors(Hypergraph * this, List ** neighbors_list, List ** in_neighbors_list) {

    INT   nv     = this->i_vertices;
    bool  is_in;
    for(INT v = 0; v < nv; v++) {
        List * cell = neighbors_list[v];
        while(cell != NULL && cell->size>0) {
            INT u = cell->i;
            is_in = false;
            List * cell_in = in_neighbors_list[u];
            while(cell_in != NULL) {
                if (cell_in->i == v) {
                    is_in = true;
                }
                cell_in = cell_in->next;
            }
            if (!is_in) {
                list_add_element(in_neighbors_list[u], v);
            }
            cell = cell->next;
        }
    }
    return (0);
}

int computeNeighbors(Hypergraph * this, Matrix * neighbors) {

    INT   nv     = this->i_vertices;
    INT   ne     = this->i_hyperedges;
    INT * e      = this->ti_hyperedges;
    INT * idx_e  = this->ti_idx_hyperedges;

    new_matrix(neighbors, nv, nv);

    for(INT i = 0; i < nv; i++) {
        neighbors->v[i].size = 0;
    }
    INT u, v;
    bool is_in;
    for(INT j = 0; j < ne; j++) {
        INT idx_j  = idx_e[j];
        INT size_j = e[idx_j+1];
        u          = e[idx_j+2];
        for(INT i = 0; i < size_j-1; ++i) {
            is_in = false;
            v     = e[idx_j+3+i];
            if (neighbors->v[u].size==0) {
                neighbors->v[u].v[neighbors->v[u].size++] = v;

            }else {
                for (INT x = 0; x < neighbors->v[u].size; x++) {

                    if (neighbors->v[u].v[x] == v) {
                        is_in = true;
                    }
                }
                if (!is_in) {
                    if (neighbors->v[u].size == neighbors->n) {
                        unsigned long n = neighbors->n;
                        for(INT x = 0; x < neighbors->m; x++) {
                            void * err = realloc(neighbors->v[x].v, sizeof(long)*2*n);
                            if (err == NULL) {
                                return (2);
                            }
                        }
                    }
                    neighbors->v[u].v[neighbors->v[u].size++] = v;

                }
            }
        }
    }
    return (0);
}

int computeInNeighbors(Hypergraph * this, Matrix * neighbors, Matrix * in_neighbors) {

    INT   nv     = this->i_vertices;

    new_matrix(in_neighbors, nv, nv);
    for(INT v = 0; v < nv; v++) {
        in_neighbors->v[v].size=0;
    }
    bool is_in;
    for(INT v = 0; v < nv; v++) {
        for (INT i = 0; i<neighbors->v[v].size; i++) {
            is_in = false;
            INT u = neighbors->v[v].v[i];
            for(INT x=0; x < in_neighbors->v[u].size; x++) {
                if (in_neighbors->v[u].v[x] == v) {
                    is_in=true;
                }
            }
            if(!is_in) {
                in_neighbors->v[u].v[in_neighbors->v[u].size++] = v;
            }
        }
    }
    return (0);
}



int topologicalSort(Hypergraph * this, List ** neighbors, List ** in_neighbors, INT * sort) {

    INT nv     = this->i_vertices;
    INT ne     = this->i_hyperedges;
    INT nr     = this->i_reds;
    INT np     = this->i_pins;
    INT nwv    = this->i_weights;
    INT * reds = this->ti_reds;

    INT i,u,v;

    INT * queue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(queue);

    INT * main_queue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(main_queue);

    INT * indeg = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(indeg);
    bool * flag      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    bool * is_red      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);

    for(i = 0; i < nv; i++) {
        indeg[i]      = in_neighbors[i]->size;
        flag[i]       = false;
        is_red[i]     = false;
        queue[i]      =-1;
        main_queue[i] =-1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
    }

    INT main_queue_start = 0;
    INT main_queue_end   = 0;
    for(u = 0; u < nv; u++) {
        if (is_red[u]) {
            main_queue[main_queue_end++] = u;
            flag[u]                      = true;
        }
    }
    INT index = 0, queue_start = 0, queue_end = 0;

    while (main_queue_start != main_queue_end) {
        queue_start=0;
        queue_end=0;
        for (i = 0; i < main_queue_end; i++) {
            // all red vertices visited are explored next
            queue[queue_end++] = main_queue[i];
        }
        main_queue_start = 0;
        main_queue_end   = 0;
        while(queue_start != queue_end) {
            u             = queue[queue_start++];
            sort[index++] = u;
            flag[u]       = true;
            List * cell   = neighbors[u];
            for(INT x = 0; x < neighbors[u]->size; x++) {
            
                v = cell->i;
                indeg[v]--;
                if (indeg[v] <= 0 && !is_red[v] && !flag[v]) {
                    flag[v]            = true;
                    queue[queue_end++] = v;
                }
                if (is_red[v] && !flag[v]) {
                    main_queue[main_queue_end++] = v;
                    flag[v]                      = true;
                }
                cell = cell->next;
            }
        }
    }

    if (index < nv) {
        // if reds vertex u,v st.
        //  u --> v --> u

        main_queue_start = 0;
        main_queue_end   = 0;
        for(u = 0; u < nv; u++) {
            if (is_red[u] && !flag[u]) {
                main_queue[main_queue_end++] = u;
                flag[u]                      = true;
            }
            if (indeg[u]==0 && !flag[u]) {
                main_queue[main_queue_end++] = u;
                flag[u]                      = true;
            }
        }
        /* if exists interdependance between u --> v --> u with v black */
        /* !flag[v] = true but indeg[v]>0  */
        for(u = 0; u < nv; u++) {
            if(flag[u] && is_red[u]) {

                List * cell   = neighbors[u];
                 for(INT x = 0; x < neighbors[u]->size; x++) {
                    v = cell->i;
                    if(!flag[v]){
                        main_queue[main_queue_end++] = v;
                        flag[v]                      = true;
                    }
                    cell = cell->next;
                }
            }

        }


        while (main_queue_start != main_queue_end) {
            queue_start=0;
            queue_end=0;
            for (i = 0; i < main_queue_end; i++) {
                // all red vertices visited are explored next
                queue[queue_end++] = main_queue[i];
            }
            main_queue_start = 0;
            main_queue_end   = 0;
            while(queue_start != queue_end) {
                u             = queue[queue_start++];
                sort[index++] = u;
                flag[u]       = true;
                List * cell   = neighbors[u];
                for(INT x = 0; x < neighbors[u]->size; x++) {
                    v = cell->i;
                    
                    if (!is_red[v] && !flag[v]) {
                        flag[v]            = true;
                        queue[queue_end++] = v;
                    }
                    if (is_red[v] && !flag[v]) {
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
    for(INT i = 0; i < nv; i++) {
    
        if(!flag[i]) {
        
           sort[index++] = i;
           
          
        
        }
    
    }


    
    
    assert(index==nv);
    free(queue);
    free(main_queue);
    free(indeg);
    free(flag);
    free(is_red);

    return(0);

}

int computeListNeighborsUnalloc(Hypergraph * this, List ** neighbors_list) {

    bool is_in;
    INT   nv     = this->i_vertices;
    INT   ne     = this->i_hyperedges;
    INT * e      = this->ti_hyperedges;
    INT * idx_e  = this->ti_idx_hyperedges;
    
    for(INT i = 0; i < nv; i++) {
        List * cell = neighbors_list[i];
        cell->size = 0;
        cell->i = 0;
    }
    
    List * cell;
    INT  idx_j, size_j, u, v, neighbors_size;
    for (INT j = 0; j < ne; j++) {
        idx_j  = idx_e[j];
        size_j = e[idx_j+1];
        u      = e[idx_j+2];

        for(INT i = 0; i < size_j-1; i++) {
            is_in = false;
            v     = e[idx_j+3+i];
            cell      = neighbors_list[u];
            
            neighbors_size = neighbors_list[u]->size;
            
            for(INT x = 0; x < neighbors_size; x++) {
                
                if (cell->i == v) {
                
                    is_in = true;
                
                }
                cell = cell->next;
            
            }
            if (!is_in) {
                
                if (cell != NULL) { 
                    
                    cell->i = v;
                    neighbors_list[u]->size++;
                    cell->size = neighbors_list[u]->size;
                    
                } else {
                
                    list_add_element(neighbors_list[u], v);
                    
                
                }
            }
        }
    }

    return (0);
}

int computeListInNeighborsUnalloc(Hypergraph * this, List ** neighbors_list, List ** in_neighbors_list) {

    INT   nv     = this->i_vertices;
    bool  is_in;
    
    
    for(INT i = 0; i < nv; i++) {
        
        in_neighbors_list[i]->size = 0;
    }
    
    INT u, v, neighbors_size, in_neighbors_size;
    
    for(v = 0; v < nv; v++) {
    
        List * cell = neighbors_list[v];
        
        neighbors_size =  neighbors_list[v]->size;
        
        for(INT i = 0; i < neighbors_size; i++) {
        
            u = cell->i;
            is_in = false;
            
            List * cell_in = in_neighbors_list[u];
            
            in_neighbors_size = in_neighbors_list[u]->size;
            
            for(INT x = 0; x < in_neighbors_size; x++) {
                
                if (cell_in->i == v) {
                
                    is_in = true;
                
                }
                cell_in = cell_in->next;
            
            }
            if (!is_in) {
                if (cell_in != NULL) { 
                    
                    cell_in->i = v;
                    in_neighbors_list[u]->size++;
                    cell_in->size = in_neighbors_list[u]->size;
                    
                } else {
                
                    list_add_element(in_neighbors_list[u], v);
                    
                }
            }
            
            cell = cell->next;
        }
    }
    return (0);
}

int computeInNeighborsUnalloc(
Hypergraph * h, 
Matrix     * neighbors, 
Matrix     * in_neighbors) {
  
  
  INT   nv = h->i_vertices;
  INT   ne = h->i_hyperedges;
  
  in_neighbors->m = nv;
  in_neighbors->n = nv;

  for(INT v = 0; v < h->i_vertices; v++) {
      
      in_neighbors->v[v].size=0;
  
  }

  bool is_in;
  for(INT v = 0; v < h->i_vertices; v++) {

    for (INT i = 0; i< neighbors->v[v].size; i++) {
       is_in = false;
       INT u = neighbors->v[v].v[i];
       for(INT x = 0; x < in_neighbors->v[u].size; x++) {
         if (in_neighbors->v[u].v[x] == v) {
             is_in=true;
         }
       }
       if(!is_in)
       {
         in_neighbors->v[u].v[in_neighbors->v[u].size++] = v;
       }

    }
  }

  return 0;
}

int computeNeighborsUnalloc(Hypergraph * h, Matrix * neighbors) {
  INT   nv = h->i_vertices;
  INT   ne = h->i_hyperedges;
  
  INT * idx_hyperedges = h->ti_idx_hyperedges;
  INT * hyperedges     = h->ti_hyperedges;

  neighbors->m=nv;
  neighbors->n=nv;
  for(INT i = 0; i < nv; i++) {
    neighbors->v[i].size = 0;
  }

  INT u, v, idx_j, weight, size_j;
  bool is_in;
  for(INT j = 0; j < ne; j++){
      idx_j = idx_hyperedges[j];
      
      weight = hyperedges[idx_j];
      size_j = hyperedges[idx_j+1];
      
      u = hyperedges[idx_j+2];
      for(INT i = 1; i < size_j; ++i){
          is_in = false;
          v     = hyperedges[idx_j+2+i];
          
          if (neighbors->v[u].size==0) {
              neighbors->v[u].v[neighbors->v[u].size++] = v;
          } else {
              for (INT k = 0; k < neighbors->v[u].size; k++) {
                 if (neighbors->v[u].v[k] == v) {
                     is_in = true;
                 }
              }
              if (!is_in) {
                  if (neighbors->v[u].size == neighbors->n) {
                      printf("realloc %lu , %lu\n", neighbors->v[u].size, neighbors->n);
                      //realloc pour tous
                      unsigned long n = neighbors->n;
                      for(unsigned int x=0; x < neighbors->m; x++) {
                          void * res = realloc(neighbors->v[x].v, sizeof(long)*2*n);
                          MEM_ERROR(res);
                      }
                      neighbors->n *= 2;
                  }
                  neighbors->v[u].v[neighbors->v[u].size++] = v;
              }
          }
      }
  }


  return 0;
}


int compute_criticality(Hypergraph * h, List ** neighbors, List **  in_neighbors, INT * sort) {
    bool verbose = false;

    INT nvp    = 0;
    INT nep    = 0;
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    bool * is_red = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);
    bool * flag = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    INT * delays_fwd = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(delays_fwd);
    INT * delays_bwd = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(delays_bwd);

    INT crit_max = 0;
    for(INT i = 0; i < nv; i++) {
        delays_fwd[i] = 0;
        delays_bwd[i] = 0;
        flag[i]       = false;
        is_red[i]     = false;
    }
    for (INT i = 0; i < nr; i++) {
        is_red[reds[i]] = true;
    }

    INT u,v;
    for(INT i = 0; i < nv; i++) {
        v       = sort[i];
        flag[v] = true;
        if (is_red[v]) {
            List * cell = neighbors[v];
            while( cell != NULL && neighbors[v]->size) {

                u = cell->i;
                if ( !is_red[u] ) {
                    delays_fwd[u] = MAX(delays_fwd[u], h->ti_delays[v]+h->ti_delays[u]);
                }
                cell = cell->next;
            }
        }
        else {
            List * cell = neighbors[v];
            while( cell != NULL && neighbors[v]->size) {

                u = cell->i;
                if ( !is_red[u] ) {
                    delays_fwd[u] = MAX(delays_fwd[u], delays_fwd[v]+h->ti_delays[u]);
                }
                cell = cell->next;
            }
        }
    }

    INT dhat;

    // retro propagation
    for (INT i = nv-1; i > 0; i--) {
        v       = sort[i];
        dhat    = 0;
        if (is_red[v] || delays_bwd[v] ==0 ){
            delays_bwd[v] = delays_fwd[v];
        }
        List * cell = in_neighbors[v];
        while( cell != NULL && in_neighbors[v]->size) {

            u = cell->i;
            if (dhat < delays_fwd[u]) {
                dhat = delays_fwd[u];
            }
            cell = cell->next;
        }
        cell = in_neighbors[v];
        while( cell != NULL && in_neighbors[v]->size) {

            u = cell->i;
            if ( !is_red[u] ) {
                delays_bwd[u] = MAX(delays_bwd[u], delays_bwd[v]-(dhat-delays_fwd[u]));
            }
            cell = cell->next;
        }
    }

    for( v = 0; v < nv; v++) {
        h->ti_criticalities_right[v] = MAX(delays_bwd[v], 1);
        h->ti_criticalities_left[v]  = MAX(delays_fwd[v], 1);
    }

    for(INT j = 0; j < ne; j++) {
        INT idx_e  = h->ti_idx_hyperedges[j];
        INT source = h->ti_hyperedges[idx_e+2];
        h->ti_hyperedges[idx_e] = h->ti_criticalities_right[source];
    }

    free(delays_fwd);
    free(delays_bwd);
    free(is_red);
    free(flag);

    return 0;
}

int compute_maxcon(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmaxcon, INT * lmaxcon, float * avg, float * stdw) {
  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;


  INT max_con=0;
  
  (*avg) = 0;
  
  float var = 0;
  
  INT * e      = h->ti_hyperedges;
  INT * idx_e  = h->ti_idx_hyperedges;

  for (INT j = 0; j < ne; j++) {
     INT idx_j  = idx_e[j];
     INT size_j = e[idx_j+1];
     (*avg) += size_j;
     var += size_j*size_j;
     
     if(max_con < size_j) {
         max_con = size_j;
     }
        
  }
  
  (*avg) /= (float)nv;
  (*stdw) = sqrt(var/(float)nv - (*avg)*(*avg)); 
   
   
   for (INT j = 0; j < ne; j++) {
     INT idx_j  = idx_e[j];
     INT size_j = e[idx_j+1];
     
     
     if(max_con == size_j) {
         lmaxcon[(*sizelmaxcon)++] = j;
     }
        
  }


  return max_con;

}

int compute_maxdeg(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmaxdeg, INT * lmaxdeg, float * avg, float * stdw) {
  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;
  INT  * degrees = (INT*)malloc(sizeof(INT)*nv);
  MEM_ERROR(degrees);

  INT max_deg=0;
  
  (*avg) = 0;
  
  float var = 0;

  for(INT   i=0; i<nv; i++)
  {
    degrees[i] = neighbors[i]->size;
    degrees[i] += in_neighbors[i]->size;
    
    (*avg) += degrees[i];
    var += degrees[i]*degrees[i];
    
    if(degrees[i]>max_deg) {
        max_deg = degrees[i];
    }

  }
  
  (*avg) /= (float)nv;
  (*stdw) = sqrt(var/(float)nv - (*avg)*(*avg)); 

  
  for(INT v=0; v<nv; v++)
  {

		if (degrees[v]==max_deg){
			lmaxdeg[(*sizelmaxdeg)++] = v;
		}
  }
  

  free(degrees);

  return max_deg;

}

int compute_pmax(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelpmax, INT * lpmax,  float * avg, float * stdw)  {
  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;
  INT  * delays = (INT*)malloc(sizeof(INT)*nv);
  MEM_ERROR(delays);
  bool * is_red   = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(is_red);
  bool * flag = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(flag);


  INT   u, v;

  INT ncuts = 0;
  INT idx_pmax = 0;

  (*avg) = 0;


  for(INT   i=0; i< nv; i++)
  {
    delays[i] =  h->ti_delays[i];
    flag[i]   = false;
    is_red[i]    = false;
  }

  for (INT i = 0; i < nr; i++) {
        is_red[reds[i]] = true;
      }

  for(INT   i=0; i<nv; i++)
  {
    v       = sort[i];
    flag[v] = true;

    if (is_red[v])
    {
      List * cell_neighbors = neighbors[v];
            while(cell_neighbors != NULL && neighbors[v]->size>0) {
                        //selection of vertex v
                        INT u = cell_neighbors->i;
                       
                                delays[u] = MAX(delays[u], h->ti_delays[v]+h->ti_delays[u]);

                        cell_neighbors = cell_neighbors->next;
                    }
     }
    else
    {
      List * cell_neighbors = neighbors[v];
      while(cell_neighbors != NULL && neighbors[v]->size>0) {
                  //selection of vertex v
                  INT u = cell_neighbors->i;

                          delays[u] = MAX(delays[u], delays[v]+h->ti_delays[u]);
                  cell_neighbors = cell_neighbors->next;
              }

    }
  }

  int max_v = -1;
  int max_delays = 0;
  float var = 0;
  for(INT v=0; v<nv; v++)
  {
  
      (*avg) += delays[v];
      var        += delays[v] * delays[v];
		if (delays[v]>max_delays){
			max_delays=delays[v];
			max_v = v;
		}
  }
  (*avg) /= (float)nv;
  (*stdw) = sqrt(var/(float)nv - (*avg)*(*avg));
  
  for(INT v=0; v<nv; v++)
  {

		if (delays[v]==max_delays){
			lpmax[idx_pmax++] = v;
		}
  }
  
  (*sizelpmax) = idx_pmax;

  free(delays);
  free(flag);
  free(is_red);

  return max_delays;

}

int compute_path_length(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, INT * sizelmax, INT * lmax, INT * depth, float * avg, float * stdw) {
  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;
  bool * is_red   = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(is_red);
  
  
  
  
  INT start;
  INT end;

  INT   u, v;

  INT ncuts = 0;
  INT idx_lmax = 0;
  
  (*avg) = 0;



  for(INT   i=0; i< nv; i++)
  {
    depth[i] = 0;
    is_red[i]    = false;
  }

  for (INT i = 0; i < nr; i++) {
        is_red[reds[i]] = true;
      }
  
  for(INT   i=0; i<nv; i++)
  {
    v       = sort[i];

    if (is_red[v])
    {
      List * cell_neighbors = neighbors[v];
            while(cell_neighbors != NULL && neighbors[v]->size>0) {
                        //selection of vertex v
                        INT u = cell_neighbors->i;
                       
                                depth[u] = MAX(depth[u], 2);

                        cell_neighbors = cell_neighbors->next;
                    }
     }
    else
    {
      List * cell_neighbors = neighbors[v];
      while(cell_neighbors != NULL && neighbors[v]->size>0) {
                  //selection of vertex v
                  INT u = cell_neighbors->i;

                          depth[u] = MAX(depth[u], depth[v]+1);
                  cell_neighbors = cell_neighbors->next;
              }

    }
  }
      

  INT max_length = 0;
  float var = 0;
  for(INT v=0; v<nv; v++)
  {
       (*avg) += depth[v];
       var += depth[v] * depth[v];
		if (depth[v]>max_length){
			max_length = depth[v] ;
			
		}
  }
  (*avg) /= (float)nv;
  (*stdw) = sqrt(var/(float)nv - (*avg)*(*avg));
  for(INT v=0; v<nv; v++)
  {

		if (depth[v]==max_length){
			lmax[idx_lmax++] = v;
			
		}
  }
  (*sizelmax) = idx_lmax;
  
  free(is_red);
  return max_length;

}

int compute_subpmax(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * sort, bool * vertices) {
  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;
  INT  * delays = (INT*)malloc(sizeof(INT)*nv);
  MEM_ERROR(delays);
  bool * is_red   = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(is_red);
  bool * flag = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(flag);


  INT   u, v;

  INT ncuts = 0;




  for(INT   i=0; i< nv; i++)
  {
    delays[i] = 0;
    if (vertices[i]){
        delays[i] =  h->ti_delays[i];
    }
    flag[i]   = false;
    is_red[i]    = false;
  }

  for (INT i = 0; i < nr; i++) {
        is_red[reds[i]] = true;
      }

  for(INT   i=0; i<nv; i++)
  {
    v       = sort[i];
    flag[v] = true;

    if (is_red[v])
    {
      List * cell_neighbors = neighbors[v];
            while(cell_neighbors != NULL && neighbors[v]->size>0) {
                        //selection of vertex v
                        INT u = cell_neighbors->i;
                        if(vertices[u] && vertices[v]){
                                delays[u] = MAX(delays[u], h->ti_delays[v]+h->ti_delays[u]);
                              }
                        cell_neighbors = cell_neighbors->next;
                    }
     }
    else
    {
      List * cell_neighbors = neighbors[v];
      while(cell_neighbors != NULL && neighbors[v]->size>0) {
                  //selection of vertex v
                  INT u = cell_neighbors->i;
                  if(vertices[u] && vertices[v]){

                          delays[u] = MAX(delays[u], delays[v]+h->ti_delays[u]);
                        }
                  cell_neighbors = cell_neighbors->next;
              }

    }
  }

  int max_v = -1;
  int max_delays = 0;
  for(INT v=0; v<nv; v++)
  {

		if (delays[v]>max_delays){
			max_delays=delays[v];
			max_v = v;
		}
  }


  free(delays);
  free(flag);
  free(is_red);

  return max_delays;

}
