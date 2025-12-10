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
/**   NAME       : fm.c                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are definitions for the     **/
/**                FM base refinement functions.           **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 01 jan 2022     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "fm.h"

INT dkfm(
Hypergraph * h,  
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance,
INT k) {

    bool verbose = false;
    bool unsize  = true;
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    /** Compute the halo of the partition **/
    
    bool  * is_in_halo = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_in_halo);
    bool * is_locked   = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_locked);
    INT halo_size      = 0;
    
    INT * part_sizes   = (INT*)calloc(k*nwv, sizeof(INT));
    MEM_ERROR(part_sizes);
    
    INT * capacity_cst  = (INT*)calloc(nwv, sizeof(INT));
    MEM_ERROR(capacity_cst);
    
    float * cost = (float*) malloc(sizeof(float)*nv);
    MEM_ERROR(cost);
    
    for (INT u = 0; u < nv; u++) {        
        
        /* compute the part_sizes */
        for(INT wi = 0; wi < nwv; wi++) {
           part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
           capacity_cst[wi] += h->ti_weights[u*nwv+wi];
        }
    
        
        is_in_halo[u] = false;
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
                    INT v = cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo[u] = true;
            }
            
            cell = cell->next;
        }
        List * in_cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
        
            INT v = in_cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo[u] = true;
            }
            
            in_cell = in_cell->next;
        }
        if (is_in_halo[u]) {
            halo_size++;
        }
        is_locked[u] = false;
    }
    for(INT wi = 0; wi < nwv; wi++) {
        
        capacity_cst[wi] = ceil((float)capacity_cst[wi]/(float)k *1.05);
    }
    
    
    /** compute the moves arrays for each part **/
    
    //compute minimum and maximum gain
    INT * partitionp = (INT*)calloc(nv, sizeof(INT));
    MEM_ERROR(partitionp);
        
    INT min_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    
    for (INT u = 0; u < nv; u++) {
        partitionp[u] = u;
    }
    
    INT max_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    
    INT cur_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
    
    INT min_gain = cur_crit_path - max_crit_path; // \in [-max_crit_path, 0]
    INT max_gain = cur_crit_path - min_crit_path; // \in [0, max_crit_path]
    
    for (INT u = 0; u < nv; u++) {
        partitionp[u] = partition[u];
        
        INT D = 0;
        
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
        
            INT v = cell->i;
            
            if (a->ti_delay[partition[u]*a->i_m+partition[v]] > D) {
            
                D = a->ti_delay[partition[u]*a->i_m+partition[v]];
            }
            
            cell = cell->next;
            
        }
        
    
        /* compute cost */
        cost[u] = ((float)tolerance/100.0) * ( ((float)max_crit_path - (float)h->ti_criticalities_right[u]) / (float)D);
    }
    
    /* relative to 0 */
    
    INT gain_size = max_gain - min_gain;
    INT relative  = gain_size - min_gain;
    
    /* gain data structure declaration */
        
    INT ** partition_moves  = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_moves);
    for (INT p = 0; p < k; p++) {
        partition_moves[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_moves[p]);
    }
    
    INT ** partition_gains = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_gains);
    for (INT p = 0; p < k; p++) {
        partition_gains[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_gains[p]);
        for (INT i = 0; i < nv; i++) {
        
            partition_gains[p][i] = min_gain;
        }
        
    }    
    INT * partition_gain_sizes = (INT*)calloc(k, sizeof(INT));
    MEM_ERROR(partition_gain_sizes);
    
    bool * is_computed = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(is_computed);
    
    INT ncomp;
    INT anc_part;
    INT gain;
    List * cell, * in_cell, * cellp;
    
    INT u, v, vp;
    
    for(INT i = 0; i < nv; i++) {
        
        if (is_in_halo[i]) {       
        
            for(INT p = 0; p < k; p++) {
                 is_computed[p] = false;
            }
            
            ncomp = 0;
            
                	    
    	    cell   = neighbors[i];
    	    for(INT x = 0; x < neighbors[i]->size && ncomp<k; x++) {
            
                v = cell->i;
                if (partition[i] != partition[v]) {
                    if(!is_computed[partition[v]]) {
                        /* add u to list of gain of partition[v] */
                        anc_part = partition[i];
                        partition[i] = partition[v];
                        
                        gain = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                        
                        // add the vertices corresponding to its gain 
                        
                        partition_gains[partition[i]][i] = gain;
                        
                                                
                        pqueue_add_element(partition_moves[partition[i]], partition_gains[partition[i]], i, partition_gain_sizes[partition[i]]);
                        partition_gain_sizes[partition[i]]++;
                                       
                        partition[i] = anc_part;
                        
                        is_computed[partition[v]] = true;
                        ncomp++;
                    }
                }
                cell = cell->next;
            }
            in_cell   = in_neighbors[i];
            for(INT x = 0; x < in_neighbors[i]->size && ncomp<k; x++) {
            
                v = in_cell->i;
                if (partition[i] != partition[v]) {
                    if(!is_computed[partition[v]]) {
                        /* add u to list of gain of partition[v] */
                        anc_part = partition[i];
                        partition[i] = partition[v];
                        
                        gain = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                       
                        // add the vertices corresponding to its gain 
                        partition_gains[partition[v]][i] = gain;
                        
                        pqueue_add_element(partition_moves[partition[v]], partition_gains[partition[v]], i, partition_gain_sizes[partition[v]]);
                        partition_gain_sizes[partition[v]]++;
                                       
                        partition[i] = anc_part;
                        
                        is_computed[partition[v]] = true;
                        ncomp++;
                     }
                }
                in_cell = in_cell->next;
            }
        
        }    	
    
    }
    
    bool * to_balance = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(to_balance);
    
    for(INT p = 0; p < k; p++) {
        
        to_balance[p] = false;
        
        for(INT wi = 0; wi < nwv; wi++){
            if(part_sizes[p] >= capacity_cst[wi]){
                
                to_balance[p] = true;
            
            }   
        }
    }
    
    /* DKFM moves */
    
    INT base_cur_crit_path = cur_crit_path;
    INT crit_path_up = 0;
    INT best_crit_path = cur_crit_path;
    
    bool exact = false;
    
    INT * ex_partition = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(ex_partition);
    
    for(INT i = 0; i < nv; i++) {
       ex_partition[i] = partition[i];
    }
    
    INT N_MOVES = ceil((float)nv*(1.0/(float)perform));
    for(INT i = 0; i < N_MOVES; i++) {
        
               
    
    	/* select part */
    	INT selected_part = -1;
    	
    	for(INT p = 0; p < k; p++) {
            
            if(!to_balance[p] && partition_gain_sizes[p]>0) {
                
                selected_part = p;
                
            }
            
         }
         
         if (selected_part == -1) {
             break;
             break;
         }
         
         
            
         INT u = pqueue_dequeue(partition_moves[selected_part], partition_gains[selected_part], partition_gain_sizes[selected_part]);
         partition_gain_sizes[selected_part]--;
         is_locked[u] = true;
         cur_crit_path = base_cur_crit_path - partition_gains[selected_part][u];
            
         for(INT wi = 0; wi < nwv; wi++){
             part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
             part_sizes[partition[u]*nwv+wi]  -= h->ti_weights[u*nwv+wi];
             if(part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]){
                
                to_balance[selected_part] = true;
            
             }   
            
          }
            
            
          if (to_balance[partition[u]]) {
              to_balance[partition[u]] = false;
          }
          
          for(INT wi = 0; wi < nwv; wi++){

             if(part_sizes[partition[u]*nwv+wi] <= capacity_cst[wi]){
                
                to_balance[partition[u]] = false;
            
             }else {
                to_balance[partition[u]] = true;
             }   
            
          }
            
          /* update the current gain */ 
          
          cell = neighbors[u];
          for(INT x = 0; x < neighbors[u]->size; x++) {
                          
             v = cell->i;
                
             if (!is_locked[v]) {
                bool is_connected = false;
                cellp = neighbors[v];
                for(INT o = 0; o < neighbors[v]->size; o++) {
                                
                    vp = cellp->i;
                    
                    if (partition[vp] == partition[u]) {
                           
                        // the future anciant part of u
                        is_connected = true;
                       
                    }
                    cellp = cellp->next;
                }
                if (!is_connected) {
                   /* delete v of anciant part of u moves */
                   partition_gains[partition[u]][v] = min_gain;
                   pqueue_siftUp(partition_moves[partition[u]], partition_gains[partition[u]], v, partition_gain_sizes[partition[u]]);

                }
                if (partition[v] != selected_part) {
                       
                    anc_part = partition[v];
                       
                    partition[u] = selected_part;
                    partition[v] = selected_part;
                       
                    partition_gains[selected_part][v] = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                                        
                    pqueue_siftUp(partition_moves[selected_part], partition_gains[selected_part], v, partition_gain_sizes[selected_part]);
                    
                     
                                
                    partition[v] = anc_part;
                       
                 }
              }
              
              cell = cell->next;
        
       } 
       partition[u] = selected_part;
       crit_path_up++;
       
       cur_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
       exact = true;
       crit_path_up = 0;
       
       if (cur_crit_path < best_crit_path) {
           if (!exact) {
               for(INT x = 0; x < nv; x++) {
               
                   partitionp[x] = partition[x];
               
               }
           }else {
               
               for(INT x = 0; x < nv; x++) {
               
                   ex_partition[x] = partition[x];
               
               }          
           
           }
       }
       exact = false;
    
    }
    
    INT crit_pathp = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    INT excrit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, ex_partition);
    
    if (crit_pathp < excrit_path) {
       for(INT i = 0; i < nv; i++) {
    
          partition[i] = partitionp[i];
    
       }
    }else{
       for(INT i = 0; i < nv; i++) {
    
          partition[i] = ex_partition[i];
    
       }
    }
    
    /* free section */
    
    
    free(is_in_halo);
    free(is_locked);
    free(partitionp);
    free(ex_partition);
    for(INT p = 0; p < k; p++){
        free(partition_moves[p]);
        free(partition_gains[p]);
    }
    free(partition_gain_sizes);
    free(to_balance);
    free(is_computed);
    free(cost);
    free(part_sizes);
    free(capacity_cst);
    
    
    return 0;   

}

INT dkfmFast(
Hypergraph * h,  
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance,
INT k) {

    bool verbose = false;
    bool unsize  = true;
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    /** Compute the halo of the partition **/
    
    bool is_in_halo; 
    INT * halo    = (INT*) malloc(sizeof(INT)*nv);
    MEM_ERROR(halo);
    INT halo_size = 0;
    
    bool * is_locked   = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_locked);
    
    
    INT * part_sizes   = (INT*)calloc(k*nwv, sizeof(INT));
    MEM_ERROR(part_sizes);
    
    INT * capacity_cst  = (INT*)calloc(nwv, sizeof(INT));
    MEM_ERROR(capacity_cst);
    
    float * cost = (float*) malloc(sizeof(float)*nv);
    MEM_ERROR(cost);
    
    for (INT u = 0; u < nv; u++) {        
        
        /* compute the part_sizes */
        for(INT wi = 0; wi < nwv; wi++) {
           part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
           capacity_cst[wi] += h->ti_weights[u*nwv+wi];
        }
    
        
        is_in_halo = false;
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
                    INT v = cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo = true;
            }
            
            cell = cell->next;
        }
        List * in_cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
        
            INT v = in_cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo = true;
            }
            
            in_cell = in_cell->next;
        }
        if (is_in_halo) {
            halo[halo_size++] = u;
        }
        is_locked[u] = false;
    }
    for(INT wi = 0; wi < nwv; wi++) {
        
        capacity_cst[wi] = ceil((float)capacity_cst[wi]/(float)k *1.05);
    }
    
    
    /** compute the moves arrays for each part **/
    
    //compute minimum and maximum gain
    INT * partitionp = (INT*)calloc(nv, sizeof(INT));
    MEM_ERROR(partitionp);
        
    INT min_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    
    for (INT u = 0; u < nv; u++) {
        partitionp[u] = u;
    }
    
    INT max_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    
    INT cur_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
    
    INT min_gain = cur_crit_path - max_crit_path; // \in [-max_crit_path, 0]
    INT max_gain = cur_crit_path - min_crit_path; // \in [0, max_crit_path]
    
    for (INT u = 0; u < nv; u++) {
        partitionp[u] = partition[u];
        
        INT D = 0;
        
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
        
            INT v = cell->i;
            
            if (a->ti_delay[partition[u]*a->i_m+partition[v]] > D) {
            
                D = a->ti_delay[partition[u]*a->i_m+partition[v]];
            }
            
            cell = cell->next;
            
        }
        
    
        /* compute cost */
        cost[u] = ((float)tolerance/100.0) * ( ((float)max_crit_path - (float)h->ti_criticalities_right[u]) / (float)D);
    }
    
    /* relative to 0 */
    
    INT gain_size = max_gain - min_gain;
    INT relative  = gain_size - min_gain;
    
    /* gain data structure declaration */
        
    INT ** partition_moves  = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_moves);
    for (INT p = 0; p < k; p++) {
        partition_moves[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_moves[p]);
    }
    
    INT ** partition_gains = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_gains);
    for (INT p = 0; p < k; p++) {
        partition_gains[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_gains[p]);
        for (INT i = 0; i < nv; i++) {
        
            partition_gains[p][i] = min_gain;
        }
        
    }    
    INT * partition_gain_sizes = (INT*)calloc(k, sizeof(INT));
    MEM_ERROR(partition_gain_sizes);
    
    bool * is_computed = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(is_computed);
    
    INT ncomp;
    INT anc_part;
    INT gain;
    List * cell, * in_cell, * cellp;
    
    INT u, v, vp;
    
    /* Compute the number of vertices evaluation in the halo */
    INT N_EVALS = ceil((float)halo_size*(1.0/(float)(perform-2)));
    
    /* Select randomly, vertices in the halo */
    for(INT x = 0; x < N_EVALS; x++) {
        
        int idx = rand() % halo_size;
        INT i = halo[idx];
        
               
        
            for(INT p = 0; p < k; p++) {
                 is_computed[p] = false;
            }
            
            ncomp = 0;
            
                	    
    	    cell   = neighbors[i];
    	    for(INT x = 0; x < neighbors[i]->size && ncomp<k; x++) {
            
                v = cell->i;
                if (partition[i] != partition[v]) {
                    if(!is_computed[partition[v]]) {
                        /* add u to list of gain of partition[v] */
                        anc_part = partition[i];
                        partition[i] = partition[v];
                        
                        gain = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                        
                        // add the vertices corresponding to its gain 
                        
                        partition_gains[partition[i]][i] = gain;
                        
                                                
                        pqueue_add_element(partition_moves[partition[i]], partition_gains[partition[i]], i, partition_gain_sizes[partition[i]]);
                        partition_gain_sizes[partition[i]]++;
                                       
                        partition[i] = anc_part;
                        
                        is_computed[partition[v]] = true;
                        ncomp++;
                    }
                }
                cell = cell->next;
            }
            in_cell   = in_neighbors[i];
            for(INT x = 0; x < in_neighbors[i]->size && ncomp<k; x++) {
            
                v = in_cell->i;
                if (partition[i] != partition[v]) {
                    if(!is_computed[partition[v]]) {
                        /* add u to list of gain of partition[v] */
                        anc_part = partition[i];
                        partition[i] = partition[v];
                        
                        gain = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                       
                        // add the vertices corresponding to its gain 
                        partition_gains[partition[v]][i] = gain;
                        
                        pqueue_add_element(partition_moves[partition[v]], partition_gains[partition[v]], i, partition_gain_sizes[partition[v]]);
                        partition_gain_sizes[partition[v]]++;
                                       
                        partition[i] = anc_part;
                        
                        is_computed[partition[v]] = true;
                        ncomp++;
                     }
                }
                in_cell = in_cell->next;
            }
        
          	
    
    }
    
    bool * to_balance = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(to_balance);
    
    for(INT p = 0; p < k; p++) {
        
        to_balance[p] = false;
        
        for(INT wi = 0; wi < nwv; wi++){
            if(part_sizes[p] >= capacity_cst[wi]){
                
                to_balance[p] = true;
            
            }   
        }
    }
    
    /* DKFM moves */
    
    INT base_cur_crit_path = cur_crit_path;
    INT crit_path_up = 0;
    INT best_crit_path = cur_crit_path;
    
    bool exact = false;
    
    INT * ex_partition = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(ex_partition);
    
    for(INT i = 0; i < nv; i++) {
       ex_partition[i] = partition[i];
    }
    
    INT N_MOVES = ceil((float)halo_size*(1.0/(float)perform));
    for(INT i = 0; i < N_MOVES; i++) {
        
               
    
    	/* select part */
    	INT selected_part = -1;
    	
    	for(INT p = 0; p < k; p++) {
            
            if(!to_balance[p] && partition_gain_sizes[p]>0) {
                
                selected_part = p;
                
            }
            
         }
         
         if (selected_part == -1) {
             break;
             break;
         }
         
         
            
         INT u = pqueue_dequeue(partition_moves[selected_part], partition_gains[selected_part], partition_gain_sizes[selected_part]);
         partition_gain_sizes[selected_part]--;
         is_locked[u] = true;
         cur_crit_path = base_cur_crit_path - partition_gains[selected_part][u];
            
         for(INT wi = 0; wi < nwv; wi++){
             part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
             part_sizes[partition[u]*nwv+wi]  -= h->ti_weights[u*nwv+wi];
             if(part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]){
                
                to_balance[selected_part] = true;
            
             }   
            
          }
            
            
          if (to_balance[partition[u]]) {
              to_balance[partition[u]] = false;
          }
          
          for(INT wi = 0; wi < nwv; wi++){

             if(part_sizes[partition[u]*nwv+wi] <= capacity_cst[wi]){
                
                to_balance[partition[u]] = false;
            
             }else {
                to_balance[partition[u]] = true;
             }   
            
          }
            
          /* update the current gain */ 
          
          cell = neighbors[u];
          for(INT x = 0; x < neighbors[u]->size && cost[u]<=0; x++) {
                          
             v = cell->i;
                
             if (!is_locked[v]) {
                bool is_connected = false;
                cellp = neighbors[v];
                for(INT o = 0; o < neighbors[v]->size; o++) {
                                
                    vp = cellp->i;
                    
                    if (partition[vp] == partition[u]) {
                           
                        // the future anciant part of u
                        is_connected = true;
                       
                    }
                    cellp = cellp->next;
                }
                if (!is_connected) {
                   /* delete v of anciant part of u moves */
                   partition_gains[partition[u]][v] = min_gain;
                   pqueue_siftUp(partition_moves[partition[u]], partition_gains[partition[u]], v, partition_gain_sizes[partition[u]]);

                }
                if (partition[v] != selected_part) {
                       
                    anc_part = partition[v];
                       
                    partition[u] = selected_part;
                    partition[v] = selected_part;
                       
                    partition_gains[selected_part][v] = cur_crit_path - compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
                                        
                    pqueue_siftUp(partition_moves[selected_part], partition_gains[selected_part], v, partition_gain_sizes[selected_part]);
                    
                     
                                
                    partition[v] = anc_part;
                       
                 }
              }
              
              cell = cell->next;
        
       } 
       partition[u] = selected_part;
       crit_path_up++;
       
       if (cost[u] <= 0 || crit_path_up >= 100) {
          cur_crit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition);
          exact = true;
          crit_path_up = 0;
       }
       
       if (cur_crit_path < best_crit_path) {
           if (!exact) {
               for(INT x = 0; x < nv; x++) {
               
                   partitionp[x] = partition[x];
               
               }
           }else {
               
               for(INT x = 0; x < nv; x++) {
               
                   ex_partition[x] = partition[x];
               
               }          
           
           }
       }
       exact = false;
    
    }
    
    INT crit_pathp = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partitionp);
    INT excrit_path = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, ex_partition);
    
    if (crit_pathp < excrit_path) {
       for(INT i = 0; i < nv; i++) {
    
          partition[i] = partitionp[i];
    
       }
    }else{
       for(INT i = 0; i < nv; i++) {
    
          partition[i] = ex_partition[i];
    
       }
    }
    
    /* free section */
    
    
    free(halo);
    free(is_locked);
    free(partitionp);
    free(ex_partition);
    for(INT p = 0; p < k; p++){
        free(partition_moves[p]);
        free(partition_gains[p]);
    }
    free(partition_gain_sizes);
    free(to_balance);
    free(is_computed);
    free(cost);
    free(part_sizes);
    free(capacity_cst);
    
    
    return 0;   
    


}

/* compute the critical path according to a partition and a target topology (partition) */
INT compute_partition_criticality(
Hypergraph * h, 
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition) {

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
    delays[i] = h->ti_delays[i];
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

    if (is_red[v]) {
      List * cell_neighbors = neighbors[v];
      for(INT x = 0; x < neighbors[v]->size; x++) {
      
          //selection of vertex v
          INT u = cell_neighbors->i;
          
          if(partition[u]!=partition[v]){
              ncuts++;
              delays[u] = MAX(delays[u], h->ti_delays[v]+h->ti_delays[u]+a->ti_delay[partition[u]*a->i_m+partition[v]]);
          } else {
              delays[u] = MAX(delays[u], h->ti_delays[v]+h->ti_delays[u]);
          }
          cell_neighbors = cell_neighbors->next;
     }
   }
   else {
      List * cell_neighbors = neighbors[v];
      for(INT x = 0; x < neighbors[v]->size; x++) {
                  //selection of vertex v
                  INT u = cell_neighbors->i;
                  
                  if(partition[u]!=partition[v]){
                      ncuts++;
                      delays[u] = MAX(delays[u], delays[v]+h->ti_delays[u]+a->ti_delay[partition[u]*a->i_m+partition[v]]);
                  } else {
                      delays[u] = MAX(delays[u], delays[v]+h->ti_delays[u]);
                  }
                  cell_neighbors = cell_neighbors->next;
              }

    }
  }

  INT max_v = -1;
  INT max_delays = 0;
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

/*compute the connectivity cut cost*/
INT compute_partition_cut(
Hypergraph * h, 
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition,
INT * lambda,
INT k){

  INT nv     = h->i_vertices;
  INT ne     = h->i_hyperedges;
  INT nr     = h->i_reds;
  INT np     = h->i_pins;
  INT nwv    = h->i_weights;
  INT * reds = h->ti_reds;
  
  INT * hyperedges     = h->ti_hyperedges;
  INT * idx_hyperedges = h->ti_idx_hyperedges; 

  INT cut_cost = 0;

  bool * is_in_part = (bool*)malloc(sizeof(bool)*k);
  MEM_ERROR(is_in_part);
  for (INT p = 0; p < k; p++) {
    is_in_part[p] = false;
  }

  INT lambda_cost;
  for (INT j = 0; j < ne; j++) {
    
    INT idx_j  = idx_hyperedges[j];
    INT wj     = hyperedges[idx_j];
    INT size_j = hyperedges[idx_j+1];
    
    lambda[j]  = 0;
    
    for (INT i = 0; i < size_j; i++) {
      INT e =hyperedges[idx_j+2+i];
      INT p = partition[e];
      is_in_part[p] = true;
    }

    lambda_cost = 0;
    for (INT p = 0; p < k; p++) {
      if(is_in_part[p]) {
        
        lambda_cost += 1;
        is_in_part[p] = false;
      }

    }

    lambda[j] = (lambda_cost-1)*MAX(1,wj);
    cut_cost += (lambda_cost-1)*MAX(1,wj);

  }
  
  free(is_in_part);

  return cut_cost;

}

INT kfm(
Hypergraph * h,  
Arch * a, 
List ** neighbors, 
List ** in_neighbors, 
INT * sort, 
INT * partition, 
INT perform, 
INT tolerance,
INT k) {

    bool verbose = false;
    bool unsize  = true;
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    INT * idx_hyperedges = h->ti_idx_hyperedges;
    INT * hyperedges     = h->ti_hyperedges;
    
    List ** vertex_to_hyper = (List**)malloc(sizeof(List*)*nv);
    MEM_ERROR(vertex_to_hyper);
    for(INT i = 0; i < nv; i++) {
        
        vertex_to_hyper[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(vertex_to_hyper[i]);
        new_list(vertex_to_hyper[i]);
    
    }
    
    
    bool is_in;
    for(INT j = 0; j < ne; j++) {
        
        INT idx_j  = idx_hyperedges[j];
        INT wj     = hyperedges[idx_j];
        INT size_j = hyperedges[idx_j+1];
        
        for (INT i = 0; i < size_j; i++) {
            INT u = hyperedges[idx_j+2+i];
            is_in = false;
            List * cell = vertex_to_hyper[u];
            while (cell != NULL) {
                if (cell->i == j) {
                    is_in = true;
                }
                cell = cell->next;
            }
            if(!is_in){
                list_add_element(vertex_to_hyper[u], j);
            }
            
        }
    }
    
    INT * neighbours_parts = (INT*)calloc(nv*k, sizeof(INT));
    MEM_ERROR(neighbours_parts);
    
    /** Compute the halo of the partition **/
    
    bool  * is_in_halo = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_in_halo);
    bool * is_locked   = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_locked);
    INT halo_size      = 0;
    
    INT * part_sizes   = (INT*)calloc(k*nwv, sizeof(INT));
    MEM_ERROR(part_sizes);
    
    INT * capacity_cst  = (INT*)calloc(nwv, sizeof(INT));
    MEM_ERROR(capacity_cst);
    
    INT * lambda = (INT*) calloc(ne, sizeof(INT));
    MEM_ERROR(lambda);
    
    for (INT u = 0; u < nv; u++) {
            
        /* compute the part_sizes */
        for(INT wi = 0; wi < nwv; wi++) {
           part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
           capacity_cst[wi] += h->ti_weights[u*nwv+wi];
        }
    
        
        is_in_halo[u] = false;
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
        
            INT v = cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo[u] = true;
            }
            
            neighbours_parts[u*k+partition[v]]++;
            cell = cell->next;            
        }
        List * in_cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            INT v = in_cell->i;
            if (partition[u] != partition[v]) {
                is_in_halo[u] = true;
            }
            
            neighbours_parts[u*k+partition[v]]++;
            in_cell = in_cell->next;
        }
        if (is_in_halo[u]) {
            halo_size++;
        }
        is_locked[u] = false;
    }
    for(INT wi = 0; wi < nwv; wi++) {
        
         capacity_cst[wi] = ceil((float)capacity_cst[wi]/(float)k *1.05);
    }
    
    /** compute the moves arrays for each part **/
    
    INT current_cut = compute_partition_cut(h, a, neighbors, in_neighbors, sort, partition, lambda, k);
                    
    /* gain data structure declaration */
        
    INT ** partition_moves  = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_moves);
    for (INT p = 0; p < k; p++) {
        partition_moves[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_moves[p]);
    }
    
    INT ** partition_gains = (INT**)malloc(sizeof(INT*)*k);
    MEM_ERROR(partition_gains);
    for (INT p = 0; p < k; p++) {
        partition_gains[p] = (INT*)calloc(nv, sizeof(INT));
        MEM_ERROR(partition_gains[p]);
        for (INT i = 0; i < nv; i++) {
        
            partition_gains[p][i] = 0 - ne;
        }
        
    }  
      
    INT * partition_gain_sizes = (INT*)calloc(k, sizeof(INT));
    MEM_ERROR(partition_gain_sizes);
    
    bool * is_in_part = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(is_in_part);
    
    for(INT p = 0; p < k; p++) {
       is_in_part[p] = false;
    }
    
    for(INT i = 0; i < nv; i++) {
        
        if (is_in_halo[i]) {
                    
            for(INT p = 0; p < k; p++) {
            
                if(p != partition[i]) {
            
                    INT current_cost = current_cut;
                    
                    INT anc_part = partition[i];
                    partition[i] = p;
                    
                    if (neighbours_parts[i*k+p]>0) {
                    
                        /* add u to list of gain of p */
                    
                     
                        List * cell = vertex_to_hyper[i];
                        while (cell != NULL && cell->size>0) {
                        
                            INT j      = cell->i;
                            INT idx_j  = idx_hyperedges[j];
                            INT wj     = hyperedges[idx_j];
                            INT size_j = hyperedges[idx_j+1];
        
                            for (INT x = 0; x < size_j; x++) {
                               is_in_part[partition[hyperedges[idx_j+2+x]]] = true;
                            }
                        
                            INT lambda_cost = 0;
                            for (INT pp = 0; pp < k; pp++) {
                                if(is_in_part[pp]) {
                        
                                    lambda_cost += 1;
                                    is_in_part[pp] = false;
                                 }
                        
                            }
                        
                            current_cost -= lambda[j];
                            current_cost += lambda_cost*wj;
                        
                            cell = cell->next;     

                        }                    
                    }
                    
                    partition[i] = anc_part;
                
                    INT gain = current_cut - current_cost;
                
                    // add the vertices corresponding to its gain 
                    partition_gains[p][i] = gain;
                    pqueue_add_element(partition_moves[p], partition_gains[p], i, partition_gain_sizes[p]);
                    partition_gain_sizes[p]++;
                }
            }	
        }
    }
    
    INT N_MOVES_MAX = 0;
    
    for(INT p = 0; p < k; p++) {
        N_MOVES_MAX += partition_gain_sizes[p];
    }
    
    bool * to_balance = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(to_balance);
    
    bool * locked_part = (bool*)malloc(sizeof(bool)*k);
    MEM_ERROR(locked_part);
    
    for(INT p = 0; p < k; p++) {
    
        to_balance[p] = false;
        locked_part[p] = false;
        
        for(INT wi = 0; wi < nwv; wi++){
            if(part_sizes[p] >= capacity_cst[wi]){
                locked_part[p] = true;
                to_balance[p] = true;
            
            }   
        }
    }
    
    INT * partitionp = (INT*)malloc(nv* sizeof(INT));
    MEM_ERROR(partitionp);
    
    memcpy(partitionp, partition, sizeof(INT)*nv);
    
    
    
    /* KFM moves */
    
    INT base_cur_cut_cost = current_cut;
    INT best_cut_cost = current_cut;
    
    
    INT N_MOVES = ceil((float)nv*(1.0/(float)perform));
    
    
    for(INT i = 0; i < N_MOVES && i < N_MOVES_MAX; i++) {
            
    	/* select part */
    	INT selected_part = -1;
    	
    	for(INT p = 0; p < k; p++) {
            
            if(!to_balance[p] && partition_gain_sizes[p]>0) {
                
                selected_part = p;
                
            }
            
         }
         if (selected_part == -1) {
             break;
             break;
         }
            
         INT u = pqueue_dequeue(partition_moves[selected_part], partition_gains[selected_part], partition_gain_sizes[selected_part]);
         partition_gain_sizes[selected_part]--;
         
         if (partition[u] != selected_part) {
         
            
            is_locked[u] = true;
            current_cut = base_cur_cut_cost - partition_gains[selected_part][u];
         
            /* update the lambda[j] values */
            INT anc_part_u = partition[u];
            partition[u]   = selected_part;
            List * cell_e  = vertex_to_hyper[i];
            while (cell_e != NULL && cell_e->size>0) {
                        
                 INT j      = cell_e->i;
                 assert(j<ne&&j>=0);
                 INT idx_j  = idx_hyperedges[j];
                 INT wj     = hyperedges[idx_j];
                 INT size_j = hyperedges[idx_j+1];
        
                 for (INT x = 0; x < size_j; x++) {
                    is_in_part[partition[hyperedges[idx_j+2+x]]] = true;
                 }
                        
                 INT lambda_cost = 0;
                 for (INT pp = 0; pp < k; pp++) {
                    if(is_in_part[pp]) {
                  
                       lambda_cost += 1;
                        is_in_part[pp] = false;
                    }
                 }
              
                 lambda[j] = lambda_cost*wj;
                        
                 cell_e = cell_e->next;     

            }
         
            partition[u] = anc_part_u;
         
            neighbours_parts[u*k+partition[u]]--;
            neighbours_parts[u*k+selected_part]++;
                     
            for(INT wi = 0; wi < nwv; wi++){
         
             
         
                part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
                part_sizes[partition[u]*nwv+wi]  -= h->ti_weights[u*nwv+wi];
             
                
                if(part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]){
                
                   to_balance[selected_part] = true;
            
                }   
             
             
             }
            
            
             if (to_balance[partition[u]]) {
                 to_balance[partition[u]] = false;
             }
            
             /* update the current gain */ 
          
             partition[u] = selected_part;
          
             List * cell = neighbors[u];
             for(INT o = 0; o < neighbors[u]->size; o++) {
                          
                INT v = cell->i;
                                
                if (!is_locked[v]) {
             
                   for(INT p = 0; p < k; p++) {
            
                      INT current_cost = current_cut;
                
                      INT anc_part = partition[v];
                      partition[v] = p;
                
                      if (neighbours_parts[v*k+p]>0) {
                    
                         /* add u to list of gain of p */
                        
                         List * cell_e = vertex_to_hyper[v];
                         while (cell_e != NULL && cell_e->size>0) {
                        
                            INT j      = cell_e->i;
                            assert(j<ne&&j>=0);
                            INT idx_j  = idx_hyperedges[j];
                            INT wj     = hyperedges[idx_j];
                            INT size_j = hyperedges[idx_j+1];
        
                            for (INT x = 0; x < size_j; x++) {
                                is_in_part[partition[hyperedges[idx_j+2+x]]] = true;
                            }
                        
                            INT lambda_cost = 0;
                            for (INT pp = 0; pp < k; pp++) {
                                if(is_in_part[pp]) {
                        
                                    lambda_cost += 1;
                                    is_in_part[pp] = false;
                                }
                            }
                        
                            current_cost -= lambda[j];
                            current_cost += lambda_cost*wj;
                        
                            cell_e = cell_e->next;     

                        }                    
                    }
                
                    INT gain = current_cut - current_cost;
                    pqueue_siftUp(partition_moves[p], partition_gains[p], partition_gain_sizes[p]-1, partition_gain_sizes[p]);
                     
                }
             }
             cell = cell->next;
          }
            
          partition[u] = selected_part;
       
       
       
          if (current_cut < best_cut_cost) {
             best_cut_cost = current_cut;
             
             printf("new best cost %d\n", best_cut_cost);
             for(INT x = 0; x < nv; x++) {
                
                partitionp[x] = partition[x];
               
             }
          }
       }
    
    }
    
    for(INT i = 0; i < nv; i++) {
    
        partition[i] = partitionp[i];
    
    }
    
   
    
    /* free section */
    
    
    free(is_in_halo);
    free(is_locked);
    free(partitionp);
    free(neighbours_parts);
    for(INT p = 0; p < k; p++){
        free(partition_moves[p]);
        free(partition_gains[p]);
    }
    free(partition_gain_sizes);
    free(to_balance);
    free(is_in_part);
    free(lambda);
    free(locked_part);
    free(part_sizes);
    free(capacity_cst);
    
    for(INT i = 0; i < nv; i++) {
        delete_list(vertex_to_hyper[i]);
    }
    free(vertex_to_hyper);
    
    
    return 0;   
    


}




