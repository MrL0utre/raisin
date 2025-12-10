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
/**   NAME       : ipart.c                                 **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are initial partitioning    **/
/**                functions definitions for red-black     **/
/**                hypergraph.                             **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/
#include "ipart.h"
int derivedBreadthFirstSearch(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * partition, 
    INT k, 
    INT epsilon) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);

    INT * main_pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(main_pqueue);

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
        pqueue[i]      =-1;
        main_pqueue[i] =-1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
    }

    INT main_pqueue_start = 0;
    INT main_pqueue_end   = 0;
    for(u = 0; u < nv; u++) {
        if (indeg[u] == 0) {

            pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
            main_pqueue_end++;
            flag[u]                      = true;
        }
    }
    INT index = 0, pqueue_start = 0, pqueue_end = 0;
   
    

    if (main_pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
        main_pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)calloc(k*nwv, sizeof(INT));
    MEM_ERROR(partition_size);
    
    INT * limit_size     = (INT*)calloc(nwv, sizeof(INT));
    MEM_ERROR(limit_size);
    
    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }

    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ;
    }

    index = 0;
    /* start exploration */
    while (main_pqueue_start != main_pqueue_end) {
        pqueue_start = 0;
        pqueue_end   = 0;
        for (i = 0; i < main_pqueue_end; i++) {
            // all red vertices visited are explored next
            pqueue[pqueue_end++] = main_pqueue[i];
        }
        main_pqueue_start = 0;
        main_pqueue_end   = 0;

        while(pqueue_end > 0) {
            u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
            pqueue_end--;
            bool new_part = false;
            for(INT wi=0; wi < nwv; wi++) {
                if (partition_size[index*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                    new_part = true;
                }
            }
            if (new_part){
                index++;
            }
            partition[u]  = index;
            for(INT wi=0; wi < nwv; wi++) {
                partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
            }
            flag[u]       = true;
            List * cell   = neighbors[u];
            for(INT x = 0; x < neighbors[u]->size; x++) {
            
                v = cell->i;
                indeg[v]--;
                if (indeg[v] <= 0 && !is_red[v] && !flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (is_red[v] && !flag[v]) {
                    /* insert by criticality */
                    pqueue_add_element(main_pqueue, criticalities, v, main_pqueue_end);
                    main_pqueue_end++;
                    flag[v]            = true;
                }
                cell = cell->next;
            }

        }
    }
    
    for(INT u = 0; u < nv; u++){
        
        assert(flag[u]);
    }
    

    /* free section */
    free(pqueue);
    free(main_pqueue);
    free(indeg);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);

    return 0;

}

int derivedBreadthFirstSearchMultilevel(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * partition, 
    INT k, 
    INT epsilon) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);

    INT * main_pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(main_pqueue);

    INT * indeg = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(indeg);
    bool * flag      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    bool * is_red      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);
    INT * connected_parts = (INT*)calloc(k, sizeof(INT));
    MEM_ERROR(connected_parts);

    for(i = 0; i < nv; i++) {
        indeg[i]      = in_neighbors[i]->size;
        flag[i]       = false;
        is_red[i]     = false;
        pqueue[i]      =-1;
        main_pqueue[i] =-1;
        partition[i]  = -1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
    }

    INT main_pqueue_start = 0;
    INT main_pqueue_end   = 0;
    for(u = 0; u < nv; u++) {
        if (indeg[u] == 0) {

            pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
            main_pqueue_end++;
            flag[u]                      = true;
        }
    }
    INT index = 0, pqueue_start = 0, pqueue_end = 0;
   
    

    if (main_pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
        main_pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)calloc(k*nwv, sizeof(INT));
    MEM_ERROR(partition_size);
    
    INT * limit_size     = (INT*)calloc(nwv, sizeof(INT));
    MEM_ERROR(limit_size);
    
    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }

    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ;
    }

    index = 0;
    /* start exploration */
    while (main_pqueue_start != main_pqueue_end) {
        pqueue_start = 0;
        pqueue_end   = 0;
        for (i = 0; i < main_pqueue_end; i++) {
            // all red vertices visited are explored next
            pqueue[pqueue_end++] = main_pqueue[i];
        }
        main_pqueue_start = 0;
        main_pqueue_end   = 0;
        
        INT smallest_part = 0;
        

        while(pqueue_end > 0) {
            u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
            pqueue_end--;
            index = -1;
            for(INT p = 0;  p < k; p++) {
                bool size_is_ok = true;
                for(INT wi=0; wi < nwv; wi++) {
                    if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                        size_is_ok = false;
                    }
                    if (partition_size[p*nwv+wi] < partition_size[smallest_part*nwv+wi]){
                        smallest_part = p;
                    }
                }
                if (size_is_ok) {
                    index = p; 
                }
            }
            if(index==-1) {
                index = smallest_part;
            }
            partition[u]  = index;
            for(INT wi=0; wi < nwv; wi++) {
                partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
            }
            flag[u]       = true;
            List * cell   = neighbors[u];
            for(INT x = 0; x < neighbors[u]->size; x++) {
            
                v = cell->i;
                indeg[v]--;
                if (indeg[v] <= 0 && !is_red[v] && !flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (is_red[v] && !flag[v]) {
                    /* insert by criticality */
                    pqueue_add_element(main_pqueue, criticalities, v, main_pqueue_end);
                    main_pqueue_end++;
                    flag[v]            = true;
                }
                cell = cell->next;
            }

        }
    }
    
    main_pqueue_start = 0;
    main_pqueue_end   = 0;
    for(u = 0; u < nv; u++) {
        if (!flag[u] && is_red[u]) {

            pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
            main_pqueue_end++;
            flag[u]                      = true;
        }
    }
    pqueue_start = 0, pqueue_end = 0;
    
    
    
    /* restart exploration */
    while (main_pqueue_start != main_pqueue_end) {
        pqueue_start = 0;
        pqueue_end   = 0;
        for (i = 0; i < main_pqueue_end; i++) {
            // all red vertices visited are explored next
            pqueue[pqueue_end++] = main_pqueue[i];
        }
        main_pqueue_start = 0;
        main_pqueue_end   = 0;
        
        INT smallest_part = 0;
        

        while(pqueue_end > 0) {
            u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
            pqueue_end--;
            index = -1;
            for(INT p = 0;  p < k; p++) {
                bool size_is_ok = true;
                for(INT wi=0; wi < nwv; wi++) {
                    if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                        size_is_ok = false;
                    }
                    if (partition_size[p*nwv+wi] < partition_size[smallest_part*nwv+wi]){
                        smallest_part = p;
                    }
                }
                if (size_is_ok) {
                    index = p; 
                }
            }
            if(index==-1) {
                index = smallest_part;
            }
            
            flag[u]       = true;
            List * cell   = neighbors[u];
            for(INT x = 0; x < neighbors[u]->size; x++) {
            
                v = cell->i;
                indeg[v]--;
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (flag[v] && partition[v] > 0) {
                    connected_parts[partition[v]] ++;
                }
                
                cell = cell->next;
            }
            
            cell   = in_neighbors[u];
            for(INT x = 0; x < in_neighbors[u]->size; x++) {
            
                v = cell->i;
                
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (flag[v] && partition[v] > 0) {
                    connected_parts[partition[v]] ++;
                }
                
                cell = cell->next;
            }
            INT max_p = -1;
            
            for(INT p = 0;  p < k; p++) {
                bool size_is_ok = true;
                for(INT wi=0; wi < nwv; wi++) {
                    if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                        size_is_ok = false;
                    }
                }
            
            
                if (connected_parts[p] > max_p && size_is_ok) {
                    index = p;
                    max_p = connected_parts[p];
                }
                connected_parts[p] = 0;
            }
           
            partition[u]  = index;
            for(INT wi=0; wi < nwv; wi++) {
                partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
            }

        }
    }
    
    
    main_pqueue_start = 0;
    main_pqueue_end   = 0;
    /*  */
    for(u = 0; u < nv; u++) {
        if (!flag[u]) {

            pqueue_add_element(main_pqueue, criticalities, u, main_pqueue_end);
            main_pqueue_end++;
            flag[u]                      = true;
        }
    }
    pqueue_start = 0, pqueue_end = 0;
 
    /* restart exploration */
    while (main_pqueue_start != main_pqueue_end) {
        pqueue_start = 0;
        pqueue_end   = 0;
        for (i = 0; i < main_pqueue_end; i++) {
            // all red vertices visited are explored next
            pqueue[pqueue_end++] = main_pqueue[i];
        }
        main_pqueue_start = 0;
        main_pqueue_end   = 0;
        
        INT smallest_part = 0;

        while(pqueue_end > 0) {
            u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
            pqueue_end--;
            index = -1;
            for(INT p = 0;  p < k; p++) {
                bool size_is_ok = true;
                for(INT wi=0; wi < nwv; wi++) {
                    if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                        size_is_ok = false;
                    }
                    if (partition_size[p*nwv+wi] < partition_size[smallest_part*nwv+wi]){
                        smallest_part = p;
                    }
                }
                if (size_is_ok) {
                    index = p; 
                }
            }
            if(index==-1) {
                index = smallest_part;
            }
            
            flag[u]       = true;
            List * cell   = neighbors[u];
            for(INT x = 0; x < neighbors[u]->size; x++) {
            
                v = cell->i;
                indeg[v]--;
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (flag[v] && partition[v] > 0) {
                    connected_parts[partition[v]] ++;
                }
                
                cell = cell->next;
            }
            
            cell   = in_neighbors[u];
            for(INT x = 0; x < in_neighbors[u]->size; x++) {
            
                v = cell->i;
                
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;

                }
                if (flag[v] && partition[v] > 0) {
                    connected_parts[partition[v]] ++;
                }
                
                cell = cell->next;
            }
            INT max_p = -1;
            
            for(INT p = 0;  p < k; p++) {
                bool size_is_ok = true;
                for(INT wi=0; wi < nwv; wi++) {
                    if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                        size_is_ok = false;
                    }
                }
            
            
                if (connected_parts[p] > max_p && size_is_ok) {
                    index = p;
                    max_p = connected_parts[p];
                }
                connected_parts[p] = 0;
            }
            partition[u]  = index;
            for(INT wi=0; wi < nwv; wi++) {
                partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
            }

        }
    }
    
    
    for(INT u = 0; u < nv; u++){
        assert(flag[u]);
    }
    

    /* free section */
    free(connected_parts);
    free(pqueue);
    free(main_pqueue);
    free(indeg);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);

    return 0;

}

int derivedDepthFirstSearch(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * partition, 
    INT k, 
    INT epsilon) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);


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
        pqueue[i]      =-1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
    }

    INT index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (indeg[u] == 0) {

            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }


    if (pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)malloc(sizeof(INT)*(k*nwv));
    MEM_ERROR(partition_size);
    for(INT i = 0; i < k; i++){
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[i*nwv+wi]=0;
        }
    }
    INT * limit_size     = (INT*)malloc(sizeof(INT)*nwv);
    MEM_ERROR(limit_size);
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi]=0;
    }

    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }

    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ;
    }

    index = 0;
    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        bool new_part = false;
        for(INT wi=0; wi < nwv; wi++) {
            if (partition_size[index*nwv+wi] >= limit_size[wi]){
                new_part = true;
            }
        }
        if (new_part){
            index++;
        }
        partition[u]  = index;
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
        }
        flag[u]       = true;
        List * cell   = neighbors[u];
        while(cell != NULL && neighbors[u]->size>0) {
            v = cell->i;
            indeg[v]--;
            if (indeg[v] <= 0 && !flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;

            }
            cell = cell->next;
        }
    }


    for(INT u = 0; u < nv; u++){

        assert(flag[u]);
    }
    /* free section */
    free(pqueue);
    free(indeg);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);

    return 0;



}

int derivedDepthFirstSearchMultilevel(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * partition, 
    INT k, 
    INT epsilon) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);


    INT * indeg = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(indeg);
    bool * flag      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    bool * is_red      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);
    
    INT * neighbors_parts = (INT*)calloc(k, sizeof(INT));
    MEM_ERROR(neighbors_parts);

    for(i = 0; i < nv; i++) {
        indeg[i]      = in_neighbors[i]->size;
        flag[i]       = false;
        is_red[i]     = false;
        pqueue[i]     = -1;
        partition[i]  = -1; 
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
    }

    INT index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (indeg[u] == 0) {

            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
        
    }


    if (pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)malloc(sizeof(INT)*(k*nwv));
    MEM_ERROR(partition_size);
    for(INT i = 0; i < k; i++){
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[i*nwv+wi]=0;
        }
    }
    INT * limit_size     = (INT*)malloc(sizeof(INT)*nwv);
    MEM_ERROR(limit_size);
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi]=0;
    }

    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }

    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ;
    }

    index = 0;
    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        bool new_part = false;
        for(INT wi=0; wi < nwv; wi++) {
            if (partition_size[index*nwv+wi] >= limit_size[wi]){
                new_part = true;
            }
        }
        if (new_part){
            index++;
        }
        partition[u]  = index;
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
        }
        flag[u]       = true;
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
        
            v = cell->i;
            indeg[v]--;
            if (indeg[v] <= 0 && !flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;

            }
            cell = cell->next;
        }
    }
    
    index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (is_red[u] && !flag[u]) {

            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }
    
    INT smallest_part = 0;

    /* restart exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        
        index = -1;
        for(INT p = 0;  p < k; p++) {
            bool size_is_ok = true;
            for(INT wi=0; wi < nwv; wi++) {
                if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                    size_is_ok = false;
                }
                if (partition_size[p*nwv+wi] < partition_size[smallest_part*nwv+wi]){
                    smallest_part = p;
                }
            }
            if (size_is_ok) {
                index = p; 
            }
        }
        if(index==-1) {
            index = smallest_part;
        }
        

        flag[u]       = true;
        List * cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
            v = cell->i;
            indeg[v]--;
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;

            }
            if (flag[v] && partition[v] > 0) {
                neighbors_parts[partition[v]] ++;
            }
            cell = cell->next;
        }
        cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v = cell->i;
            
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;

            }
            if (flag[v] && partition[v] > 0) {
                neighbors_parts[partition[v]] ++;
            }
            cell = cell->next;
            
        }
        INT max_p = -1;
            
        for(INT p = 0;  p < k; p++) {
            bool size_is_ok = true;
            for(INT wi=0; wi < nwv; wi++) {
                if (partition_size[p*nwv+wi]+h->ti_weights[u*nwv+wi] >= limit_size[wi]){
                    size_is_ok = false;
                }
            }
            if (neighbors_parts[p] > max_p && size_is_ok) {
                index = p;
                max_p = neighbors_parts[p];
            }
            neighbors_parts[p] = 0;
        }
        partition[u]  = index;
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[index*nwv+wi] += h->ti_weights[u*nwv+wi];
        }
        
        
    }

     for(INT u = 0; u < nv; u++){

        assert(flag[u]);
    }
    
    /* free section */
    free(neighbors_parts);
    free(pqueue);
    free(indeg);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);

    return 0;



}

int criticalConnectedComponentPartitioning(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors, 
    INT * umap, 
    INT k, 
    INT epsilon, 
    INT bound) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);


    INT * outdeg = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(outdeg);
    INT * partition = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(partition);
    bool * flag      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    bool * is_red      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);
    
    INT * r_candidates = (INT*)malloc(sizeof(INT)*nr);
    MEM_ERROR(r_candidates);
    INT size_r_candidates = 0;

    for(i = 0; i < nv; i++) {
        outdeg[i]     = neighbors[i]->size;
        flag[i]       = false;
        is_red[i]     = false;
        pqueue[i]     =-1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
        if (criticalities[reds[i]] > bound) {
        
            r_candidates[size_r_candidates++] = reds[i];
        
        }
    }


    INT index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (outdeg[u] == 0 && criticalities[u] > bound) {

            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }


    if (pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        while(criticalities[u] < bound && size_r_candidates > 0) {
            i_start       = rand() % size_r_candidates;
            u             = r_candidates[i_start];
        }
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)malloc(sizeof(INT)*(nv*nwv));
    MEM_ERROR(partition_size);
    for(INT i = 0; i < nv; i++){
        partition[i] = i;
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[i*nwv+wi] = h->ti_weights[i*nwv+wi];
        }
    }
    INT * limit_size     = (INT*)malloc(sizeof(INT)*nwv);
    MEM_ERROR(limit_size);
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi]=0;
    }

    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ; 
    }
    
    Union_find ** connected_components = (Union_find**) malloc(sizeof(Union_find*)*nv);
    MEM_ERROR(connected_components);
    for(INT i = 0 ; i < nv; i++) {
        connected_components[i] = new_set(i);
    }

    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];

        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            if (criticalities[v] > bound) {
                bool fusionable = true;
                INT partu = find(connected_components[u])->value;
                INT partv = find(connected_components[v])->value;

                if (partu != partv) {
                    // merge partition[u] with partition[v]
                    for(INT wi = 0; wi < nwv; wi++) {
                        if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                            fusionable = false;
                       }
                    }
                    if (fusionable) {
                     
                        Union_find * new_root = merge(connected_components[u], connected_components[v]);
                        if(new_root->value == partu) {
                            for(INT wi = 0; wi < nwv; wi++) {
                                partition_size[new_root->value*nwv+wi] += partition_size[partv*nwv+wi];
                            }
                        } else {
                            for(INT wi = 0; wi < nwv; wi++) {
                                partition_size[new_root->value*nwv+wi] += partition_size[partu*nwv+wi];
                            }
                        }
                        
                    }
                }
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;
                }
            }
            cell = cell->next;
        }
    }

    for(INT i = 0; i < nv; i++) {
        flag[i] = false;
    }

    index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (outdeg[u] == 0) {
            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }

    if (pqueue_end == 0) {
        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;
    }

    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                        }
                    }
                }
            }
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;
            }
            cell = cell->next;
        }
    }

    /* refinement */
    for(INT u = 0; u < nv; u++) {

        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {          
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                            
                        }
                   }
                   
               }
                   
            }
            cell = cell->next;
        }
        cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
            v               = cell->i;
            bool fusionable = true;
            
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                            
                        }
                   }
                   
               }
                   
            }
            cell = cell->next;
        }
    }

    for(INT u = 0; u < nv; u++){
        is_red[u] = false;
        assert(flag[u]);
        partition[u] = find(connected_components[u])->value;
        for(INT wi = 0; wi < nwv; wi++) {
            assert(partition_size[partition[u]*nwv+wi] <= limit_size[wi]);
        }
        
    }
    /* using is_red as bool array to re index clusters */
    for(INT u = 0; u < nv; u++){
        is_red[partition[u]] = true;
        
    }
    INT cpt = 0;
    INT max_p = 0;
    for(INT u = 0; u < nv; u++){
        if (is_red[partition[u]]) {
            INT p = partition[u];
            for(INT v = 0; v < nv; v++){
                if(partition[v] == p){
                    umap[v] = cpt;
                }
            }
            cpt++;
            is_red[p] = false;
        }
    }


    /* free section */

    for(INT i = 0 ; i < nv; i++) {
        free(connected_components[i]);
    }
    free(connected_components);
    free(pqueue);
    free(outdeg);
    free(partition);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);
    free(r_candidates);

    return 0;



}

int criticalConnectedComponentPartitioningMultilevel(
    Hypergraph * h, 
    List ** neighbors, 
    List ** in_neighbors,
    INT * umap, 
    INT k, 
    INT epsilon, 
    INT bound) {

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;

    INT * criticalities = h->ti_criticalities_right;

    INT i,u,v;

    /* priority queue */
    INT * pqueue = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(pqueue);


    INT * outdeg = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(outdeg);
    INT * partition = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(partition);
    bool * flag      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(flag);
    bool * is_red      = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_red);
    
    INT * r_candidates = (INT*)malloc(sizeof(INT)*nr);
    MEM_ERROR(r_candidates);
    INT size_r_candidates = 0;

    for(i = 0; i < nv; i++) {
        outdeg[i]     = neighbors[i]->size;
        flag[i]       = false;
        is_red[i]     = false;
        pqueue[i]     =-1;
    }
    for(i = 0; i < nr; i++) {
        is_red[reds[i]]     = true;
        if (criticalities[reds[i]] > bound) {
        
            r_candidates[size_r_candidates++] = reds[i];
        
        }
    }


    INT index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (outdeg[u] == 0 && criticalities[u] > bound) {

            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }


    if (pqueue_end == 0) {

        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        while(criticalities[u] < bound && size_r_candidates > 0) {
            i_start       = rand() % size_r_candidates;
            u             = r_candidates[i_start];
        }
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;

    }

    INT * partition_size = (INT*)malloc(sizeof(INT)*(nv*nwv));
    MEM_ERROR(partition_size);
    for(INT i = 0; i < nv; i++){
        partition[i] = i;
        for(INT wi=0; wi < nwv; wi++) {
            partition_size[i*nwv+wi] = h->ti_weights[i*nwv+wi];
        }
    }
    INT * limit_size     = (INT*)malloc(sizeof(INT)*nwv);
    MEM_ERROR(limit_size);
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi]=0;
    }

    for(INT i = 0; i < nv; i++) {

        for(INT wi=0; wi < nwv; wi++) {
            limit_size[wi] += h->ti_weights[i*nwv+wi];
        }

    }
    for(INT wi = 0; wi < nwv; wi++){
        limit_size[wi] = ceil(((float)limit_size[wi]/(float)k )*(1+(float)epsilon/100)) ; 
    }
    
    Union_find ** connected_components = (Union_find**) malloc(sizeof(Union_find*)*nv);
    MEM_ERROR(connected_components);
    for(INT i = 0 ; i < nv; i++) {
        connected_components[i] = new_set(i);
    }

    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];

        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            if (criticalities[v] > bound) {
                bool fusionable = true;
                INT partu = find(connected_components[u])->value;
                INT partv = find(connected_components[v])->value;

                if (partu != partv) {
                    // merge partition[u] with partition[v]
                    for(INT wi = 0; wi < nwv; wi++) {
                        if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                            fusionable = false;
                       }
                    }
                    if (fusionable) {
                     
                        Union_find * new_root = merge(connected_components[u], connected_components[v]);
                        if(new_root->value == partu) {
                            for(INT wi = 0; wi < nwv; wi++) {
                                partition_size[new_root->value*nwv+wi] += partition_size[partv*nwv+wi];
                            }
                        } else {
                            for(INT wi = 0; wi < nwv; wi++) {
                                partition_size[new_root->value*nwv+wi] += partition_size[partu*nwv+wi];
                            }
                        }
                        
                    }
                }
                if (!flag[v]) {
                    flag[v]            = true;
                    /* insert by criticality */
                    pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                    pqueue_end++;
                }
            }
            cell = cell->next;
        }
    }

    for(INT i = 0; i < nv; i++) {
        flag[i] = false;
    }

    index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (outdeg[u] == 0) {
            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }

    if (pqueue_end == 0) {
        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;
    }

    /* start exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                        }
                    }
                }
            }
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;
            }
            cell = cell->next;
        }
    }
    
    index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (is_red[u] && !flag[u]) {
            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }

    if (pqueue_end == 0) {
        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;
    }

    /* restart exploration */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                        }
                    }
                }
            }
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;
            }
            cell = cell->next;
        }
    }
    
    index = 0, pqueue_start = 0, pqueue_end = 0;
    for(u = 0; u < nv; u++) {
        if (!flag[u]) {
            pqueue_add_element(pqueue, criticalities, u, pqueue_end);
            pqueue_end++;
            flag[u]                      = true;
        }
    }

    if (pqueue_end == 0) {
        // thorus form red-black hypergraph
        // select a random red vertex
        INT i_start       = rand() % nr;
        u                 = reds[i_start];
        pqueue_add_element(pqueue, criticalities, u, pqueue_end);
        pqueue_end++;
        flag[u]           = true;
    }

    /*  */
    while(pqueue_end > 0) {
        u             = pqueue_dequeue(pqueue, criticalities, pqueue_end);
        pqueue_end--;
        flag[u]       = true;
        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                        }
                    }
                }
            }
            if (!flag[v]) {
                flag[v]            = true;
                /* insert by criticality */
                pqueue_add_element(pqueue, criticalities, v, pqueue_end);
                pqueue_end++;
            }
            cell = cell->next;
        }
    }

    /* refinement */
    for(INT u = 0; u < nv; u++) {

        List * cell   = in_neighbors[u];
        for(INT x = 0; x < in_neighbors[u]->size; x++) {
            v               = cell->i;
            bool fusionable = true;
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {          
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                            
                        }
                   }
                   
               }
                   
            }
            cell = cell->next;
        }
        cell   = neighbors[u];
        for(INT x = 0; x < neighbors[u]->size; x++) {
            v               = cell->i;
            bool fusionable = true;
            
            INT partu = find(connected_components[u])->value;
            INT partv = find(connected_components[v])->value;

            if (partu != partv) {
                // merge partition[u] with partition[v]
                for(INT wi = 0; wi < nwv; wi++) {
                    if (partition_size[partu*nwv+wi] + partition_size[partv*nwv+wi] > limit_size[wi]) {
                        fusionable = false;
                    }
                }
                if (fusionable) {
                    Union_find * new_root = merge(connected_components[u], connected_components[v]);
                    if(new_root->value == partu) {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partv*nwv+wi];
                        }
                    } else {
                        for(INT wi = 0; wi < nwv; wi++) {
                            partition_size[new_root->value *nwv+wi] += partition_size[partu*nwv+wi];
                            
                        }
                   }
                   
               }
                   
            }
            cell = cell->next;
        }
    }

    for(INT u = 0; u < nv; u++){
        is_red[u] = false;
        assert(flag[u]);
        partition[u] = find(connected_components[u])->value;
        for(INT wi = 0; wi < nwv; wi++) {
            assert(partition_size[partition[u]*nwv+wi] <= limit_size[wi]);
        }
        
    }
    /* using is_red as bool array to re index clusters */
    for(INT u = 0; u < nv; u++){
        is_red[partition[u]] = true;
        
    }
    INT cpt = 0;
    INT max_p = 0;
    for(INT u = 0; u < nv; u++){
        if (is_red[partition[u]]) {
            INT p = partition[u];
            for(INT v = 0; v < nv; v++){
                if(partition[v] == p){
                    umap[v] = cpt;
                }
            }
            cpt++;
            is_red[p] = false;
        }
    }


    /* free section */

    for(INT i = 0 ; i < nv; i++) {
        free(connected_components[i]);
    }
    free(connected_components);
    free(pqueue);
    free(outdeg);
    free(partition);
    free(flag);
    free(is_red);
    free(partition_size);
    free(limit_size);
    free(r_candidates);

    return 0;



}
