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


#include"../include/rbh.h"
#include"../include/a.h"
#include"../include/crbh.h"
#include"../include/pqueue.h"
#include"../include/list.h"
#include <stdlib.h>
#include <time.h>

int main(void){

    Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h);

    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);


    Arch * a = (Arch*)malloc(sizeof(Arch));
    MEM_ERROR(a);

    char * file_path = (char*)malloc(sizeof(char)*25);
    MEM_ERROR(file_path);
    strcpy(file_path, "../hypergraphs/b14.rzn2");

    char * arch_path = (char*)malloc(sizeof(char)*25);
    MEM_ERROR(arch_path);
    strcpy(arch_path, "../targets/arch0.arch");




    rbhLoad(h, file_path, 0, false);

    archLoad(a, arch_path, true);

    h->s_rbh_name = (char*)malloc(sizeof(char)*strlen(file_path));
    MEM_ERROR(h->s_rbh_name);

    a->s_arch_name = (char*)malloc(sizeof(char)*strlen(arch_path));
    MEM_ERROR(a->s_arch_name);


    strcpy(h->s_rbh_name,  "b14");
    strcpy(a->s_arch_name, "arch0");


    rbhSave(h, h->s_rbh_name, 0, true);

    INT epsilon = 5;
    INT k       = 4;
    INT * map   = (INT*)malloc(sizeof(INT)*h->i_vertices);
    MEM_ERROR(map);


    List ** neighbors_list = (List**)malloc(sizeof(List*)*h->i_vertices);
    MEM_ERROR(neighbors_list);
    
    List ** base_neighbors_list = (List**)malloc(sizeof(List*)*h->i_vertices);
    MEM_ERROR(base_neighbors_list);

    for (INT i = 0; i < h->i_vertices; i++) {
        neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(neighbors_list[i]);
        new_list(neighbors_list[i]);
        
        base_neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(base_neighbors_list[i]);
        new_list(base_neighbors_list[i]);
    }

    List ** in_neighbors_list = (List**)malloc(sizeof(List*)*h->i_vertices);
    MEM_ERROR(in_neighbors_list);
    
    List ** base_in_neighbors_list = (List**)malloc(sizeof(List*)*h->i_vertices);
    MEM_ERROR(base_in_neighbors_list);

    for (INT i = 0; i < h->i_vertices; i++) {
        
        in_neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(in_neighbors_list[i]);
        new_list(in_neighbors_list[i]);
    
    
        base_in_neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(base_in_neighbors_list[i]);
        new_list(base_in_neighbors_list[i]);
    }


    computeListNeighbors(h, neighbors_list);
    printf("neighbours\n");
    computeListInNeighbors(h, neighbors_list, in_neighbors_list);
    printf("in_neighbours\n");
    INT * sort = (INT*)malloc(sizeof(INT)*h->i_vertices);
    MEM_ERROR(sort);

    topologicalSort(h, neighbors_list, in_neighbors_list, sort);

    compute_criticality(h, neighbors_list, in_neighbors_list, sort);


    for (INT i = 0; i < h->i_vertices; i++) {
        map[i]=i;
    }

    int pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

     printf("\nmax pmax;%d\n#clusters;%d\n", pmax, h->i_vertices);

    compute_hypergraph_clustering(h, h2, map, h->i_vertices);

    assert(h->i_vertices == h2->i_vertices);

    assert(h->i_hyperedges == h2->i_hyperedges);
    assert(h->i_reds == h2->i_reds);
    assert(h->i_pins == h2->i_pins);
    assert(h->i_weights == h2->i_weights);


    for(INT i = 0; i < h->i_vertices; i++){
        assert(h->ti_criticalities_right[i] == h2->ti_criticalities_right[i]);
        assert(h->ti_criticalities_left[i] == h2->ti_criticalities_left[i]);
        if(h->ti_delays[i] != h2->ti_delays[i]) {
            printf("error : %d, delays= %d : %d\n", i, h->ti_delays[i], h2->ti_delays[i]);
        }
        assert(h->ti_delays[i] == h2->ti_delays[i]);
        for(INT wi = 0; wi < h->i_weights; wi++){
            assert(h->ti_weights[i*h->i_weights+wi] == h2->ti_weights[i*h->i_weights+wi]);
        }

    }
    for(INT i = 0; i < h->i_reds; i++){

        assert(h->ti_reds[i] == h2->ti_reds[i]);
    }
    for(INT j = 0; j < h->i_hyperedges; j++){
        assert(h->ti_idx_hyperedges[j] == h2->ti_idx_hyperedges[j]);
        INT idx_j = h->ti_idx_hyperedges[j];
        INT size_j = h-> ti_hyperedges[idx_j+1];
        assert(h->ti_hyperedges[idx_j] == h2->ti_hyperedges[idx_j]);
        assert(h->ti_hyperedges[idx_j+1] == h2->ti_hyperedges[idx_j+1]);
        for(INT i = 0; i < size_j; i++){
            assert(h->ti_hyperedges[idx_j+2+i] == h2->ti_hyperedges[idx_j+2+i]);
        }
    }

    rbhFree(h2);

    heavyEdgeMatching(h, h2, neighbors_list, in_neighbors_list, map, k, epsilon);
    
    pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

    printf("\nHEM pmax;%d\n#clusters;%d\n", pmax, h2->i_vertices);
    
    List ** neighbors_list2 = (List**)malloc(sizeof(List*)*h2->i_vertices);
    MEM_ERROR(neighbors_list2);

    for (INT i = 0; i < h2->i_vertices; i++) {
        neighbors_list2[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(neighbors_list2[i]);
        new_list(neighbors_list2[i]);
    }

    List ** in_neighbors_list2 = (List**)malloc(sizeof(List*)*h2->i_vertices);
    MEM_ERROR(in_neighbors_list2);

    for (INT i = 0; i < h2->i_vertices; i++) {
        in_neighbors_list2[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(in_neighbors_list2[i]);
        new_list(in_neighbors_list2[i]);
    }
    
    computeListNeighbors(h2, neighbors_list2);
    printf("neighbours2\n");
    computeListInNeighbors(h2, neighbors_list2, in_neighbors_list2);
    printf("in_neighbours2\n");
    
    bool * is_in = (bool*)malloc(sizeof(bool)*h->i_vertices);
    MEM_ERROR(is_in);
    
       
    
    
    computeListNeighborsUnalloc(h2, neighbors_list);
    printf("neighbours2 unalloc\n");
    computeListInNeighborsUnalloc(h2, neighbors_list, in_neighbors_list);
    printf("in_neighbours2 unalloc\n");
    
    for (INT i = 0; i < h2->i_vertices; i++) {
        
        List * cell  = neighbors_list2[i];
        List * cellu = neighbors_list[i];
        for(INT x = 0; x < neighbors_list[i]->size; x++) {
            assert(cell->i == cellu->i);
            cell  = cell->next;
            cellu = cellu->next;
        }
        cell  = in_neighbors_list2[i];
        cellu = in_neighbors_list[i];
        for(INT x = 0; x < in_neighbors_list[i]->size; x++) {
            assert(cell->i == cellu->i);
            cell  = cell->next;
            cellu = cellu->next;
        }
    }
    
    computeListNeighbors(h, base_neighbors_list);
    printf("neighbours\n");
    computeListInNeighbors(h, base_neighbors_list, base_in_neighbors_list);
    printf("in_neighbours\n");
    
    computeListNeighborsUnalloc(h, neighbors_list);
    printf("neighbours unalloc\n");
    computeListInNeighborsUnalloc(h, neighbors_list, in_neighbors_list);
    printf("in_neighbours unalloc\n");
    
    for (INT i = 0; i < h->i_vertices; i++) {
        
        List * cell  = neighbors_list[i];
        List * cellu = base_neighbors_list[i];
        for(INT x = 0; x < neighbors_list[i]->size; x++) {
            assert(cell->i == cellu->i);
            cell  = cell->next;
            cellu = cellu->next;
        }
        cell  = base_in_neighbors_list[i];
        cellu = in_neighbors_list[i];
        for(INT x = 0; x < in_neighbors_list[i]->size; x++) {
            assert(cell->i == cellu->i);
            cell  = cell->next;
            cellu = cellu->next;
        }
    }
    
    for (INT i = 0; i < h2->i_vertices; i++) {
        delete_list(neighbors_list2[i]);
        delete_list(in_neighbors_list2[i]);
    }
    
    rbhFree(h2);
    
    
    INT cluster_size = 2;
    epsilon = 5;
    INT relax = 999;
    bestPhiClustering(h, h2, a, neighbors_list, in_neighbors_list, sort, map, cluster_size, epsilon, relax);
    pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

    printf("\nBSC pmax;%d\n#clusters;%d\n", pmax, h2->i_vertices);

    free(sort);
    free(map);

    for(INT i = 0; i < h->i_vertices; i++){

        delete_list(neighbors_list[i]);
        delete_list(in_neighbors_list[i]);
        delete_list(base_neighbors_list[i]);
        delete_list(base_in_neighbors_list[i]);
        

    }

    free(neighbors_list);
    free(in_neighbors_list);
    free(base_neighbors_list);
    free(base_in_neighbors_list);
    free(neighbors_list2);
    free(in_neighbors_list2);

    free(is_in);

    rbhFree(h2);
    rbhFree(h);

    archFree(a);



    free(h2);
    free(h);


    free(arch_path);
    free(file_path);

    return 0;
}
