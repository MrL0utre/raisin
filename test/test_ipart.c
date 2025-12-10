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


#include "../include/rbh.h"
#include "../include/a.h"
#include "../include/crbh.h"
#include "../include/pqueue.h"
#include "../include/list.h"
#include "../include/ipart.h"
#include <stdlib.h>
#include <time.h>

int main(void){

    Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h);

    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);


    Arch * a = (Arch*)malloc(sizeof(Arch));
    MEM_ERROR(a);

    char * file_path = (char*)malloc(sizeof(char)*20);
    MEM_ERROR(file_path);
    strcpy(file_path, "../hypergraphs/b14.rzn2");

    char * arch_path = (char*)malloc(sizeof(char)*20);
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

    for (INT i = 0; i < h->i_vertices; i++) {
        neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(neighbors_list[i]);
        new_list(neighbors_list[i]);
    }

    List ** in_neighbors_list = (List**)malloc(sizeof(List*)*h->i_vertices);
    MEM_ERROR(in_neighbors_list);

    for (INT i = 0; i < h->i_vertices; i++) {
        in_neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(in_neighbors_list[i]);
        new_list(in_neighbors_list[i]);
    }

    computeListNeighbors(h, neighbors_list);
    printf("neighbours\n");
    computeListInNeighbors(h, neighbors_list, in_neighbors_list);
    printf("in_neighbours\n");
    INT * sort = (INT*)malloc(sizeof(INT)*h->i_vertices);
    MEM_ERROR(sort);

    topologicalSort(h, neighbors_list, in_neighbors_list, sort);

    compute_criticality(h, neighbors_list, in_neighbors_list, sort);

    derivedBreadthFirstSearch(h, neighbors_list, in_neighbors_list, map, k, epsilon);

    int pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

    printf("\ndbfs  pmax;%d\n#clusters;%d\n", pmax, k);
    
    INT * partition_sizes = (INT*)calloc(k*h->i_weights, sizeof(INT));
    MEM_ERROR(partition_sizes);
    
     for(INT i = 0; i < h->i_vertices; i++) {
    
        for(INT wi = 0; wi < h->i_weights; wi++) {
        
            partition_sizes[map[i]*h->i_weights+wi] += h->ti_weights[i*h->i_weights+wi];
        
        }
    
    } 
    
    
    
    for(INT p = 0; p < k; p++) {
    
        printf("[pmax] final part size of %d = %d\n", p, partition_sizes[p]); 
        partition_sizes[p]=0;
    
    }   
    
    
    

    derivedDepthFirstSearch(h, neighbors_list, in_neighbors_list, map, k, epsilon);

    pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

    printf("\nddfs pmax;%d\n#clusters;%d\n", pmax, k);
   
    
    
     for(INT i = 0; i < h->i_vertices; i++) {
    
        for(INT wi = 0; wi < h->i_weights; wi++) {
        
            partition_sizes[map[i]*h->i_weights+wi] += h->ti_weights[i*h->i_weights+wi];
        
        }
    
    } 
    
    
    
    for(INT p = 0; p < k; p++) {
    
        printf("[pmax] final part size of %d = %d\n", p, partition_sizes[p]); 
        partition_sizes[p]=0;
    
    }   
    
    INT crit_max = 0;
    for(INT i = 0; i < h->i_vertices; i++) {
        if(h->ti_criticalities_right[i]> crit_max){
            crit_max = h->ti_criticalities_right[i];
        }
        map[i] = 0;
    }
    
    


    criticalConnectedComponentPartitioning(h, neighbors_list, in_neighbors_list, map, k, epsilon,  crit_max-800);

    bool * is_in = (bool*)malloc(sizeof(bool)*h->i_vertices);
    MEM_ERROR(is_in);
    
    INT * map_final = (INT*) malloc(sizeof(INT)*h->i_vertices);
    MEM_ERROR(map_final);
    INT * partition = (INT*) malloc(sizeof(INT)*h->i_vertices);
    MEM_ERROR(partition);

    for(INT i = 0; i < h->i_vertices; i++) {
        is_in[i] = false;
        partition[i] = -1;
    }

    for(INT i = 0; i < h->i_vertices; i++) {
        is_in[map[i]] = true;
    }

    INT nvp = 0;
    for(INT i = 0; i < h->i_vertices; i++) {
        if(is_in[map[i]]) {
            nvp++;
            is_in[map[i]]=false;
        }
    }
    printf("nvp : %d\n", nvp);
   
    
    INT * part_size = (INT*) calloc(nvp*h->i_weights, sizeof(INT));
    MEM_ERROR(part_size);
    
    for(INT i = 0; i < h->i_vertices; i++) {
        for(INT wi = 0; wi < h->i_weights; wi++) {
            part_size[map[i]*h->i_weights+wi] += h->ti_weights[i*h->i_weights+wi];
        }
    }
    
    for(INT i = 0; i < nvp; i++) {
        for(INT wi = 0; wi < h->i_weights; wi++) {
            
             part_size[i*h->i_weights+wi] = 0;
        }
       
    }

    pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);

    printf("\ncccp pmax1;%d\n#clusters;%d\n", pmax, nvp);
    
    compute_subhypergraph(h, h2, map, nvp);
    computeListNeighborsUnalloc(h2, neighbors_list);
    computeListInNeighborsUnalloc(h2, neighbors_list, in_neighbors_list);
    derivedBreadthFirstSearchMultilevel(h2, neighbors_list, in_neighbors_list, partition, k, epsilon);
    
    /* map[u] = v, u \in h, v \in h2*/
    for(INT i = 0; i < h->i_vertices; i++) {
        map_final[i] = partition[map[i]]; 
    }

    computeListNeighborsUnalloc(h, neighbors_list);
    computeListInNeighborsUnalloc(h, neighbors_list, in_neighbors_list);
    
    pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map_final, 800);

    printf("\ncccp+dbfs pmax1;%d\n#clusters;%d\n", pmax, k);

     for(INT i = 0; i < h->i_vertices; i++) {
    
        for(INT wi = 0; wi < h->i_weights; wi++) {
        
            partition_sizes[map_final[i]*h->i_weights+wi] += h->ti_weights[i*h->i_weights+wi];
        
        }
    
    } 
    
    
    
    for(INT p = 0; p < k; p++) {
    
        printf("[pmax] final part size of %d = %d\n", p, partition_sizes[p]); 
        partition_sizes[p]=0;
    
    }   
    
    
    free(partition_sizes);
    free(sort);
    free(map);
    free(map_final);
    free(partition);
    for(INT i = 0; i < h->i_vertices; i++){

        delete_list(neighbors_list[i]);
        delete_list(in_neighbors_list[i]);

    }
    free(neighbors_list);
    free(in_neighbors_list);
    
    rbhFree(h);
    rbhFree(h2);

    archFree(a);


    free(h2);
    free(h);


    free(arch_path);
    free(file_path);


    return 0;
}
