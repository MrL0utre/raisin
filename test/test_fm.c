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

/** 
 * @file test_fm.c
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  These lines are definitions for testing  
 *        refinement algorithms. 
 *      
 */
#include "../include/rbh.h"
#include "../include/a.h"
#include "../include/crbh.h"
#include "../include/pqueue.h"
#include "../include/list.h"
#include "../include/ipart.h"
#include "../include/fm.h"
#include "../include/dlist.h"
#include <stdlib.h>
#include <time.h>

int 
main(void)
  {
    Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h);

    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);

    Arch * a = (Arch*)malloc(sizeof(Arch));
    MEM_ERROR(a);

    char * file_path = (char*)malloc(sizeof(char) * 25);
    MEM_ERROR(file_path);

    strcpy(file_path, "../hypergraphs/b14.rzn2");

    char * arch_path = (char*)malloc(sizeof(char) * 23);
    MEM_ERROR(arch_path);

    strcpy(arch_path, "../targets/arch0.arch");

    rbhLoad(h, file_path, 0, false);

    arch_load(a, arch_path, true);

    h->s_rbh_name = (char*)malloc(sizeof(char) * strlen(file_path));
    MEM_ERROR(h->s_rbh_name);

    a->s_arch_name = (char*)malloc(sizeof(char) * strlen(arch_path));
    MEM_ERROR(a->s_arch_name);

    strcpy(h->s_rbh_name,  "b14");
    
    strcpy(a->s_arch_name, "arch0");

    rbh_save(h, h->s_rbh_name, 0, true);

    INT epsilon = 5;
    
    INT k       = 4;
    
    INT * map   = (INT*)malloc(sizeof(INT) * h->i_vertices);
    MEM_ERROR(map);
    
    INT * lambda   = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
    MEM_ERROR(lambda);
    
    List ** neighbors_list = (List**)malloc(sizeof(List*) * h->i_vertices);
    MEM_ERROR(neighbors_list);

    for (INT i = 0; i < h->i_vertices; i++) 
      {
        neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(neighbors_list[i]);
        
        new_list(neighbors_list[i]);
      }

    List ** in_neighbors_list = (List**)malloc(sizeof(List*) * h->i_vertices);
    MEM_ERROR(in_neighbors_list);

    for (INT i = 0; i < h->i_vertices; i++) 
      {
        in_neighbors_list[i] = (List*)malloc(sizeof(List));
        MEM_ERROR(in_neighbors_list[i]);
        
        new_list(in_neighbors_list[i]);
      }
    compute_list_neighbors(h, neighbors_list);
    
    printf("neighbours\n");
    
    compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
    
    printf("in_neighbours\n");
    
    INT * sort = (INT*)malloc(sizeof(INT) * h->i_vertices);
    MEM_ERROR(sort);

    INT * map2 = (INT*)malloc(sizeof(INT) * h->i_vertices);
    MEM_ERROR(map2);
    
    clock_t ta, tb;
    double time_taken;
    
    ta = clock();    
    topological_sort(h, neighbors_list, in_neighbors_list, sort);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    
    printf("Time taken by topological sort %2.f s.\n", time_taken);

    ta = clock();
    compute_criticality(h, neighbors_list, in_neighbors_list, sort);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by compute_criticality %2.f s.\n", time_taken);

    ta = clock();
    derived_breadth_first_search(h, neighbors_list, in_neighbors_list, map, k, epsilon);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by derived_breadth_first_search %2.f s.\n", time_taken);
    
    for(INT i = 0; i < h->i_vertices; i++) 
      {
        map2[i] = map[i];
      }
    int cut=0;
    
    int pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, map);   
    
    cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, map, lambda, k);

    printf("\ndbfs  pmax;%d\n#clusters;%d\n", pmax, k);
    
    printf("\ndbfs  cut;%d\n#clusters;%d\n", cut, k);
    
    INT perform = 10;
    INT tolerance = 0;
    
    ta = clock();
    dkfm(h, a, neighbors_list, in_neighbors_list, sort, map, perform, tolerance, k);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by dkfm %2.f s.\n", time_taken);    
    
    pmax =  compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, map);
    printf("\ndbfs+dkfm pmax;%d\n#clusters;%d\n", pmax, k);
    
    ta = clock();
    kfm(h, a, neighbors_list, in_neighbors_list, sort, map2, perform, tolerance, k);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by kfm %2.f s.\n", time_taken);    
    
    cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, map2, lambda, k);
    printf("\ndbfs+dkfm  cut;%d\n#clusters;%d\n", cut, k);
    
    ta = clock();
    derived_depth_first_search(h, neighbors_list, in_neighbors_list, map, k, epsilon);
    tb = clock();
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by derived_depth_first_search %2.f s.\n", time_taken);
    
    for(INT i = 0; i < h->i_vertices; i++) 
      {
        map2[i] = map[i];
      }
    pmax =  compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, map);
    
    cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, map, lambda, k); 

    printf("\nddfs pmax;%d\n#clusters;%d\n", pmax, k);
    
    printf("\nddfs  cut;%d\n#clusters;%d\n", cut, k);
    
    ta = clock();
    dkfm(h, a, neighbors_list, in_neighbors_list, sort, map, perform, tolerance, k);
    tb = clock();
    
    pmax =  compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, map); 
    printf("\nddfs+dkfm pmax;%d\n#clusters;%d\n", pmax, k);
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by dkfm %2.f s.\n", time_taken);
    
    ta = clock();
    kfm(h, a, neighbors_list, in_neighbors_list, sort, map2, perform, tolerance, k);
    tb = clock();
    
    cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, map, lambda, k);
    printf("\nddfs+dkfm  cut;%d\n#clusters;%d\n", cut, k);
    
    time_taken = ((double)tb - (double)ta) / CLOCKS_PER_SEC;
    printf("Time taken by kfm %2.f s.\n", time_taken);
     
    
    INT crit_max = 0;
    for(INT i = 0; i < h->i_vertices; i++) 
      {
        if(h->ti_criticalities_right[i] > crit_max)
          {
            crit_max = h->ti_criticalities_right[i];
          }
        map[i] = 0;
    }
    
    free(sort);
    free(map);
    free(map2);
    free(lambda);

    for(INT i = 0; i < h->i_vertices; i++)
      {
        delete_list(neighbors_list[i]);
        delete_list(in_neighbors_list[i]);
      }
    free(neighbors_list);
    free(in_neighbors_list);
    
    rbh_free(h);

    arch_free(a);

    free(h2);
    free(h);

    free(arch_path);
    free(file_path);

    return 0;
}
