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
 * @file test_m.c
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief These lines are the functions    
 *        for testing the multilevel algorithms. 
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
#include "../include/m.h"
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

    char * file_path = (char*)malloc(sizeof(char)*25);
    MEM_ERROR(file_path);
    
    strcpy(file_path, "../hypergraphs/b14.rzn2");

    char * arch_path = (char*)malloc(sizeof(char)*25);
    MEM_ERROR(arch_path);
    
    strcpy(arch_path, "../targets/arch0.arch");

    rbhLoad(h, file_path, 0, false);

    arch_load(a, arch_path, true);

    h->s_rbh_name = (char*)malloc(sizeof(char)*strlen(file_path));
    MEM_ERROR(h->s_rbh_name);

    a->s_arch_name = (char*)malloc(sizeof(char)*strlen(arch_path));
    MEM_ERROR(a->s_arch_name);

    strcpy(h->s_rbh_name,  "b14");
    
    strcpy(a->s_arch_name, "arch0");

    rbh_save(h, h->s_rbh_name, 0, true);

    INT epsilon = 5;
    INT k       = 4;
    
    INT * partition   = (INT*)malloc(sizeof(INT ) * h->i_vertices);
    MEM_ERROR(partition);
    
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
    
    INT perform = 10;
    
    INT tolerance = 0;

    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    INT div    = nv / 10; 
    INT levels = log(nv);/// 2 + 1;    
    
    /**  **/
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*) * levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*) * (levels - 1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT) * nv);
    MEM_ERROR(map_final);
    
    INT * partitionp = (INT*)malloc(sizeof(INT) * nv);
    MEM_ERROR(partitionp);
    
    INT * partition_sizes = (INT*)calloc(k * h->i_weights, sizeof(INT));
    MEM_ERROR(partition_sizes);
    
    INT * partition_crit = (INT*)malloc(sizeof(INT) * nv);
    MEM_ERROR(partition_crit);
    
    INT * partition_critp = (INT*)malloc(sizeof(INT) * nv);
    MEM_ERROR(partition_critp);
    
    hypergraphs[0] = h; 
        
    for(INT i = 0; i < nv; i++) 
      {
        map_final[i] = i;
      }
  
    for(INT i = 1; i < levels; i++) 
      {  
        printf("level %d\n", i);
                        
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
        
        printf("nv %d\n", hypergraphs[i - 1]->i_vertices);
      
        maps[i - 1] = (INT*)malloc(sizeof(INT)*hypergraphs[i - 1]->i_vertices);
        MEM_ERROR(maps[i - 1]);
        
        heavy_edge_matching(hypergraphs[i - 1], hypergraphs[i], neighbors_list, in_neighbors_list, maps[i - 1], k, epsilon);

        for(INT u = 0; u < h->i_vertices; u++) 
          {
            map_final[u] = maps[i - 1][map_final[u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i], neighbors_list);
        
        compute_list_in_neighbors_unalloc(hypergraphs[i], neighbors_list, in_neighbors_list);       
      }
    
    derived_depth_first_searchMultilevel(hypergraphs[levels - 1], neighbors_list, in_neighbors_list, partition, k, epsilon);
    
    memcpy(partition_crit, partition, sizeof(INT) * h->i_vertices);
    
    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        partition_crit[i] = partition[i];
      }
    
    INT cut_cost = compute_partition_cut(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
    
    topological_sort(hypergraphs[levels - 1], neighbors_list, in_neighbors_list, sort);
    
    INT pmax = compute_partition_criticality(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition_crit);
    
    for(INT p = 0; p < k; p++) 
      {
        partition_sizes[p] = 0;
      }

    printf("dbfs cut cost = %d\n", cut_cost); 
    
    printf("dbfs pmax cost = %d\n", pmax); 
    
    kfm(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k);
    
    dkfm(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition_crit, perform, tolerance, k);
    
    cut_cost = compute_partition_cut(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
    
    pmax = compute_partition_criticality(hypergraphs[levels - 1], a, neighbors_list, in_neighbors_list, sort, partition_crit);
    
    printf(" cut cost = %d\n", cut_cost); 
    printf(" pmax cost = %d\n", pmax); 

    for(INT i = levels - 1; i > 0; i--) 
      {  
        printf("level %d\n", i);
                
        printf("nv %d\n", hypergraphs[i - 1]->i_vertices);
        
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partitionp[u] = partition[maps[i - 1][u]];
            
            partition_critp[u] = partition_crit[maps[i - 1][u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i - 1], neighbors_list);
        
        compute_list_in_neighbors_unalloc(hypergraphs[i - 1], neighbors_list, in_neighbors_list);
        
        topological_sort(hypergraphs[i - 1], neighbors_list, in_neighbors_list, sort);
        
        kfm(hypergraphs[i - 1], a, neighbors_list, in_neighbors_list, sort, partitionp, perform, tolerance, k); 
        
        dkfm(hypergraphs[i - 1], a, neighbors_list, in_neighbors_list, sort, partition_critp, perform, tolerance, k);
        
        cut_cost = compute_partition_cut(hypergraphs[i - 1], a, neighbors_list, in_neighbors_list, sort, partitionp, lambda, k); 
        
           pmax = compute_partition_criticality(hypergraphs[i - 1], a, neighbors_list, in_neighbors_list, sort, partition_critp);
        
        printf("cut cost = %d\n", cut_cost);
        printf(" pmax cost = %d\n", pmax);       
        
        for(INT u = 0; u < hypergraphs[0]->i_vertices; u++) 
          {  
            partition[u] = partitionp[u];
            partition_crit[u] = partition_critp[u];
          } 
      }
    cut_cost = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
    
    pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition_crit);
    
    for(INT i = 0; i < hypergraphs[0]->i_vertices; i++) 
      {
        for(INT wi = 0; wi < hypergraphs[0]->i_weights; wi++) 
          {
            partition_sizes[partition[i] * hypergraphs[0]->i_weights + wi] += hypergraphs[0]->ti_weights[i * hypergraphs[0]->i_weights + wi];
          }
      }
    
    for(INT p = 0; p < k; p++) 
      {
        printf("[cut] final part size of %d = %d\n", p, partition_sizes[p]); 
        
        partition_sizes[p] = 0;
      }   
    
    for(INT i = 0; i < hypergraphs[0]->i_vertices; i++) 
      {
        for(INT wi = 0; wi < hypergraphs[0]->i_weights; wi++) 
          {
            partition_sizes[partition_crit[i] * hypergraphs[0]->i_weights + wi] += hypergraphs[0]->ti_weights[i * hypergraphs[0]->i_weights + wi];
          }
      } 
    
    printf("final cut cost = %d\n", cut_cost);
    
    for(INT p = 0; p < k; p++) 
      {
        printf("[pmax] final part size of %d = %d\n", p, partition_sizes[p]); 
        
        partition_sizes[p] = 0;
      }   
    printf("final pmax cost = %d\n", pmax);
     
    multilevel_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, 0);	
     
    cut_cost = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
      
    printf("final multilevel function cut cost = %d\n", cut_cost);
      
    multilevel_pmax(h, a, neighbors_list, in_neighbors_list, sort, partition_crit, perform, tolerance, k, epsilon, 0);
     
    pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition_crit);
      
    printf("final multilevel function pmax cost = %d\n", pmax);
    
    /* free section */
    for(INT level = 1; level < levels; level++) 
      {   
        rbh_free(hypergraphs[level]);
        
        free(hypergraphs[level]);            
        
        free(maps[level - 1]);
      }
    
    free(partition_crit);
    free(partition_critp);
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);
    free(partition_sizes);
    free(sort);
    free(partition);
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
