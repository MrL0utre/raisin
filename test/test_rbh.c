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
 * @file test_rbh.c
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief These lines are the declaration of a main 
 *   to test the red-black hypergraph functions: 
 *       - rbh_save, 
 *       - compute_list_neighbors,
 *       - compute_list_in_neighbors,
 *       - topological_sort,
 *       - compute_criticality
 *       - compute_clustering_criticality.
 */
#include"../include/rbh.h"
#include"../include/a.h"
#include"../include/crbh.h"
#include"../include/pqueue.h"
#include"../include/list.h"

#include <stdlib.h>
#include <time.h>

int 
main(void)
{
    Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h);

    Arch * a = (Arch*)malloc(sizeof(Arch));
    MEM_ERROR(a);

    char * file_path = (char*)malloc(sizeof(char) * 25);
    MEM_ERROR(file_path);

    strcpy(file_path, "../hypergraphs/b14.rzn2");

    char * arch_path = (char*)malloc(sizeof(char) * 22);
    MEM_ERROR(arch_path);

    strcpy(arch_path, "../targets/arch0.arch");

    h->s_rbh_name = (char*)malloc(sizeof(char) * strlen(file_path));
    MEM_ERROR(h->s_rbh_name);

    a->s_arch_name = (char*)malloc(sizeof(char) * strlen(arch_path));
    MEM_ERROR(a->s_arch_name);

    strcpy(h->s_rbh_name,  "b14");
    strcpy(a->s_arch_name, "arch0");

    rbhLoad(h, file_path, 0, false);

    arch_load(a, arch_path, true);

    rbh_save(h, h->s_rbh_name, 0, true);

    INT * map   = (INT*)malloc(sizeof(INT) * h->i_vertices);
    MEM_ERROR(map);

    List ** neighbors_list = (List**)malloc(sizeof(List*) * h->i_vertices);
    MEM_ERROR(neighbors_list);

    for(INT i = 0; i < h->i_vertices; i++) 
      {
        neighbors_list[i] = (List*)malloc(sizeof(List));
        
        MEM_ERROR(neighbors_list[i]);
        
        new_list(neighbors_list[i]);
      }

    List ** in_neighbors_list = (List**)malloc(sizeof(List*) * h->i_vertices);
    MEM_ERROR(in_neighbors_list);

    for(INT i = 0; i < h->i_vertices; i++) 
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

    /* test in neighbour */
    for(INT u = 0; u < h->i_vertices; u++)
      {
        bool is_in;
        List * cell = neighbors_list[u];
        
        while (cell != NULL && neighbors_list[u]->size > 0) 
          {
            /* arc u --> v */
            /* for each v, \exist u in \Gamma^-(v)*/
            INT v = cell->i;
            
            is_in = false;
            
            List * in_cell = in_neighbors_list[v];
            
            while (in_cell != NULL && in_neighbors_list[v]->size > 0) 
              {
                /* arc u --> v */
                INT w = in_cell->i;
                
                if (w == u)
                  {
                    is_in = true;
                  }
                in_cell = in_cell->next;
            }
            cell = cell->next;
            
            assert(is_in);
          }
      }
    topological_sort(h, neighbors_list, in_neighbors_list, sort);

    /* test topological sort */
    bool * visited = (bool*)malloc(sizeof(bool) * h->i_vertices);
    MEM_ERROR(visited);

    bool * is_red  = (bool*)malloc(sizeof(bool) * h->i_vertices);
    MEM_ERROR(is_red);

    for(INT i = 0; i < h->i_vertices; i++)
      { 
        visited[i] = false;
        is_red[i] = false; 
      }

    for(INT i = 0; i < h->i_reds; i++)
      {
        is_red[h->ti_reds[i]] = true;
      }
    
    for(INT i = 0; i < h->i_vertices; i++)
      {
        INT u = sort[i];
        
        visited[u] = true;
        
        List * cell = in_neighbors_list[u];
        
        while (cell != NULL && in_neighbors_list[u]->size > 0) 
          {
            INT v = cell->i;

            assert(visited[v] || is_red[v]);

            cell = cell->next;
          }
        cell = neighbors_list[u];
        
        while (cell != NULL && neighbors_list[u]->size > 0) 
          {
            INT v = cell->i;
            
            if(visited[v] && !is_red[v])
              {
                printf("visited : (%d, %d), %d\n", u, v, v);
                
                printf("problem : (%d, %d), %d\n", u, v, v);
                    
                List * aux_cell = neighbors_list[v];
                            
                while (aux_cell != NULL && neighbors_list[u]->size > 0) 
                  {
                    INT w = aux_cell->i;
                                
                    printf(" %d", w);
                                
                    aux_cell = aux_cell->next;
                  }
                printf("\n");
              }
            if(is_red[v] || !visited[v]) 
              {
                printf("This hypergraph contains multiple unconnected components (%d).\n", v);
              }
                    
            cell = cell->next;
          }
      }

    free(visited);
    
    free(is_red);

    compute_criticality(h, neighbors_list, in_neighbors_list, sort);

    INT * partition = (INT*) malloc(sizeof(INT) * h->i_vertices);

    for(INT i = 0; i < h->i_vertices; i++) 
      {
        partition[i] = 0;
      }

    int pmaxa = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, partition, 800);
    
    for(INT i = 0; i < h->i_vertices; i++) 
      {
        partition[i] = i;
      }
    
    int pmaxb = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, partition, 800);

    printf("static pmax : %d\n, max pmax : %d\n " , pmaxa, pmaxb);
    
    INT * pqueue = (INT*)malloc(sizeof(INT) * h->i_vertices);

    MEM_ERROR(pqueue);
    
    INT size = 0;
    
    INT x;

    for(INT i = 0; i < h->i_vertices; i++)
      {
        x = rand() % 1000;
        
        pqueue_add_element(pqueue, h->ti_criticalities_right, i, size);
        
        size++;

        map[i] = i;
    }
    srand(time(NULL));

    for(INT i = 0; i < 10; i++)
      {
        printf("heap[%d] : h->ti_criticalities_right[%d](%d) / delay : %d\n", i, pqueue[i], h->ti_criticalities_right[pqueue[i]], h->ti_delays[pqueue[i]]);
      }

    free(sort);
    free(pqueue);

    for(INT i = 0; i < h->i_vertices; i++)
      {
        delete_list(neighbors_list[i]);
        delete_list(in_neighbors_list[i]);
      }

    free(neighbors_list);
    free(in_neighbors_list);

    rbh_free(h);

    arch_free(a);

    free(h);

    free(map);
    free(arch_path);
    free(file_path);

    return 0;
}
