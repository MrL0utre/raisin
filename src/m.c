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
/**   NAME       : m.c                                     **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are declarations for the    **/
/**                multi-level schemes.                    **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 01 jan 2022     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


/*
**  The defines and includes.
*/
#include "m.h"

INT multilevel_cut(
Hypergraph  * h,
Arch        * a,
List       ** neighbors, 
List       ** in_neighbors, 
INT         * sort, 
INT         * partition, 
INT           perform, 
INT           tolerance,
INT           k,
INT           epsilon,
INT           algo_cluster) {

    bool verbose = false;
    
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    INT relax  = 999;
    INT div    = nv / 10; 
    INT levels = log(nv); 
    
    if (algo_cluster==10) {
        levels /= 2;
    }   
    
    /**  **/
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*)*levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*)*(levels-1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    INT * partitionp = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(partitionp);
    
    INT * lambda   = (INT*)malloc(sizeof(INT)*h->i_hyperedges);
    MEM_ERROR(lambda);
    
    
    hypergraphs[0]    = h; 
        
    for(INT i = 0; i < nv; i++) {
    
        map_final[i] = i;
    
    }

        
    for(INT i = 1; i < levels; i++) {

        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
        
        maps[i-1]   = (INT*)malloc(sizeof(INT)*hypergraphs[i-1]->i_vertices);
        MEM_ERROR(maps[i-1]);
        
        if (algo_cluster==0 || algo_cluster==10) {
        
            heavyEdgeMatching(hypergraphs[i-1], hypergraphs[i], neighbors, in_neighbors, maps[i-1], k, epsilon);
        
        }
        
        if (algo_cluster==1 || algo_cluster==11) {
        
            bestPhiClustering(hypergraphs[i-1], hypergraphs[i], a, neighbors, in_neighbors, sort, maps[i-1], levels, epsilon, relax);

        }

        for(INT u = 0; u < h->i_vertices; u++) {
            
            map_final[u] = maps[i-1][map_final[u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i], neighbors, in_neighbors);
        
               
    }
    
    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);
           
    INT crit_max = 0;
    for(INT i = 0; i < h->i_vertices; i++) {
         if(h->ti_criticalities_right[i]> crit_max){
              crit_max = h->ti_criticalities_right[i];
         }
    }
    
    INT * map = (INT*)malloc(sizeof(INT)*hypergraphs[levels-1]->i_vertices);
    MEM_ERROR(map);
           
    INT D = ceil((float)crit_max*0.8); // default value ~ 20 % of CP 
    
    criticalConnectedComponentPartitioningMultilevel(hypergraphs[levels-1], neighbors, in_neighbors, map, k, epsilon,  crit_max-D);

    bool * is_in = (bool*) malloc(sizeof(bool)*hypergraphs[levels-1]->i_vertices);
    MEM_ERROR(is_in);

    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        is_in[i] = false;
    }

    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        is_in[map[i]] = true;
    }

    INT nvp = 0;
    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        if(is_in[map[i]]) {
            nvp++;
            is_in[map[i]]=false;
        }
    }      
           
    compute_subhypergraph(hypergraphs[levels-1], h2, map, nvp);
    computeListNeighborsUnalloc(h2, neighbors);
    computeListInNeighborsUnalloc(h2, neighbors, in_neighbors);
    derivedBreadthFirstSearchMultilevel(h2, neighbors, in_neighbors, partitionp, k, epsilon);  
    
    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        partition[i] = partitionp[map[i]];
    }
    
    computeListNeighborsUnalloc(hypergraphs[levels-1], neighbors);
    computeListInNeighborsUnalloc(hypergraphs[levels-1], neighbors, in_neighbors);
    
    INT cut_cost = compute_partition_cut(hypergraphs[levels-1], a, neighbors, in_neighbors, sort, partition, lambda, k);
    
    kfm(hypergraphs[levels-1], a, neighbors, in_neighbors, sort, partition, perform, tolerance, k);

    for(INT i = levels-1; i > 0; i--) {
              
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partitionp[u] = partition[maps[i-1][u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i-1], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i-1], neighbors, in_neighbors);
        
        kfm(hypergraphs[i-1], a, neighbors, in_neighbors, sort, partitionp, perform, tolerance, k);        
        
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partition[u] = partitionp[u];
        } 
        
    }   
    

    /* free section */
    for(INT level = 1; level < levels; level++) {
        
        
        rbhFree(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level-1]);
    
    }
    free(lambda);
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);
    free(map);
    free(is_in);
    
    
}


INT multilevel_pmax(
Hypergraph  * h,
Arch        * a,
List       ** neighbors, 
List       ** in_neighbors, 
INT         * sort, 
INT         * partition, 
INT           perform, 
INT           tolerance,
INT           k,
INT           epsilon,
INT           algo_cluster) {

    bool verbose = false;
    
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    INT relax  = 999;
    INT div    = nv / 10; 
    INT levels = log(nv);    
    
    if (algo_cluster==10) {
        levels /= 2;
    }
    
    /**  **/
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*)*levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*)*(levels-1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    INT * partitionp = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(partitionp);   
    
    hypergraphs[0]    = h; 
        
    for(INT i = 0; i < nv; i++) {
    
        map_final[i] = i;
    
    }

        
    for(INT i = 1; i < levels; i++) {
                
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
      
        maps[i-1]   = (INT*)malloc(sizeof(INT)*hypergraphs[i-1]->i_vertices);
        MEM_ERROR(maps[i-1]);
        
        if (algo_cluster==0 || algo_cluster==10) {
        
            heavyEdgeMatching(hypergraphs[i-1], hypergraphs[i], neighbors, in_neighbors, maps[i-1], k, epsilon);
        
        }
        
        if (algo_cluster==1 || algo_cluster==11) {
        
            bestPhiClustering(hypergraphs[i-1], hypergraphs[i], a, neighbors, in_neighbors, sort, maps[i-1], levels, epsilon, relax);

        }

        for(INT u = 0; u < h->i_vertices; u++) {
            
            map_final[u] = maps[i-1][map_final[u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i], neighbors, in_neighbors);
        
               
    }
    
    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);
           
    INT crit_max = 0;
    for(INT i = 0; i < h->i_vertices; i++) {
         if(h->ti_criticalities_right[i]> crit_max){
              crit_max = h->ti_criticalities_right[i];
         }
    }
    
    INT * map = (INT*)malloc(sizeof(INT)*hypergraphs[levels-1]->i_vertices);
    MEM_ERROR(map);
           
    INT D = ceil((float)crit_max*0.8); // default value ~ 20 % of CP 
    
    criticalConnectedComponentPartitioningMultilevel(hypergraphs[levels-1], neighbors, in_neighbors, map, k, epsilon,  crit_max-D);

    bool * is_in = (bool*) malloc(sizeof(bool)*hypergraphs[levels-1]->i_vertices);
    MEM_ERROR(is_in);

    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        is_in[i] = false;
    }

    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        is_in[map[i]] = true;
    }

    INT nvp = 0;
    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        if(is_in[map[i]]) {
            nvp++;
            is_in[map[i]]=false;
        }
    }      
           
    compute_subhypergraph(hypergraphs[levels-1], h2, map, nvp);
    computeListNeighborsUnalloc(h2, neighbors);
    computeListInNeighborsUnalloc(h2, neighbors, in_neighbors);
    derivedBreadthFirstSearchMultilevel(h2, neighbors, in_neighbors, partitionp, k, epsilon);  
    
    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) {
        partition[i] = partitionp[map[i]];
    }
    
    computeListNeighborsUnalloc(hypergraphs[levels-1], neighbors);
    computeListInNeighborsUnalloc(hypergraphs[levels-1], neighbors, in_neighbors);
    
    topologicalSort(hypergraphs[levels-1], neighbors, in_neighbors, sort);
    
    INT pmax = compute_partition_criticality(h, a, neighbors, in_neighbors, sort, partition, k);
    
    if (algo_cluster==11 || algo_cluster==10) {
        dkfmFast(hypergraphs[levels-1], a, neighbors, in_neighbors, sort, partition, perform, tolerance, k);
    } else {
        dkfm(hypergraphs[levels-1], a, neighbors, in_neighbors, sort, partition, perform, tolerance, k);
    }
    
    

    for(INT i = levels-1; i > 0; i--) {
    
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partitionp[u] = partition[maps[i-1][u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i-1], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i-1], neighbors, in_neighbors);
        
        topologicalSort(hypergraphs[i-1], neighbors, in_neighbors, sort);
        if (algo_cluster==11 || algo_cluster==10) {
            dkfmFast(hypergraphs[i-1], a, neighbors, in_neighbors, sort, partitionp, perform, tolerance, k); 
        } else {
            dkfm(hypergraphs[i-1], a, neighbors, in_neighbors, sort, partitionp, perform, tolerance, k); 
        }      
               
        
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partition[u] = partitionp[u];
        } 
        
    }


    /* free section */
    for(INT level = 1; level < levels; level++) {
        
        
        rbhFree(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level-1]);
    
    }
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);
    
    
}

INT multilevel(
Hypergraph  * h,
Arch        * a,
List       ** neighbors, 
List       ** in_neighbors, 
INT         * sort, 
INT         * partition, 
INT           perform, 
INT           tolerance,
INT           k,
INT           epsilon,
INT           algo_cluster,
int (*intialPart)(Hypergraph * h, List ** neighbors, List ** in_neighbors, INT * partition, INT k, INT epsilon),
int (*refinement)(Hypergraph * h, Arch * a, List ** neighbors, List ** in_neighbors, INT * sort, INT * partition, INT perform, INT tolerance, INT k)) {

    bool verbose = false;
    
    
    INT nv     = h->i_vertices;
    INT ne     = h->i_hyperedges;
    INT nr     = h->i_reds;
    INT np     = h->i_pins;
    INT nwv    = h->i_weights;
    INT * reds = h->ti_reds;
    
    INT relax  = 999;
    INT div    = nv / 10; 
    INT levels = log(nv);   
    
    if (algo_cluster==10) {
        levels /= 2;
    } 
    
    /**  **/
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*)*levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*)*(levels-1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    INT * partitionp = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(partitionp);
    
    INT * lambda   = (INT*)malloc(sizeof(INT)*h->i_hyperedges);
    MEM_ERROR(lambda);
    
    
    
    hypergraphs[0]    = h; 
        
    for(INT i = 0; i < nv; i++) {
    
        map_final[i] = i;
    
    }

        
    for(INT i = 1; i < levels; i++) {
    
         
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
      
        maps[i-1]   = (INT*)malloc(sizeof(INT)*hypergraphs[i-1]->i_vertices);
        MEM_ERROR(maps[i-1]);
        
        if (algo_cluster==0 || algo_cluster==10) {
        
            heavyEdgeMatching(hypergraphs[i-1], hypergraphs[i], neighbors, in_neighbors, maps[i-1], k, epsilon);
        
        }
        
        if (algo_cluster==11 || algo_cluster==1) {
        
            bestPhiClustering(hypergraphs[i-1], hypergraphs[i], a, neighbors, in_neighbors, sort, maps[i-1], levels, epsilon, relax);

        }
        
        

        for(INT u = 0; u < h->i_vertices; u++) {
            
            map_final[u] = maps[i-1][map_final[u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i], neighbors, in_neighbors);
        
               
    }
    
    intialPart(hypergraphs[levels-1], neighbors, in_neighbors, partition, k, epsilon);
    
    topologicalSort(hypergraphs[levels-1], neighbors, in_neighbors, sort);
    
    refinement(hypergraphs[levels-1], a, neighbors, in_neighbors, sort, partition, perform, tolerance, k);

    for(INT i = levels-1; i > 0; i--) {
        
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partitionp[u] = partition[maps[i-1][u]];
        }      
        
        computeListNeighborsUnalloc(hypergraphs[i-1], neighbors);
        computeListInNeighborsUnalloc(hypergraphs[i-1], neighbors, in_neighbors);
        topologicalSort(hypergraphs[i-1], neighbors, in_neighbors, sort);
        refinement(hypergraphs[i-1], a, neighbors, in_neighbors, sort, partitionp, perform, tolerance, k);        
        
        for(INT u = 0; u < hypergraphs[i-1]->i_vertices; u++) {
            
            partition[u] = partitionp[u];
        } 
        
    }   
    

    /* free section */
    for(INT level = 1; level < levels; level++) {
        
        
        rbhFree(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level-1]);
    
    }
    free(lambda);
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);
    
    
}
