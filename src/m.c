/*
 * Copyright (C) 2025 Commissariat à l'énergie atomique et aux énergies
 * alternatives (CEA) and Institut national de recherche en sciences et
 * technologies du numérique (INRIA)
 * Contributor(s): Julien Rodriguez <julien.ro34@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
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
#include "rbh.h"
#include "m.h"

/**
 * @brief Function implementing a multilevel schemme approach using
 * coarsening for red-black hypergraph defined in crbh.c, 
 * partitioning algorithms defined in ipart.c and refinement algorithm
 * defined in fm.c for minimizing the cut cost.  
 *
 * @param h                Input hypergraph.
 * @param a                Architecture (target topology).
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Partition array.
 * @param perform          Percentage of moves. 
 * @param tolerance        Tolerance for cost degradation. 
 * @param k                Number of part. 
 * @param epsilon          Tolerance for balancing cost. 
 * @param algo_cluster     Clustering algorithm.
 *
 * @return Integer.
 */
INT 
multilevel_cut(Hypergraph  * h,
               Arch        * a,
               List       ** out_neighbors, 
               List       ** in_neighbors, 
               INT         * sort, 
               PART * partition, 
               INT           perform, 
               INT           tolerance,
               INT           k,
               INT           epsilon,
               INT           algo_cluster) 
{
    INT nv     = h->i_vertices;    /* number of vertices                        */

    INT relax  = 999;
    INT levels = log(nv); 
    
    if(algo_cluster==10) {
        levels /= 2;
    }   
    
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*) * levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*) * (levels - 1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    PART * partitionp = (PART*)malloc(sizeof(PART)*nv);
    MEM_ERROR(partitionp);
    
    INT * lambda   = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
    MEM_ERROR(lambda);
    
    hypergraphs[0]    = h; 
        
    for(INT i = 0; i < nv; i++) 
      {
        map_final[i] = i;
      }
        
    for(INT i = 1; i < levels; i++) 
      {
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
        
        hypergraphs[i]->is_red = NULL; hypergraphs[i]->vth = NULL;
        
        maps[i - 1] = (INT*)malloc(sizeof(INT) * hypergraphs[i - 1]->i_vertices);
        MEM_ERROR(maps[i - 1]);
        
        if(algo_cluster == 0 || algo_cluster == 10) 
          {
            heavy_edge_matching(hypergraphs[i - 1], hypergraphs[i], out_neighbors, in_neighbors, maps[i - 1], k, epsilon);        
          }
        
        if(algo_cluster == 1 || algo_cluster == 11) 
          {
            best_phi_clustering(hypergraphs[i - 1], hypergraphs[i], a, out_neighbors, in_neighbors, sort, maps[i - 1], levels, epsilon, relax);
            rbh_build_precomputed(hypergraphs[i]);
          }

        for(INT u = 0; u < h->i_vertices; u++) 
          {  
            map_final[u] = maps[i - 1][map_final[u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i], out_neighbors);
        compute_list_in_neighbors_unalloc(hypergraphs[i], out_neighbors, in_neighbors);         
      }
    
    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);           
    INT crit_max = 0;

    for(INT i = 0; i < h->i_vertices; i++) 
      {
         if(h->ti_criticalities_right[i] > crit_max)
           {
              crit_max = h->ti_criticalities_right[i];
           }
      }
    
    INT * map = (INT*)malloc(sizeof(INT) * hypergraphs[levels - 1]->i_vertices);
    MEM_ERROR(map);
           
    INT D = ceil((float)crit_max * 0.8); // default value ~ 20 % of CP 
    
    critical_connected_component_partitioning_multilevel(hypergraphs[levels - 1], out_neighbors, in_neighbors, map, k, epsilon,  crit_max-D);

    bool * is_in = (bool*)malloc(sizeof(bool) * hypergraphs[levels - 1]->i_vertices);
    MEM_ERROR(is_in);

    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        is_in[i] = false;
      }

    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        is_in[map[i]] = true;
      }

    INT nvp = 0;
    
    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        if(is_in[map[i]]) 
          {
            nvp++;
            is_in[map[i]]=false;
          }
      }      
           
    compute_subhypergraph(hypergraphs[levels-1], h2, map, nvp);
    
    compute_list_neighbors_unalloc(h2, out_neighbors);
    
    compute_list_in_neighbors_unalloc(h2, out_neighbors, in_neighbors);
    
    derived_breadth_first_search_multilevel(h2, out_neighbors, in_neighbors, partitionp, k, epsilon);  
    
    for(INT i = 0; i < hypergraphs[levels-1]->i_vertices; i++) 
      {
        partition[i] = partitionp[map[i]];
      }
    
    compute_list_neighbors_unalloc(hypergraphs[levels - 1], out_neighbors);
    
    compute_list_in_neighbors_unalloc(hypergraphs[levels - 1], out_neighbors, in_neighbors);
    
    INT cut_cost = compute_partition_cut(hypergraphs[levels - 1], a, out_neighbors, in_neighbors, sort, partition, lambda, k);
    
    kfm(hypergraphs[levels - 1], a, out_neighbors, in_neighbors, sort, partition, perform, tolerance, k);

    for(INT i = levels - 1; i > 0; i--) 
      {       
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partitionp[u] = partition[maps[i - 1][u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i - 1], out_neighbors);
        
        compute_list_in_neighbors_unalloc(hypergraphs[i - 1], out_neighbors, in_neighbors);
        
        kfm(hypergraphs[i - 1], a, out_neighbors, in_neighbors, sort, partitionp, perform, tolerance, k);        
        
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partition[u] = partitionp[u];
          }  
      }   

    /* free section */
    for(INT level = 1; level < levels; level++) 
      {        
        rbh_free(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level - 1]);
    }

    free(lambda);
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);
    free(map);
    free(is_in);

    return cut_cost;
}

/**
 * @brief Function implementing a multilevel schemme approach using
 * coarsening for red-black hypergraph defined in crbh.c, 
 * partitioning algorithms defined in ipart.c and refinement algorithm
 * defined in fm.c for minimizing the degradation of the critical cost.
 *
 * @param h                Input hypergraph.
 * @param a                Architecture (target topology).
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Partition array.
 * @param perform          Percentage of moves. 
 * @param tolerance        Tolerance for cost degradation. 
 * @param k                Number of part. 
 * @param epsilon          Tolerance for balancing cost. 
 * @param algo_cluster     Clustering algorithm.
 *
 * @return Integer.
 */
INT 
multilevel_pmax(Hypergraph  * h,
                Arch        * a,
                List       ** out_neighbors, 
                List       ** in_neighbors, 
                INT         * sort, 
                PART * partition, 
                INT           perform, 
                INT           tolerance,
                INT           k,
                INT           epsilon,
                INT           algo_cluster) 
{
    INT nv     = h->i_vertices;    /* number of vertices                        */

    INT relax  = 999;
    INT levels = log(nv);    
    
    if(algo_cluster == 10) 
      {
        levels /= 2;
      }
    
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*) * levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*) * (levels - 1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    PART * partitionp = (PART*)malloc(sizeof(PART)*nv);
    MEM_ERROR(partitionp);   
    
    hypergraphs[0] = h; 
        
    for(INT i = 0; i < nv; i++) 
      {
        map_final[i] = i;
      }
  
    for(INT i = 1; i < levels; i++) 
      {          
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
      
        hypergraphs[i]->is_red = NULL; hypergraphs[i]->vth = NULL;
      
        maps[i - 1] = (INT*)malloc(sizeof(INT) * hypergraphs[i - 1]->i_vertices);
        MEM_ERROR(maps[i - 1]);
        
        if(algo_cluster == 0 || algo_cluster == 10) 
          {
            heavy_edge_matching(hypergraphs[i - 1], hypergraphs[i], out_neighbors, in_neighbors, maps[i - 1], k, epsilon);
        rbh_build_precomputed(hypergraphs[i]);
          }
        
        if(algo_cluster == 1 || algo_cluster == 11) 
          {
            best_phi_clustering(hypergraphs[i - 1], hypergraphs[i], a, out_neighbors, in_neighbors, sort, maps[i - 1], levels, epsilon, relax);
            rbh_build_precomputed(hypergraphs[i]);
          }

        for(INT u = 0; u < h->i_vertices; u++) 
          {  
            map_final[u] = maps[i - 1][map_final[u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i], out_neighbors);
        
        compute_list_in_neighbors_unalloc(hypergraphs[i], out_neighbors, in_neighbors);       
    }
    
    Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
    MEM_ERROR(h2);
           
    h2->is_red = NULL; h2->vth = NULL;
           
    INT crit_max = 0;

    for(INT i = 0; i < h->i_vertices; i++) 
      {
         if(h->ti_criticalities_right[i] > crit_max)
           {
              crit_max = h->ti_criticalities_right[i];
           }
      }
    
    INT * map = (INT*)malloc(sizeof(INT) * hypergraphs[levels - 1]->i_vertices);
    MEM_ERROR(map);
           
    INT D = ceil((float)crit_max * 0.8); // default value ~ 20 % of CP 
    
    critical_connected_component_partitioning_multilevel(hypergraphs[levels -1 ], out_neighbors, in_neighbors, map, k, epsilon,  crit_max - D);

    bool * is_in = (bool*)malloc(sizeof(bool)*hypergraphs[levels - 1]->i_vertices);
    MEM_ERROR(is_in);

    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        is_in[i] = false;
      }

    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        is_in[map[i]] = true;
      }

    INT nvp = 0;
    
    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        if(is_in[map[i]]) 
          {
            nvp++;
            is_in[map[i]] = false;
          }
      }      
           
    compute_subhypergraph(hypergraphs[levels-1], h2, map, nvp);
    
    compute_list_neighbors_unalloc(h2, out_neighbors);
    
    compute_list_in_neighbors_unalloc(h2, out_neighbors, in_neighbors);
    
    derived_breadth_first_search_multilevel(h2, out_neighbors, in_neighbors, partitionp, k, epsilon);  
    
    for(INT i = 0; i < hypergraphs[levels - 1]->i_vertices; i++) 
      {
        partition[i] = partitionp[map[i]];
      }
    
    compute_list_neighbors_unalloc(hypergraphs[levels - 1], out_neighbors);
    
    compute_list_in_neighbors_unalloc(hypergraphs[levels - 1], out_neighbors, in_neighbors);
    
    topological_sort(hypergraphs[levels - 1], out_neighbors, in_neighbors, sort);
    
    INT pmax = compute_partition_criticality(h, a, out_neighbors, in_neighbors, sort, partition);
    
    if(algo_cluster == 11 || algo_cluster == 10) 
      {
        dkfm_fast(hypergraphs[levels - 1], a, out_neighbors, in_neighbors, sort, partition, perform, tolerance, k);
      } 
    else 
      {
        dkfm(hypergraphs[levels - 1], a, out_neighbors, in_neighbors, sort, partition, perform, tolerance, k);
      }
    
    for(INT i = levels - 1; i > 0; i--) 
      {
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partitionp[u] = partition[maps[i - 1][u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i - 1], out_neighbors);
        
        compute_list_in_neighbors_unalloc(hypergraphs[i - 1], out_neighbors, in_neighbors);
        
        topological_sort(hypergraphs[i - 1], out_neighbors, in_neighbors, sort);
        
        if(algo_cluster == 11 || algo_cluster == 10) 
          {
            dkfm_fast(hypergraphs[i - 1], a, out_neighbors, in_neighbors, sort, partitionp, perform, tolerance, k); 
          } 
        else 
          {
            dkfm(hypergraphs[i - 1], a, out_neighbors, in_neighbors, sort, partitionp, perform, tolerance, k); 
          }      
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partition[u] = partitionp[u];
          } 
      }

    /* free section */
    for(INT level = 1; level < levels; level++) 
      {  
        rbh_free(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level - 1]);
      }
    
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);

    return pmax;
}

/**
 * @brief Function implementing a multilevel schemme approach using
 * coarsening for red-black hypergraph defined in crbh.c, 
 * partitioning algorithms defined in ipart.c and refinement algorithm
 * defined in fm.c.
 *
 * @param h                Input hypergraph.
 * @param a                Architecture (target topology).
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Partition array.
 * @param perform          Percentage of moves. 
 * @param tolerance        Tolerance for cost degradation. 
 * @param k                Number of part. 
 * @param epsilon          Tolerance for balancing cost. 
 * @param algo_cluster     Clustering algorithm.
 *
 * @return Integer.
 */
INT 
multilevel(Hypergraph  * h,
           Arch        * a,
           List       ** out_neighbors, 
           List       ** in_neighbors, 
           INT         * sort, 
           PART * partition, 
           INT           perform, 
           INT           tolerance,
           INT           k,
           INT           epsilon,
           INT           algo_cluster,
           int (*intialPart)(Hypergraph * h, List ** out_neighbors, List ** in_neighbors, PART * partition, INT k, INT epsilon),
           int (*refinement)(Hypergraph * h, Arch * a, List ** out_neighbors, List ** in_neighbors, INT * sort, PART * partition, INT perform, INT tolerance, INT k)) 
{
    INT nv     = h->i_vertices;    /* number of vertices                        */
    
    INT relax  = 999;
    INT levels = log(nv);   
    
    if(algo_cluster==10) {
        levels /= 2;
    } 
    
    Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*) * levels);
    MEM_ERROR(hypergraphs);
    
    INT ** maps = (INT**)malloc(sizeof(INT*) * (levels - 1));
    MEM_ERROR(maps);
    
    INT * map_final = (INT*)malloc(sizeof(INT)*nv);
    MEM_ERROR(map_final);
    
    PART * partitionp = (PART*)malloc(sizeof(PART)*nv);
    MEM_ERROR(partitionp);
    
    INT * lambda = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
    MEM_ERROR(lambda);
    
    
    
    hypergraphs[0] = h; 
        
    for(INT i = 0; i < nv; i++) 
      {
        map_final[i] = i;
      }
 
    for(INT i = 1; i < levels; i++) 
      {  
        hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
        MEM_ERROR(hypergraphs[i]);
      
        hypergraphs[i]->is_red = NULL; hypergraphs[i]->vth = NULL;
      
        maps[i - 1]   = (INT*)malloc(sizeof(INT) * hypergraphs[i - 1]->i_vertices);
        MEM_ERROR(maps[i - 1]);
        
        if(algo_cluster == 0 || algo_cluster == 10) 
          {
            heavy_edge_matching(hypergraphs[i - 1], hypergraphs[i], out_neighbors, in_neighbors, maps[i - 1], k, epsilon);
        rbh_build_precomputed(hypergraphs[i]);
          }
        
        if(algo_cluster == 11 || algo_cluster == 1) 
          {
            best_phi_clustering(hypergraphs[i - 1], hypergraphs[i], a, out_neighbors, in_neighbors, sort, maps[i - 1], levels, epsilon, relax);
            rbh_build_precomputed(hypergraphs[i]);
          }

        for(INT u = 0; u < h->i_vertices; u++) 
          {  
            map_final[u] = maps[i - 1][map_final[u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i], out_neighbors);

        compute_list_in_neighbors_unalloc(hypergraphs[i], out_neighbors, in_neighbors);      
    }
    
    intialPart(hypergraphs[levels - 1], out_neighbors, in_neighbors, partition, k, epsilon);
    
    topological_sort(hypergraphs[levels - 1], out_neighbors, in_neighbors, sort);
    
    refinement(hypergraphs[levels - 1], a, out_neighbors, in_neighbors, sort, partition, perform, tolerance, k);

    for(INT i = levels - 1; i > 0; i--) 
      {  
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {  
            partitionp[u] = partition[maps[i - 1][u]];
          }      
        
        compute_list_neighbors_unalloc(hypergraphs[i - 1], out_neighbors);
         
        compute_list_in_neighbors_unalloc(hypergraphs[i - 1], out_neighbors, in_neighbors);
        
        topological_sort(hypergraphs[i - 1], out_neighbors, in_neighbors, sort);
        
        refinement(hypergraphs[i - 1], a, out_neighbors, in_neighbors, sort, partitionp, perform, tolerance, k);        
        
        for(INT u = 0; u < hypergraphs[i - 1]->i_vertices; u++) 
          {
            partition[u] = partitionp[u];
          }  
      }   

    /* free section */
    for(INT level = 1; level < levels; level++) 
      {   
        rbh_free(hypergraphs[level]);
        free(hypergraphs[level]);            
        free(maps[level - 1]);
      }
    free(lambda);
    free(partitionp);
    free(hypergraphs);
    free(maps);
    free(map_final);

    return 0;
}
