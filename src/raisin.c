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
/**   NAME       : raisin.c                                **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are the functions           **/
/**                definition of the application.          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/



#include "../include/rbh.h"
#include "../include/a.h"
#include "../include/crbh.h"
#include "../include/pqueue.h"
#include "../include/list.h"
#include "../include/ipart.h"
#include "../include/fm.h"
#include "../include/dlist.h"
#include "../include/m.h"
#include "../include/capacity.h"

/* Standard library import */
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <limits.h>

static char *
default_partition_path(const char *graph_path)
{
    static const char suffix[] = ".part";
    size_t length = strlen(graph_path) + sizeof(suffix);
    char *path = (char *)malloc(length);

    MEM_ERROR(path);
    snprintf(path, length, "%s%s", graph_path, suffix);
    return path;
}

static INT
parse_int_arg(const char *value, const char *name, long minimum, long maximum)
{
    char *end = NULL;
    long parsed;

    errno = 0;
    parsed = strtol(value, &end, 10);
    if (errno == ERANGE || end == value || *end != '\0' ||
        parsed < minimum || parsed > maximum) {
        fprintf(stderr, "Invalid %s: '%s' (expected %ld..%ld)\n",
                name, value, minimum, maximum);
        exit(EXIT_FAILURE);
    }
    return (INT)parsed;
}

static void
require_cli(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void
require_partition_range(INT vertices, const PART *partition, INT part_count)
{
    for (INT vertex = 0; vertex < vertices; vertex++) {
        if ((INT)partition[vertex] >= part_count) {
            fprintf(stderr,
                    "Partition value %d for vertex %d is outside 0..%d\n",
                    (INT)partition[vertex], vertex, part_count - 1);
            exit(EXIT_FAILURE);
        }
    }
}

static void
require_architecture_parts(const Arch *arch, INT part_count)
{
    if (part_count > arch->i_n) {
        fprintf(stderr,
                "part_number %d exceeds architecture size %d\n",
                part_count, arch->i_n);
        exit(EXIT_FAILURE);
    }
}

static void
require_communication_feasible(const Hypergraph *hypergraph,
                               const Arch *arch,
                               const PART *partition,
                               INT part_count)
{
    CommunicationUsage usage;
    int status = communication_usage_compute(&usage, hypergraph, arch,
                                             partition, part_count);

    require_cli(status == 0, "Unable to compute communication capacity usage");
    printf("communication signal hops;%d\n", usage.i_routed_signal_hops);
    printf("communication overloaded links;%d\n", usage.i_overloaded_links);
    printf("communication max overload;%d\n", usage.i_max_overload);
    printf("communication feasible;%s\n",
           communication_usage_is_feasible(&usage) ? "yes" : "no");

    if (!communication_usage_is_feasible(&usage)) {
        for (INT u = 0; u < arch->i_n; u++) {
            for (INT v = u + 1; v < arch->i_n; v++) {
                INT load;
                INT capacity;
                if (!arch_has_link(arch, u, v))
                    continue;
                load = communication_usage_link_load(&usage, u, v);
                capacity = arch->ti_capacity[u * arch->i_n + v];
                if (load > capacity) {
                    fprintf(stderr,
                            "Communication link %d-%d exceeds capacity: %d > %d\n",
                            u, v, load, capacity);
                }
            }
        }
        communication_usage_free(&usage);
        require_cli(false, "Communication capacity constraints are not satisfied");
    }
    communication_usage_free(&usage);
}

static void
require_resource_feasible(const Hypergraph *hypergraph,
                          const Arch *arch,
                          const PART *partition,
                          INT part_count)
{
    ResourceUsage usage;
    int status;

    if (arch->i_resources == 0) {
        printf("resource capacity mode;legacy-balance\n");
        return;
    }
    if (arch->i_resources != hypergraph->i_weights) {
        fprintf(stderr,
                "Architecture has %d resource dimensions but circuit has %d\n",
                arch->i_resources, hypergraph->i_weights);
        exit(EXIT_FAILURE);
    }

    status = resource_usage_compute(&usage, hypergraph, arch,
                                    partition, part_count);
    require_cli(status == 0, "Unable to compute resource capacity usage");
    printf("resource capacity mode;explicit\n");
    printf("resource overloaded dimensions;%d\n",
           usage.i_overloaded_dimensions);
    printf("resource overloaded parts;%d\n", usage.i_overloaded_parts);
    printf("resource max overload;%d\n", usage.i_max_overload);
    printf("resource feasible;%s\n",
           resource_usage_is_feasible(&usage) ? "yes" : "no");

    if (!resource_usage_is_feasible(&usage)) {
        for (INT part = 0; part < part_count; part++) {
            for (INT resource = 0; resource < arch->i_resources; resource++) {
                INT load = resource_usage_part_load(&usage, part, resource);
                INT capacity = arch_part_capacity(arch, part, resource);
                if (load > capacity) {
                    fprintf(stderr,
                            "Part %d resource %d exceeds capacity: %d > %d\n",
                            part, resource, load, capacity);
                }
            }
        }
        resource_usage_free(&usage);
        require_cli(false, "Resource capacity constraints are not satisfied");
    }
    resource_usage_free(&usage);
}

static void
repair_resource_partition(const Hypergraph *hypergraph,
                          const Arch *arch,
                          PART *partition,
                          INT part_count)
{
    INT moves = 0;
    int status;

    if (arch->i_resources == 0)
        return;
    if (arch->i_resources != hypergraph->i_weights) {
        fprintf(stderr,
                "Architecture has %d resource dimensions but circuit has %d\n",
                arch->i_resources, hypergraph->i_weights);
        exit(EXIT_FAILURE);
    }
    status = resource_partition_repair(hypergraph, arch, partition,
                                       part_count, &moves);
    if (status != 0) {
        fprintf(stderr, "Unable to repair resource capacities (status %d)\n",
                status);
        exit(EXIT_FAILURE);
    }
    printf("resource repair moves;%d\n", moves);
}

static void
repair_communication_partition(const Hypergraph *hypergraph,
                               const Arch *arch,
                               PART *partition,
                               INT part_count)
{
    INT moves = 0;
    int status = communication_partition_repair(
        hypergraph, arch, partition, part_count, &moves);

    if (status != 0) {
        fprintf(stderr,
                "Unable to repair communication capacities (status %d)\n",
                status);
        exit(EXIT_FAILURE);
    }
    printf("communication repair moves;%d\n", moves);
}

static void
print_usage(FILE *stream)
{
    fprintf(stream,
        "Usage:\n"
        "  raisin <graph.rzn2> cluster <hem|bsc> size <n> [options]\n"
        "  raisin <graph.rzn2> part <dbfs|ddfs|ccp> part_number <1..%d> [options]\n"
        "  raisin <graph.rzn2> refine <kfm|dkfm> part_file <file.sol> "
        "part_number <1..%d> [options]\n"
        "  raisin <graph.rzn2> multilevel cluster <hem|bsc> "
        "part <dbfs|ddfs|ccp> refine <kfm|dkfm|dkfmfast> "
        "part_number <1..%d> [options]\n"
        "  raisin <graph.rzn2> eval part_number <1..%d> "
        "partfile <file.sol> [options]\n"
        "  raisin <graph.rzn2> stats [seed <n>]\n"
        "\n"
        "Common options:\n"
        "  bfactor <n>    balance factor (default: 5)\n"
        "  partfile <p>   output prefix; RaiSin appends .sol\n"
        "  archfile <p>   target architecture (default: targets/arch0.arch)\n"
        "  seed <n>       random seed (default: 1)\n"
        "\n"
        "Refinement options:\n"
        "  perform <n>    number of refinement passes (default: 10)\n"
        "  tolerance <n>  refinement tolerance (default: 0)\n",
        RAISIN_PART_MAX, RAISIN_PART_MAX, RAISIN_PART_MAX, RAISIN_PART_MAX);
}


/**
 * @brief Function launching application.  
 *
 * @param argv         Number of parameters.
 * @param argc         Array of parameters. 
 *
 * @return Void.
 */
int main(int argv, char ** argc){

    INT seed = 1;

    if (argv < 3 || strcmp(argc[1], "help") == 0 ||
        strcmp(argc[1], "--help") == 0) {
        bool explicit_help = argv >= 2 &&
                             (strcmp(argc[1], "help") == 0 ||
                              strcmp(argc[1], "--help") == 0);
        print_usage(explicit_help ? stdout : stderr);
        return explicit_help ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    for (int i = 3; i < argv; i++) {
        if (strcmp(argc[i], "seed") == 0) {
            require_cli(i + 1 < argv, "Missing value after 'seed'");
            seed = parse_int_arg(argc[i + 1], "seed", 0, INT_MAX);
        }
    }
    srand((unsigned int)seed);
    printf("seed;%d\n", seed);

   if(strcmp(argc[2], "cluster") == 0) 
     {  
       if(argv < 6)
         {  
           print_usage(stderr);
           
           return EXIT_FAILURE;
         }
       
       printf("mode cluster\n");
       
       /* default values */
       INT epsilon = 5;
       INT k       = 0;
       INT size    = 0;
       
       char * arch_path  = NULL;
       char * graph_path = NULL;
       char * part_path  = NULL;
       char * algo       = NULL;
       
       for(int i = 6; i < argv; i++) 
         {
           if(strcmp(argc[i], "bfactor") == 0 &&  i + 1 < argv) 
             {
               epsilon = parse_int_arg(argc[i + 1], "bfactor", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "partfile") == 0 &&  i + 1 < argv) 
             {
               part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(part_path);
               strcpy(part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "archfile") == 0 && i + 1 < argv)
             {
               arch_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(arch_path);
               strcpy(arch_path, argc[i + 1]);
             }       
       }
       
       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);
       strcpy(graph_path, argc[1]);
       
       algo       = (char*)malloc(sizeof(char) * (strlen(argc[3]) + 1));
       MEM_ERROR(algo);
       strcpy(algo, argc[3]);

       require_cli(strcmp(algo, "hem") == 0 || strcmp(algo, "bsc") == 0,
                   "Unknown clustering algorithm (expected hem or bsc)");
       require_cli(strcmp(argc[4], "size") == 0,
                   "Missing required 'size' option");
       size = parse_int_arg(argc[5], "size", 1, INT_MAX);
       
       if(arch_path == NULL)
         {
          /* load default architecture */
          arch_path = (char*)malloc(sizeof(char) * 25);
          MEM_ERROR(arch_path);
          strcpy(arch_path, "targets/arch0.arch");
         }
       
       if(part_path == NULL)
         {
           /* set default partition file */
           part_path = default_partition_path(graph_path);
         }
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);

       Hypergraph * h2 = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h2);

       h2->is_red = NULL; h2->vth = NULL;

       Arch * a = (Arch*)malloc(sizeof(Arch));
       MEM_ERROR(a);       

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       require_cli(arch_load(a, arch_path, false) == 0,
                   "Unable to load target architecture");

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       a->s_arch_name = (char*)malloc(sizeof(char) * (strlen(arch_path) + 1));
       MEM_ERROR(a->s_arch_name);

       strcpy(h->s_rbh_name,  graph_path);
       strcpy(a->s_arch_name, arch_path);
       
       /* compute corresponding part number */
       INT * total_weight = (INT*)calloc(h->i_weights, sizeof(INT));
       MEM_ERROR(total_weight);
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {  
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
                total_weight[wi] += h->ti_weights[i * h->i_weights + wi];
             }
         }
       
       INT max_weights = 0;
       
       for(INT wi = 0; wi < h->i_weights; wi++)
         {
           if(max_weights < total_weight[wi])
             {
               max_weights = total_weight[wi];
             }
         }
       
       k = ceil((float)max_weights / (float)size);      
       
       /* compute neighbors and in neighbors */
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
       
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT) * h->i_vertices);
       MEM_ERROR(sort);
       
       topological_sort(h, neighbors_list, in_neighbors_list, sort);
       
       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT relax = 999; // relaxation for BSC D cost [constant]
       
       INT * map = (INT*)malloc(sizeof(INT)*h->i_vertices);
       MEM_ERROR(map);
       for(INT i = 0; i < h->i_vertices; i++)
         {
           map[i] = i;
         }       
       
       if(strcmp(algo, "hem") == 0)
         {
           INT levels = log(size);
           
           if(levels > 1){
           
               Hypergraph ** hypergraphs = (Hypergraph**)malloc(sizeof(Hypergraph*) * levels);
               MEM_ERROR(hypergraphs);
           
               INT ** maps = (INT**)malloc(sizeof(INT*) * (levels - 1));
               MEM_ERROR(maps);
           
               hypergraphs[0] = h;            
               
               for(INT i = 1; i < levels; i++) 
                 {
                        
                   hypergraphs[i] = (Hypergraph*)malloc(sizeof(Hypergraph));
                   MEM_ERROR(hypergraphs[i]);
                   hypergraphs[i]->is_red = NULL; hypergraphs[i]->vth = NULL;       
      
                   maps[i - 1]   = (INT*)malloc(sizeof(INT) * hypergraphs[i - 1]->i_vertices);
                   MEM_ERROR(maps[i - 1]);
                  
                   heavy_edge_matching(hypergraphs[i - 1], hypergraphs[i], neighbors_list, in_neighbors_list, maps[i - 1], k, epsilon);

                   for(INT u = 0; u < h->i_vertices; u++)
                     {
                       map[u] = maps[i - 1][map[u]];
                     }

                   compute_list_neighbors_unalloc(hypergraphs[i], neighbors_list);
                   compute_list_in_neighbors_unalloc(hypergraphs[i], neighbors_list, in_neighbors_list);

                   /* Build is_red and VtH after the contracted graph is fully populated */
                   rbh_build_precomputed(hypergraphs[i]);
        
                 }
           
               for(INT level = 1; level < levels; level++) 
                 {
        
                   rbh_free(hypergraphs[level]);
                   
                   free(hypergraphs[level]);            
                   
                   free(maps[level - 1]);
                 }
            }else{
                heavy_edge_matching(h, h2, neighbors_list, in_neighbors_list, map, k, epsilon);
            rbh_build_precomputed(h2);
            }
         }
       
        if(strcmp(algo, "bsc") == 0) 
          {
            best_phi_clustering(h, h2, a, neighbors_list, in_neighbors_list, sort, map, size, epsilon, relax);
          }
       
       INT pmax = compute_clustering_criticality(h, neighbors_list, in_neighbors_list, sort, map, 800);
       
       printf("%s cost;%d\n", algo, pmax);
       
       require_cli(write_int_partition(h->i_vertices, map, part_path) == 0,
                   "Unable to write clustering partition");
       
       bool * is_cluster = (bool*)malloc(h->i_vertices * sizeof(bool));
       MEM_ERROR(is_cluster);
       
       INT cluster_number = 0;
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           is_cluster[i] = false;
         }
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           is_cluster[map[i]] = true;
         }
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           if(is_cluster[map[i]]) 
             {
               cluster_number++;
               is_cluster[map[i]] = false;
             }
         }
       
       printf("%s clusters;%d\n", algo, cluster_number);
              
       /* free section */
       free(sort);
       free(map);
       free(total_weight);
       
       for(INT i = 0; i < h->i_vertices; i++)
         {

          delete_list(neighbors_list[i]);
          delete_list(in_neighbors_list[i]);        
         }

       free(neighbors_list);
       free(in_neighbors_list);
       free(is_cluster);
       rbh_free(h2);
       rbh_free(h);
       arch_free(a);
       free(h2);
       free(h);
       free(algo);
       free(arch_path);
       free(graph_path);
       free(part_path);
       
       return EXIT_SUCCESS;
      }
    
    if(strcmp(argc[2], "part") == 0) 
      {
       
       if(argv < 6) 
         {
           print_usage(stderr);
           return EXIT_FAILURE;
         }
       
       /* default values */
       INT epsilon = 5;
       INT k       = 0;
       
       char * arch_path  = NULL;
       char * graph_path = NULL;
       char * part_path  = NULL;
       char * algo       = NULL;
       
       for(int i = 2; i < argv; i++) 
         {
           if(strcmp(argc[i], "part_number") == 0 && i + 1 < argv) 
             {
               k = parse_int_arg(argc[i + 1], "part_number", 1, RAISIN_PART_MAX);
             }
       
           if(strcmp(argc[i], "bfactor") == 0 && i + 1 < argv) 
             {
               epsilon = parse_int_arg(argc[i + 1], "bfactor", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "partfile") == 0 && i + 1 < argv) 
             {
               part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(part_path);
               strcpy(part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "archfile") == 0 && i + 1 < argv) 
             {
               arch_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(arch_path);
               strcpy(arch_path, argc[i + 1]);
             }
         }

       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);
       
       strcpy(graph_path, argc[1]);
       
       algo       = (char*)malloc(sizeof(char) * (strlen(argc[3]) + 1));
       MEM_ERROR(algo);
       
       strcpy(algo, argc[3]);

       require_cli(k > 0, "Missing required 'part_number' option");
       require_cli(strcmp(algo, "dbfs") == 0 || strcmp(algo, "ddfs") == 0 ||
                   strcmp(algo, "ccp") == 0,
                   "Unknown partitioning algorithm (expected dbfs, ddfs or ccp)");

       if(arch_path == NULL) 
         {
          /* load default architecture */
          arch_path = (char*)malloc(sizeof(char) * 25);
          MEM_ERROR(arch_path);

          strcpy(arch_path, "targets/arch0.arch");
         }
       
       if(part_path == NULL) 
         {
           /* set default partition file */
           part_path = default_partition_path(graph_path);
         }
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);

       Arch * a = (Arch*)malloc(sizeof(Arch));
       MEM_ERROR(a);       

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       require_cli(arch_load(a, arch_path, false) == 0,
                   "Unable to load target architecture");
       require_architecture_parts(a, k);

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       a->s_arch_name = (char*)malloc(sizeof(char) * (strlen(arch_path) + 1));
       MEM_ERROR(a->s_arch_name);

       strcpy(h->s_rbh_name,  graph_path);
       strcpy(a->s_arch_name, arch_path);
       
       /* compute neighbors and in neighbors */
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
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT) * h->i_vertices);
       MEM_ERROR(sort);

       topological_sort(h, neighbors_list, in_neighbors_list, sort);
       
       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT * lambda = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
       MEM_ERROR(lambda);
       
       PART * partition = (PART*)malloc(sizeof(PART)*h->i_vertices);
       MEM_ERROR(partition);

       for(INT i = 0; i < h->i_vertices; i++) 
         {
           partition[i] = i;
         }       
       
       if(strcmp(algo, "dbfs") == 0)
         {
           derived_breadth_first_search(h, neighbors_list, in_neighbors_list, partition, k, epsilon);
           
           INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
           
           INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
           printf("dbfs pmax;%d\n", pmax);
           printf("dbfs cut;%d\n", cut);
         }
       
       if(strcmp(algo, "ddfs") == 0) 
         {
           derived_depth_first_search(h, neighbors_list, in_neighbors_list, partition, k, epsilon);
           
           INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
           
           INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
           printf("ddfs pmax;%d\n", pmax);
           printf("ddfs cut;%d\n", cut);
         }
       
       if(strcmp(algo, "ccp") == 0)
         {
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
           
           INT * map = (INT*)malloc(sizeof(INT) * h->i_vertices);
           MEM_ERROR(map);
           
           INT D = ceil((float)crit_max * 0.8); // default value ~ 20 % of CP 

           critical_connected_component_partitioning(h, neighbors_list, in_neighbors_list, map, k, epsilon,  crit_max-D);

           bool * is_in = (bool*)malloc(sizeof(bool) * h->i_vertices);
           MEM_ERROR(is_in);

           for(INT i = 0; i < h->i_vertices; i++) 
             {
               is_in[i] = false;
             }

           for(INT i = 0; i < h->i_vertices; i++) 
             {
               is_in[map[i]] = true;
             }

           INT nvp = 0;
           
           for(INT i = 0; i < h->i_vertices; i++) 
             {
               if(is_in[map[i]]) 
                 {
                   nvp++;
                   is_in[map[i]] = false;
                 }
             }  
           
           PART * partitionp = (PART*)malloc(sizeof(PART)*nvp);
           MEM_ERROR(partitionp);         
       
           compute_subhypergraph(h, h2, map, nvp);
           
           compute_list_neighbors_unalloc(h2, neighbors_list);
           
           compute_list_in_neighbors_unalloc(h2, neighbors_list, in_neighbors_list);
           
           derived_breadth_first_search_multilevel(h2, neighbors_list, in_neighbors_list, partitionp, k, epsilon);
           
            for(INT i = 0; i < h->i_vertices; i++)
              {
                partition[i] = partitionp[map[i]];
              }

            compute_list_neighbors_unalloc(h, neighbors_list);
            compute_list_in_neighbors_unalloc(h, neighbors_list,
                                               in_neighbors_list);
            topological_sort(h, neighbors_list, in_neighbors_list, sort);

            INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition);
           
           INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
           printf("ccp pmax;%d\n", pmax);
           printf("ccp cut;%d\n", cut);
           
           rbh_free(h2);
           free(h2);
           free(map);
           free(is_in);
         }

       repair_resource_partition(h, a, partition, k);
       require_resource_feasible(h, a, partition, k);
       repair_communication_partition(h, a, partition, k);
       require_communication_feasible(h, a, partition, k);
       require_cli(write_partition(h->i_vertices, partition, part_path) == 0,
                   "Unable to write partition");
       
       INT * size_parts = (INT*)calloc(k * h->i_weights, sizeof(INT));
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               size_parts[partition[i] * h->i_weights+wi] += h->ti_weights[i * h->i_weights + wi];
             }
         }
       
       double avg = 0;
       double var = 0;
       double sigma;
       
       INT max_dif = 0;
       
       for(INT p = 0; p < k; p++)  
         {
           printf("size of part %d", p);
           
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               printf(";%d", size_parts[p * h->i_weights + wi]);
               
               avg += size_parts[p * h->i_weights + wi];
               
               var += size_parts[p * h->i_weights + wi] * size_parts[p * h->i_weights + wi];
             }
           
           printf("\n");
         }
       
       avg /= (k * h->i_weights);
       var /= (k * h->i_weights);
       
       sigma = sqrt(var - avg*avg);
       
       for(INT p = 0; p < k; p++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               if(max_dif < size_parts[p * h->i_weights + wi] - avg) 
                 {
                   max_dif = size_parts[p * h->i_weights + wi] - avg;
                 }
             }
         }

       printf("weight std deviation;%.2f\n", sigma);

       printf("balance;%.2f\n", (float)max_dif * 100.0 / (float)h->i_vertices);

       /* free section */
       free(size_parts);
       free(sort);
       free(partition);

       for(INT i = 0; i < h->i_vertices; i++)
         {

          delete_list(neighbors_list[i]);
          
          delete_list(in_neighbors_list[i]);        

         }
       
       free(neighbors_list);
       free(in_neighbors_list);
       free(algo);
       rbh_free(h);
       arch_free(a);
       free(h);
       free(lambda);
       free(arch_path);
       free(graph_path);
       free(part_path);
       
       return EXIT_SUCCESS;
       
    }
    
    if(strcmp(argc[2], "refine") == 0) 
      {
        if(argv < 8) 
          {
            print_usage(stderr);
            
            return EXIT_FAILURE;
          }
       
       printf("mode refine\n");
       
       /* default values */
       INT epsilon   = 5;
       INT k         = 0;
       INT perform   = 10;
       INT tolerance = 0;
       
       char * arch_path      = NULL;
       char * graph_path     = NULL;
       char * part_path      = NULL;
       char * algo           = NULL;
       char * init_part_path = NULL;
       
       for(int i = 2; i < argv; i++) 
         {
           if(strcmp(argc[i], "part_number") == 0 && i + 1 < argv)
             {
               k = parse_int_arg(argc[i + 1], "part_number", 1, RAISIN_PART_MAX);
             }
       
           if(strcmp(argc[i], "bfactor") == 0 && i + 1 < argv) 
             {
               epsilon = parse_int_arg(argc[i + 1], "bfactor", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "perform") == 0 && i + 1 < argv) 
             {
               perform = parse_int_arg(argc[i + 1], "perform", 1, INT_MAX);
             }
           
           if(strcmp(argc[i], "tolerance") == 0 && i + 1 < argv) 
             {
               tolerance = parse_int_arg(argc[i + 1], "tolerance", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "part_file") == 0 && i + 1 < argv) 
             {
               init_part_path = (char*)malloc(sizeof(char) * (strlen(argc[i+1]) + 1));
               MEM_ERROR(init_part_path);

               strcpy(init_part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "partfile") == 0 && i + 1 < argv) 
             {
               part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(part_path);
               strcpy(part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "archfile") == 0 && i + 1 < argv) 
             {
               arch_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(arch_path);

               strcpy(arch_path, argc[i + 1]);
             }
           }

       RAISIN_UNUSED(epsilon);

       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);
       strcpy(graph_path, argc[1]);
       
       algo       = (char*)malloc(sizeof(char) * (strlen(argc[3]) + 1));
       MEM_ERROR(algo);
       strcpy(algo, argc[3]);

       require_cli(k > 0, "Missing required 'part_number' option");
       require_cli(init_part_path != NULL, "Missing required 'part_file' option");
       require_cli(strcmp(algo, "kfm") == 0 || strcmp(algo, "dkfm") == 0,
                   "Unknown refinement algorithm (expected kfm or dkfm)");

       if(arch_path == NULL) 
         {
          /* load default architecture */
          arch_path = (char*)malloc(sizeof(char) * 25);
          MEM_ERROR(arch_path);
          
          strcpy(arch_path, "targets/arch0.arch");
         }
       
       if(part_path == NULL) 
         {
           /* set default partition file */
           part_path = default_partition_path(graph_path);
         }
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);

       Arch * a = (Arch*)malloc(sizeof(Arch));
       MEM_ERROR(a);       

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       require_cli(arch_load(a, arch_path, false) == 0,
                   "Unable to load target architecture");
       require_architecture_parts(a, k);

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       a->s_arch_name = (char*)malloc(sizeof(char) * (strlen(arch_path) + 1));
       MEM_ERROR(a->s_arch_name);

       strcpy(h->s_rbh_name,  graph_path);
       strcpy(a->s_arch_name, arch_path);
            
       /* compute neighbors and in neighbors */
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
       
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT)*h->i_vertices);
       MEM_ERROR(sort);
       
       topological_sort(h, neighbors_list, in_neighbors_list, sort);
       
       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT * lambda = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
       
       PART * partition = (PART*)malloc(sizeof(PART)*h->i_vertices);
       MEM_ERROR(partition);
       
       require_cli(load_partition(h->i_vertices, partition, init_part_path) == 0,
                   "Unable to load initial partition");
       require_partition_range(h->i_vertices, partition, k);

       if(strcmp(algo, "kfm") == 0) 
         {
       
           kfm(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k);
       
           INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
           
           INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
           printf("\nkfm pmax;%d\n", pmax);
           printf("\nkfm cut;%d\n", cut);
       
         }
       
       if(strcmp(algo, "dkfm") == 0) 
         {
           dkfm(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k);
           
           INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
           
           INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
           printf("\ndkfm pmax;%d\n", pmax);
           printf("\ndkfm cut;%d\n", cut);
         }
       
       repair_resource_partition(h, a, partition, k);
       require_resource_feasible(h, a, partition, k);
       repair_communication_partition(h, a, partition, k);
       require_communication_feasible(h, a, partition, k);
       require_cli(write_partition(h->i_vertices, partition, part_path) == 0,
                   "Unable to write refined partition");
       
       INT * size_parts = (INT*)calloc(k*h->i_weights, sizeof(INT));
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               size_parts[partition[i] * h->i_weights + wi] += h->ti_weights[i * h->i_weights + wi];          
             }
         }
       
       double avg = 0;
       double var = 0;
       double sigma;
       
       INT max_dif = 0;
       
       for(INT p = 0; p < k; p++) 
         {
           printf("size of part %d", p);
           
           for(INT wi = 0; wi < h->i_weights; wi++)  
             {
               printf(";%d", size_parts[p * h->i_weights + wi]);
               
               avg += size_parts[p * h->i_weights + wi];
               
               var += size_parts[p * h->i_weights + wi] * size_parts[p * h->i_weights + wi];
             }
           
           printf("\n");
         }
       
       avg /= (k * h->i_weights);
       var /= (k * h->i_weights);
       
       sigma = sqrt(var - avg * avg);
       
       for(INT p = 0; p < k; p++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               if(max_dif < size_parts[p * h->i_weights + wi] - avg) 
                 {
                   max_dif = size_parts[p * h->i_weights + wi] - avg;
                 }
             }
         }

       printf("weight std deviation;%.2f\n", sigma);

       printf("balance;%.2f\n", (float)max_dif * 100.0 / (float)h->i_vertices);
             
       /* free section */
       free(size_parts);
       free(sort);
       free(partition);
       
       for(INT i = 0; i < h->i_vertices; i++)
         {
           delete_list(neighbors_list[i]);
           delete_list(in_neighbors_list[i]);        
         }
       
       free(neighbors_list);
       free(in_neighbors_list);
       free(lambda);
       rbh_free(h);
       arch_free(a);
       free(h);
       free(init_part_path);
       free(algo);
       free(arch_path);
       free(graph_path);
       free(part_path);
       return EXIT_SUCCESS;
       }
    
    if(strcmp(argc[2], "multilevel") == 0) 
      { 
        if(argv < 11) 
          {
            print_usage(stderr);
            return EXIT_FAILURE;
          }
       
       printf("mode multilevel\n");
       
       /* default values */
       INT epsilon   = 5;
       INT k         = 0;
       INT perform   = 10;
       INT tolerance = 0;
       INT algo_mode = 0;
       
       char * arch_path      = NULL;
       char * graph_path     = NULL;
       char * part_path      = NULL;
       char * algo_cluster   = NULL;
       char * algo_part      = NULL;
       char * algo_refine    = NULL;
       char * init_part_path = NULL;
       
       for(int i = 2; i < argv; i++) 
         {
           if(strcmp(argc[i], "cluster") == 0 && i + 1 < argv) 
             {
               algo_cluster = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(algo_cluster);

               strcpy(algo_cluster, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "part") == 0 && i + 1 < argv) 
             {
               algo_part = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(algo_part);
               
               strcpy(algo_part, argc[i+1]);
             }
           
           if(strcmp(argc[i], "refine") == 0 && i + 1 < argv) 
             {
               algo_refine = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(algo_refine);

               strcpy(algo_refine, argc[i + 1]);
             }

           if(strcmp(argc[i], "part_number") == 0 && i + 1 < argv) 
             {
               k = parse_int_arg(argc[i + 1], "part_number", 1, RAISIN_PART_MAX);
             }
       
           if(strcmp(argc[i], "bfactor") == 0 && i + 1 < argv) 
             {
               epsilon = parse_int_arg(argc[i + 1], "bfactor", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "perform") == 0 && i + 1 < argv) 
             {
               perform = parse_int_arg(argc[i + 1], "perform", 1, INT_MAX);
             }
           
           if(strcmp(argc[i], "tolerance") == 0 && i + 1 < argv) 
             {
               tolerance = parse_int_arg(argc[i + 1], "tolerance", 0, INT_MAX);
             }
           
           if(strcmp(argc[i], "part_file") == 0 && i + 1 < argv) 
             {
               init_part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(init_part_path);
               
               strcpy(init_part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "partfile") == 0 && i + 1 < argv) 
             {
               part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(part_path);
               
               strcpy(part_path, argc[i + 1]);
             }
           
           if(strcmp(argc[i], "archfile") == 0 && i + 1 < argv) 
             {
               arch_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(arch_path);
               
               strcpy(arch_path, argc[i + 1]);
             }
         }
       
         if(algo_refine == NULL || algo_cluster == NULL || algo_part == NULL) 
           {
              print_usage(stderr);
             return EXIT_FAILURE;
           }

       require_cli(k > 0, "Missing required 'part_number' option");
       require_cli(strcmp(algo_cluster, "hem") == 0 || strcmp(algo_cluster, "bsc") == 0,
                   "Unknown clustering algorithm (expected hem or bsc)");
       require_cli(strcmp(algo_part, "dbfs") == 0 || strcmp(algo_part, "ddfs") == 0 ||
                   strcmp(algo_part, "ccp") == 0,
                   "Unknown partitioning algorithm (expected dbfs, ddfs or ccp)");
       require_cli(strcmp(algo_refine, "kfm") == 0 || strcmp(algo_refine, "dkfm") == 0 ||
                   strcmp(algo_refine, "dkfmfast") == 0,
                   "Unknown refinement algorithm (expected kfm, dkfm or dkfmfast)");

       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);
       
       strcpy(graph_path, argc[1]);
       
       if(arch_path == NULL)
         {
           /* load default architecture */
           arch_path = (char*)malloc(sizeof(char) * 25);
           MEM_ERROR(arch_path);
          
           strcpy(arch_path, "targets/arch0.arch");
         }
       
       if(part_path == NULL) 
         {
           /* set default partition file */
           part_path = default_partition_path(graph_path);
         }
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);

       Arch * a = (Arch*)malloc(sizeof(Arch));
       MEM_ERROR(a);       

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       require_cli(arch_load(a, arch_path, false) == 0,
                   "Unable to load target architecture");
       require_architecture_parts(a, k);

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       a->s_arch_name = (char*)malloc(sizeof(char) * (strlen(arch_path) + 1));
       MEM_ERROR(a->s_arch_name);

       strcpy(h->s_rbh_name,  graph_path);

       strcpy(a->s_arch_name, arch_path);
            
       /* compute neighbors and in neighbors */
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
       
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT)*h->i_vertices);
       MEM_ERROR(sort);
       
       topological_sort(h, neighbors_list, in_neighbors_list, sort);
       
       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT * lambda = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
       
       PART * partition = (PART*)malloc(sizeof(PART)*h->i_vertices);
       MEM_ERROR(partition);
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
            algo_mode = 0;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, kfm);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
            algo_mode = 0;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, kfm);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
             algo_mode = 0;
             multilevel_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
            algo_mode = 1;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, kfm);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
            algo_mode = 1;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, kfm);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "kfm") == 0) 
         {    
            algo_mode = 1;
            multilevel_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 0;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, dkfm);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 0;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, dkfm);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 0;
            multilevel_pmax(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 1;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, dkfm);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 1;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, dkfm);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "dkfm") == 0) 
         {    
            algo_mode = 1;
            multilevel_pmax(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 10;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, dkfm_fast);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 10;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, dkfm_fast);
         }
       
       if(strcmp(algo_cluster, "hem") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 10;
            multilevel_pmax(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "dbfs") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 11;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_breadth_first_search_multilevel, dkfm_fast);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ddfs") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 11;
            multilevel(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode, derived_depth_first_searchMultilevel, dkfm_fast);
         }
       
       if(strcmp(algo_cluster, "bsc") == 0 && strcmp(algo_part, "ccp") == 0 && strcmp(algo_refine, "dkfmfast") == 0) 
         {    
            algo_mode = 11;
            multilevel_pmax(h, a, neighbors_list, in_neighbors_list, sort, partition, perform, tolerance, k, epsilon, algo_mode);
         }
       
       INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
       
       INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
       printf("\n%s+%s+%s pmax;%d\n", algo_cluster, algo_part, algo_refine, pmax);
       
       printf("\n%s+%s+%s cut;%d\n", algo_cluster, algo_part, algo_refine, cut);
       
       INT * size_parts = (INT*)calloc(k * h->i_weights, sizeof(INT));
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               size_parts[partition[i] * h->i_weights + wi] += h->ti_weights[i * h->i_weights + wi];         
             }
         }
       
       double avg = 0;
       double var = 0;
       double sigma;
       
       INT max_dif = 0;
       
       for(INT p = 0; p < k; p++) 
         {
           printf("size of part %d", p);
           
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               printf(";%d", size_parts[p * h->i_weights + wi]);
               
               avg += size_parts[p * h->i_weights + wi];

               var += size_parts[p * h->i_weights + wi] * size_parts[p * h->i_weights + wi];
             }
           printf("\n");
         }
       
       avg /= (k * h->i_weights);
       var /= (k * h->i_weights);
       
       sigma = sqrt(var - avg * avg);
       
       for(INT p = 0; p < k; p++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               if(max_dif < size_parts[p * h->i_weights + wi] - avg) 
                 {
                   max_dif = size_parts[p * h->i_weights + wi] - avg;
                 }
             }
         }

       printf("weight std deviation;%.2f\n", sigma);
       
       printf("balance;%.2f\n", (float)max_dif * 100.0 / (float)h->i_vertices);
       
       repair_resource_partition(h, a, partition, k);
       require_resource_feasible(h, a, partition, k);
       repair_communication_partition(h, a, partition, k);
       require_communication_feasible(h, a, partition, k);
       require_cli(write_partition(h->i_vertices, partition, part_path) == 0,
                   "Unable to write multilevel partition");
       
       /* free section */
       free(size_parts);
       free(sort);
       free(partition);
       
       for(INT i = 0; i < h->i_vertices; i++)
         {
          delete_list(neighbors_list[i]);
          delete_list(in_neighbors_list[i]);        
         }

       free(neighbors_list);
       free(in_neighbors_list);
       free(lambda);
       rbh_free(h);
       arch_free(a);
       free(h);
       free(algo_cluster);
       free(algo_part); 
       free(algo_refine); 
       free(arch_path);
       free(graph_path);
       free(part_path);
       return EXIT_SUCCESS;
       }

    if(strcmp(argc[2], "eval") == 0) 
      { 
        if(argv < 7) 
          { 
             print_usage(stderr);
            return EXIT_FAILURE;
          }
       
       printf("mode eval\n");
       
       /* default values */
       INT k = 0;
       
       char * arch_path      = NULL;
       char * graph_path     = NULL;
       char * part_path      = NULL;
       char * algo           = NULL;
       char * init_part_path = NULL;
       
       for(int i = 2; i < argv; i++) 
         {
           if(strcmp(argc[i], "part_number") == 0 && i + 1 < argv) 
             {
               k = parse_int_arg(argc[i + 1], "part_number", 1, RAISIN_PART_MAX);
              }
       
           if(strcmp(argc[i], "partfile") == 0 && i + 1 < argv) 
             {
               init_part_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(init_part_path);
               
               strcpy(init_part_path, argc[i + 1]);
             }

           if(strcmp(argc[i], "archfile") == 0 && i + 1 < argv) 
             {
               arch_path = (char*)malloc(sizeof(char) * (strlen(argc[i + 1]) + 1));
               MEM_ERROR(arch_path);

               strcpy(arch_path, argc[i + 1]);
             }
           }

       require_cli(k > 0, "Missing required 'part_number' option");
       require_cli(init_part_path != NULL, "Missing required 'partfile' option");

       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);

       strcpy(graph_path, argc[1]);
       
       algo       = (char*)malloc(sizeof(char) * (strlen(argc[3]) + 1));
       MEM_ERROR(algo);

       strcpy(algo, argc[3]);
       
       if(arch_path == NULL) 
         {
          /* load default architecture */
          arch_path = (char*)malloc(sizeof(char) * 25);
          MEM_ERROR(arch_path);

          strcpy(arch_path, "targets/arch0.arch");
         }
       
       if(part_path == NULL) 
         {
           /* set default partition file */
           part_path = default_partition_path(graph_path);
         }
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);

       Arch * a = (Arch*)malloc(sizeof(Arch));
       MEM_ERROR(a);       

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       require_cli(arch_load(a, arch_path, false) == 0,
                   "Unable to load target architecture");
       require_architecture_parts(a, k);

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       a->s_arch_name = (char*)malloc(sizeof(char) * (strlen(arch_path) + 1));
       MEM_ERROR(a->s_arch_name);

       strcpy(h->s_rbh_name,  graph_path);
       strcpy(a->s_arch_name, arch_path);
            
       /* compute neighbors and in neighbors */
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
       
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT) * h->i_vertices);
       MEM_ERROR(sort);

       topological_sort(h, neighbors_list, in_neighbors_list, sort);
       
       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT * lambda = (INT*)malloc(sizeof(INT) * h->i_hyperedges);
       
       PART * partition = (PART*)malloc(sizeof(PART)*h->i_vertices);
       MEM_ERROR(partition);
       
       require_cli(load_partition(h->i_vertices, partition, init_part_path) == 0,
                   "Unable to load partition for evaluation");
       require_partition_range(h->i_vertices, partition, k);

       INT pmax = compute_partition_criticality(h, a, neighbors_list, in_neighbors_list, sort, partition); 
       
       INT cut = compute_partition_cut(h, a, neighbors_list, in_neighbors_list, sort, partition, lambda, k);
           
       printf("\neval pmax;%d\n", pmax);
       printf("\neval cut;%d\n", cut);
       require_resource_feasible(h, a, partition, k);
       require_communication_feasible(h, a, partition, k);
       
       INT * size_parts = (INT*)calloc(k * h->i_weights, sizeof(INT));
       
       for(INT i = 0; i < h->i_vertices; i++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               size_parts[partition[i] * h->i_weights + wi] += h->ti_weights[i * h->i_weights + wi];
             }
         }
       
       double avg = 0;
       double var = 0;
       double sigma;
       
       INT max_dif = 0;
       
       for(INT p = 0; p < k; p++) 
         {
           printf("size of part %d", p);

           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               printf(";%d", size_parts[p * h->i_weights + wi]);

               avg += size_parts[p * h->i_weights + wi];
               var += size_parts[p*h->i_weights+wi] * size_parts[p * h->i_weights + wi];
           }

           printf("\n");
         }
       
       avg /= (k * h->i_weights);
       var /= (k * h->i_weights);
       
       sigma = sqrt(var - avg * avg);
       
       for(INT p = 0; p < k; p++) 
         {
           for(INT wi = 0; wi < h->i_weights; wi++) 
             {
               if(max_dif < size_parts[p * h->i_weights + wi] - avg) 
                 {
                   max_dif = size_parts[p * h->i_weights + wi] - avg;
                 }
             }
         }

       printf("weight std deviation;%.2f\n", sigma);
       
       printf("balance;%.2f\n", (float)max_dif * 100.0 / (float)h->i_vertices);
           
       /* free section */
       free(size_parts);
       free(sort);
       free(partition);

       for(INT i = 0; i < h->i_vertices; i++)
         {
          delete_list(neighbors_list[i]);
          
          delete_list(in_neighbors_list[i]);        
         }

       free(neighbors_list);
       free(in_neighbors_list);
       free(lambda);
       rbh_free(h);
       arch_free(a);
       free(h);
       free(init_part_path);
       free(algo);
       free(arch_path);
       free(graph_path);
       free(part_path);
       return EXIT_SUCCESS;
       }
    
    if(strcmp(argc[2], "stats") == 0) 
      { 
        if(argv < 3) 
          { 
             print_usage(stderr);
            return EXIT_FAILURE;
          }
       
       printf("mode stats\n");
       
       /* default values */
       char * graph_path     = NULL;
       
       graph_path = (char*)malloc(sizeof(char) * (strlen(argc[1]) + 1));
       MEM_ERROR(graph_path);

       strcpy(graph_path, argc[1]);
       
       Hypergraph * h = (Hypergraph*)malloc(sizeof(Hypergraph));
       MEM_ERROR(h);      

       require_cli(rbh_load(h, graph_path, 0, false) == 0,
                   "Unable to load hypergraph");
       require_cli(rbh_validate(h) == 0, "Invalid hypergraph structure");

       h->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(graph_path) + 1));
       MEM_ERROR(h->s_rbh_name);

       strcpy(h->s_rbh_name,  graph_path);
       
       /* compute neighbors and in neighbors */
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
       
       compute_list_in_neighbors(h, neighbors_list, in_neighbors_list);
       
       /* Compute topological sort */
       INT * sort = (INT*)malloc(sizeof(INT) * h->i_vertices);
       MEM_ERROR(sort);
       topological_sort(h, neighbors_list, in_neighbors_list, sort);

       compute_criticality(h, neighbors_list, in_neighbors_list, sort);
       
       INT * lpmax = (INT *) malloc(sizeof(INT) * (h->i_vertices + 1));
       MEM_ERROR(lpmax);
       
       INT * depth = (INT *) calloc(h->i_vertices + 1, sizeof(INT));
       MEM_ERROR(depth);
       
       INT * lmax = (INT *) malloc(sizeof(INT) * (h->i_vertices + 1));
       MEM_ERROR(lmax);
       INT * lmaxdeg = (INT *) malloc(sizeof(INT) * (h->i_vertices + 1));
       MEM_ERROR(lmaxdeg);
       
       INT * lmaxcon = (INT *) malloc(sizeof(INT) * (h->i_vertices + 1));
       MEM_ERROR(lmaxcon);
       
       INT  sizelmax    = 0;
       INT  sizelpmax  = 0;
       INT  sizelmaxdeg = 0;
       INT  sizelmaxcon = 0;
       
       float avg_deg = 0;
       
       float stdw_deg = 0;
       float avg_con  = 0;
       float stdw_con = 0;
       
       float avg_pmax  = 0;
       float stdw_pmax = 0;
       
       float avg_length  = 0;
       float stdw_length = 0;
       
       printf("#vertices;%d\n", h->i_vertices);
       printf("#hyperedges;%d\n", h->i_hyperedges);
       printf("#reds;%d\n", h->i_reds);
       
       /* critical path value */
       INT pmax = compute_pmax(h, neighbors_list, in_neighbors_list, sort, &sizelpmax, lpmax, &avg_pmax, &stdw_pmax);
       INT max_deg = compute_maxdeg(h, neighbors_list, in_neighbors_list, sort, &sizelmaxdeg, lmaxdeg, &avg_deg, &stdw_deg);
       INT max_con =
       compute_maxcon(h, neighbors_list, in_neighbors_list, sort, &sizelmaxcon, lmaxcon, &avg_con, &stdw_con);
       
       INT max_length = compute_path_length(h, neighbors_list, in_neighbors_list, sort, &sizelmax, lmax, depth, &avg_length, &stdw_length);
       
       printf("maxdeg;%d\n", max_deg);
       
       printf("#maxdeg vertex;%d\n", sizelmaxdeg);
       
       printf("avg deg;%.2f\n", avg_deg);
       
       printf("stdw deg;%.2f\n", stdw_deg);
       
       printf("maxcon;%d\n", max_con);
       printf("#maxcon hyperedge;%d\n", sizelmaxcon);
       
       printf("avg con;%0.2f\n", avg_con);
       
       printf("stdw con;%.2f\n", stdw_con);
       
       printf("pmax;%d\n", pmax);
       
       printf("#critical vertex;%d\n", sizelpmax);
       
       printf("avg pmax;%.2f\n", avg_pmax);
       
       printf("stdw pmax;%.2f\n", stdw_pmax);
       
       printf("length max;%d\n", max_length);
       
       printf("#longest paths;%d\n", sizelmax);
       
       printf("avg length;%.2f\n", avg_length);
       
       printf("stdw length;%.2f\n", stdw_length);
              
       printf("hyperedges per vertex;%.2f\n", (float)h->i_hyperedges/(float)h->i_vertices);
       
       printf("reds prop;%.2f\n", (float)h->i_reds/(float)h->i_vertices);
       
       printf("critical vertex prop;%.2f\n", (float)sizelpmax/(float)h->i_vertices);
       
       /* free section */
       
       free(sort);
       
       for(INT i = 0; i < h->i_vertices; i++)
         {
           delete_list(neighbors_list[i]);
           delete_list(in_neighbors_list[i]);         
         }

       free(lmaxcon);
       free(lmaxdeg);
       free(lpmax);
       free(lmax);
       free(depth);
       free(neighbors_list);
       free(in_neighbors_list);
       rbh_free(h);
       free(h);
       free(graph_path);
       return EXIT_SUCCESS;
     }
    
    fprintf(stderr, "Unknown mode: %s\n", argc[2]);
    return EXIT_FAILURE;
}
