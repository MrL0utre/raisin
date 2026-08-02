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
/**   NAME       : crbh.c                                  **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are definitions for the     **/
/**                clustering functions.                   **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 01 jan 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/************************************************************/

#include "crbh.h"
#include "vth.h"
#include "rng.h"

/**
 * @brief Function to swap elements.
 *
 * @param a            Pointer.
 * @param b            Pointer.
 *
 * @return void.
 */
void 
swap(int * a, 
     int * b) 
{
  int t = * a;
  * a = * b;
  * b = t;
}

/**
 * @brief Function to find the partition position for quicksort.
 *
 * @param array       Array of integers.
 * @param keys        Arrays of integers(keys)
 * @param low         Lowerbound.
 * @param high        Upperbound.
 *
 * @return Partition position.
 */
static inline void _pm3(int a[],int k[],int lo,int hi){
  int mid=lo+(hi-lo)/2;
  if(a[k[lo]]>a[k[mid]])  swap(&k[lo],&k[mid]);
  if(a[k[lo]]>a[k[hi]])   swap(&k[lo],&k[hi]);
  if(a[k[mid]]>a[k[hi]])  swap(&k[mid],&k[hi]);
  swap(&k[mid],&k[hi]);
}
int 
partition(int array[], 
          int keys[], 
          int low, 
          int high) 
{
  if(high-low>=2) _pm3(array,keys,low,high);
  int pivot=array[keys[high]];
  int i=(low-1);
  for(int j=low;j<high;j++){if(array[keys[j]]<=pivot){i++;swap(&keys[i],&keys[j]);}}
  swap(&keys[i+1],&keys[high]);
  return(i+1);
}

/**
 * @brief Sort function (quicksort) of an array of integers.
 *
 * @param array       Array of integers.
 * @param keys        Arrays of integers(keys)
 * @param low         Lowerbound.
 * @param high        Upperbound.
 *
 * @return Partition position.
 */
void 
quick_sort(int array[], 
           int keys[], 
           int low, 
           int high) 
{
  if (low < high) 
    {

      int pi = partition(array, keys, low, high);

      // recursive call on the left of pivot
      quick_sort(array, keys, low, pi - 1);

      // recursive call on the right of pivot
      quick_sort(array, keys, pi + 1, high);
    }
}

/**
 * @brief Compute the critical path length of a clustering partition
 *        in a red-black hypergraph.
 *
 * @param h                Input hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 * @param D                Delay factor.
 *
 * @return The critical path length (maximum accumulated delay).
 */
int 
compute_clustering_criticality(Hypergraph  * h, 
                               List       ** neighbors, 
                               List       ** in_neighbors, 
                               INT * sort, 
                               INT         * partition, 
                               INT           D) 
{
  RAISIN_UNUSED(in_neighbors);
  INT nv     = h->i_vertices; /* number of vertices */
  INT nr     = h->i_reds;     /* number of red vertices */
  INT * reds = h->ti_reds;    /* array of red vertices */
  /* Allocate memory for delays */
  INT  * delays = (INT *) malloc(sizeof(INT) * nv);
  MEM_ERROR(delays);

  bool * is_red   = (bool *) malloc(sizeof(bool) * nv);
  MEM_ERROR(is_red);
  bool * flag = (bool *) malloc(sizeof(bool) * nv);
  MEM_ERROR(flag);

  for (INT i = 0; i < nv; i++)
    {
      delays[i] = h->ti_delays[i];
      flag[i]   = false;
      is_red[i] = false;
    }

  for (INT i = 0; i < nr; i++) 
    {
      is_red[reds[i]] = true;
    }
  
  INT v;

  for (INT i = 0; i < nv; i++) 
    {
      v       = sort[i];
      flag[v] = true;

      if (is_red[v])
        {
          List * cell_neighbors = neighbors[v];
          while (cell_neighbors != NULL && neighbors[v]->size > 0) 
            {
              //selection of vertex v
              INT u = cell_neighbors->i;
              
              if (partition[u] != partition[v])
                {
                  delays[u] = MAX(delays[u], h->ti_delays[v] + h->ti_delays[u] + D);
                }
                
                cell_neighbors = cell_neighbors->next;
            }
         }
        else
          {
            List * cell_neighbors = neighbors[v];
            while (cell_neighbors != NULL && neighbors[v]->size > 0) 
              {
                //selection of vertex v
                INT u = cell_neighbors->i;
                if (partition[u] != partition[v])
                  {
                    delays[u] = MAX(delays[u], delays[v] + h->ti_delays[u] + D);
                  }
                
                cell_neighbors = cell_neighbors->next;
              }
          }
    }

  int max_delays = 0;

  for (int v = 0; v < nv; v++)
    {

		  if (delays[v] > max_delays)
        {
			    max_delays = delays[v];
		   	}
    }


  free(delays);
  free(flag);

  return max_delays;

}

/**
 * @brief Compute a matching (map array) of a red-black hypergraph (h1) vertices according 
 * to criticality weighting and construct a new contracted 
 * red-black hypergraph (h2).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 * @param D                Delay factor.
 *
 * @return Return a flag.
 */
int 
heavy_edge_matching(Hypergraph  * h1, 
                    Hypergraph  * h2, 
                    List       ** out_neighbors, 
                    List       ** in_neighbors, 
                    INT         * map, 
                    INT           k, 
                    INT           epsilon) 
{

  INT nvp    = 0;                 /* number of vertices in the contracted RBH  */
  INT nep    = 0;                 /* number of hyperarcs in the contracted RBH */
  INT nv     = h1->i_vertices;    /* number of vertices                        */
  INT ne     = h1->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h1->i_reds;        /* number of red vertices                    */
  INT np     = h1->i_pins;        /* number of pins                            */
  INT nwv    = h1->i_weights;     /* number of vertices weight                 */
  INT * reds = h1->ti_reds;       /* arrays of red vertices                    */

  // A limit for each weight to constraint vertices contration.
  INT * limit_weight = (INT *) malloc(sizeof(INT) * nwv);
  MEM_ERROR(limit_weight);

  // Initialization of limit weigth array to 0.
  for (INT wi = 0; wi < nwv; wi++) 
    {
      limit_weight[wi] = 0;
    }
  
  bool *is_red = h1->is_red;

  // Array of weight for each vertices.
  INT * weights = malloc(sizeof(INT) * nv * nwv);
  MEM_ERROR(weights);  

  // Array of sizes for each vertices.
  INT * sizes   =  (INT *) malloc(sizeof(INT) * nv);
  MEM_ERROR(sizes);

  // Boolean array : true iff a vertex is locked.
  bool * locked =  (bool*) malloc(sizeof(bool)*nv);
  MEM_ERROR(locked);

  INT dmax_ = 0;
  for (INT i = 0; i < nv; i++) {

    map[i]    = i;
    sizes[i]  = 0;
    locked[i] = false;

    for (INT wi=0; wi<nwv; wi++) {

      limit_weight[wi] += h1->ti_weights[i*nwv+wi];
      weights[i*nwv+wi]   = 0;

    }

    if (out_neighbors[i]->size+in_neighbors[i]->size>dmax_) {
      dmax_=out_neighbors[i]->size+in_neighbors[i]->size;
    }
  }

  for (INT wi=0; wi<nwv; wi++) {

    limit_weight[wi] = (limit_weight[wi]/k) *(1+epsilon/100);

  }

  INT *in_idx = (INT*)malloc(sizeof(INT)*(nv+1));
  MEM_ERROR(in_idx);
  in_idx[0]=0;
  for(INT u=0;u<nv;u++) in_idx[u+1]=in_idx[u]+in_neighbors[u]->size;
  INT *r_star_flat=(INT*)calloc(in_idx[nv]?in_idx[nv]:1, sizeof(INT));
  MEM_ERROR(r_star_flat);
  List *cell_in_neighbors;
  for (INT u = 0; u < nv; u++) 
    {
      INT w_hat = 0;
      
      // take the in neighbour
      cell_in_neighbors = in_neighbors[u];
      for (INT x = 0; x < in_neighbors[u]->size; x++) 
        {
        
          //selection of vertex v
          INT v = cell_in_neighbors->i;
          if (h1->ti_criticalities_left[v] > w_hat) 
            {
              w_hat=h1->ti_criticalities_left[v];
            }
            cell_in_neighbors = cell_in_neighbors->next;
        }

    cell_in_neighbors = in_neighbors[u];
    for (INT x = 0; x < in_neighbors[u]->size; x++) {
        INT v   = cell_in_neighbors->i;
        INT pos = in_idx[u] + x;
        r_star_flat[pos] = MAX(h1->ti_criticalities_right[v]-(w_hat-h1->ti_criticalities_left[v]), r_star_flat[pos]);
        if (is_red[v] || is_red[u])
            r_star_flat[pos] = MAX(h1->ti_criticalities_right[v], w_hat);
        cell_in_neighbors = cell_in_neighbors->next;
    }
  }

    INT * sorted_by_criticality = (INT *)malloc(sizeof(INT) * h1->i_vertices);
    
    for (INT i = 0; i < nv; i++) 
      {
        sorted_by_criticality[i] = i;
      }

    quick_sort(h1->ti_criticalities_right, sorted_by_criticality, 0, nv-1);

    bool change;
    bool correct_size;
    INT candidate;
    INT crit_candidate;

    for (INT i = 0; i < nv; i++) 
      {
      
        /* heavy edge matching */
        /* neighbor fusion     */

        INT u          = sorted_by_criticality[nv - i - 1];
        candidate      = -1;
        crit_candidate = -1;
        change         = false;

        if (!locked[u]) 
          {
            locked[u]   = true;
            map[u]      = nvp;
            change      = true;
            sizes[nvp] += 1;
            for (INT wi = 0; wi < nwv; wi++) 
              {
                weights[nvp * nwv + wi] += h1->ti_weights[u * nwv + wi];
              }

            List *mc=in_neighbors[u];
            for(INT x=0;x<in_neighbors[u]->size;x++){
              INT v=mc->i, r_v=r_star_flat[in_idx[u]+x];
              correct_size=true;
              for(INT wi=0;wi<nwv;wi++)
                if(weights[nvp*nwv+wi]+h1->ti_weights[v*nwv+wi]>=limit_weight[wi]) correct_size=false;
              if(!locked[v]&&correct_size&&r_v>crit_candidate){crit_candidate=r_v;candidate=v;}
              else if(!locked[v]&&correct_size&&r_v>=crit_candidate&&raisin_rng_bounded(2)>=1){crit_candidate=r_v;candidate=v;}
              mc=mc->next;
            }

            if (candidate == -1) 
              {
                List * cell = out_neighbors[u];
                
                for (INT x = 0; x < out_neighbors[u]->size; x++) 
                  {
                    //selection of vertex v
                    INT v   = cell->i;
                    correct_size = true;
                    
                    for (INT wi = 0; wi < nwv; wi++) 
                      {
                        if (weights[nvp * nwv + wi] + h1->ti_weights[v * nwv + wi] >= limit_weight[wi]) 
                          {
                            correct_size = false;
                          }
                      }
                    
                    if (!locked[v] && correct_size && h1->ti_criticalities_right[v] > crit_candidate) 
                      {
                        crit_candidate = h1->ti_criticalities_right[v];
                        candidate      = v;
                      } 
                    else if (!locked[v] && correct_size &&
                             h1->ti_criticalities_right[v] >= crit_candidate &&
                             raisin_rng_bounded(2) >= 1)
                      {
                        crit_candidate = h1->ti_criticalities_right[v];
                        candidate      = v;
                      }
                    
                    cell = cell->next;
                  }
              }
          }
        
        if (candidate > -1) 
          {

            map[candidate]      = nvp;
            locked[candidate]   = true;
            sizes[nvp]         += 1;
            
            for (INT wi = 0; wi < nwv; wi++) 
              {
                weights[nvp * nwv + wi] += h1->ti_weights[candidate * nwv + wi];
              }
          }
        
        if (change) 
          nvp++;
      

    }

  INT * ti_hyperedges_p = (INT*) malloc(sizeof(INT)*(np+2*ne));
  MEM_ERROR(ti_hyperedges_p);

  INT * ti_idx_hyperedges_p = (INT*) malloc(sizeof(INT)*ne);
  MEM_ERROR(ti_idx_hyperedges_p);
  
  INT * pins_jp = (INT*) calloc(nvp, sizeof(INT));
  MEM_ERROR(pins_jp);
  
  INT idx_jp  = 0;
  INT i_redsp = 0;
  INT i_pinsp = 0;
  INT idx_j;
  INT weight_j;
  INT size_jp, size_j;
  INT u, v;
  bool is_in_pins_j;

  for (INT j = 0; j < ne; j++) 
    {

      idx_j            = h1->ti_idx_hyperedges[j];
      weight_j         = h1->ti_hyperedges[idx_j];
      size_j           = h1->ti_hyperedges[idx_j+1];

      size_jp          = 0;
    
      for (INT i = 0; i < size_j; i++) 
        {

          is_in_pins_j  = false;
          u             = h1->ti_hyperedges[i + idx_j + 2];
      
          for (INT x = 0; x < size_jp; x++) 
            {
              v           = pins_jp[x];
        
              if (v == map[u]) 
                {
                  is_in_pins_j = true;
                }
            }
      
          if (!is_in_pins_j) 
            {
              pins_jp[size_jp++] = map[u];
            }
        }

    if (size_jp > 1) 
    {
      i_pinsp  += size_jp;
      ti_idx_hyperedges_p[nep++] = idx_jp;
      ti_hyperedges_p[idx_jp++] = weight_j;
      ti_hyperedges_p[idx_jp++] = size_jp;

      for (INT x = 0; x < size_jp; x++) 
        {
          ti_hyperedges_p[idx_jp++] = pins_jp[x];
        }
    }
  }

  for (INT i = 0; i < nvp; i++) 
    {
      is_red[i] = false;
    }

  for (INT i = 0; i < nr; i++) 
    {
      is_red[map[reds[i]]] = true;
    }

  for (INT i = 0; i < nvp; i++) 
    {
      if (is_red[i]) 
        {
          i_redsp++;
        }
    }
  
  h2->i_vertices       = nvp;
  h2->i_hyperedges     = nep;
  h2->i_reds           = i_redsp;
  h2->i_pins           = i_pinsp;
  h2->i_weights        = nwv;

  rbh_init(h2);

  INT idx_r = 0;

  for (INT v = 0; v < nv; v++) 
    {
      INT u = map[v];

      for (INT wi = 0; wi < nwv; wi++) 
        {
          h2->ti_weights[u * nwv + wi] += h1->ti_weights[v * nwv + wi];
        }
    
      if (h1->ti_delays[v] > h2->ti_delays[u]) 
        {
          h2->ti_delays[u] = h1->ti_delays[v];
        }
    
      if (h1->ti_criticalities_right[v] > h2->ti_criticalities_right[u]) 
        {
          h2->ti_criticalities_right[u] = h1->ti_criticalities_right[v];
        }
    
      if (is_red[u]) 
        {
          h2->ti_reds[idx_r++] = u ;
          is_red[u]            = false;
        }
    }
  
  for (INT j = 0; j < idx_jp; j++) 
    {
      h2->ti_hyperedges[j] = ti_hyperedges_p[j];
    }
  
  for (INT j = 0; j < nep; j++) 
    {
      h2->ti_idx_hyperedges[j] = ti_idx_hyperedges_p[j];
    }
  
  h2->s_rbh_name=(char*)malloc(sizeof(char)*(strlen(h1->s_rbh_name)+1));
  MEM_ERROR(h2->s_rbh_name);
  strcpy(h2->s_rbh_name,h1->s_rbh_name);
  h2->is_red=(bool*)calloc(h2->i_vertices,sizeof(bool));
  MEM_ERROR(h2->is_red);
  for(INT _i=0;_i<h2->i_reds;_i++) h2->is_red[h2->ti_reds[_i]]=true;
  h2->vth=vth_build(h2); MEM_ERROR(h2->vth);

  /* free section */
  free(r_star_flat);
  free(in_idx);
  free(pins_jp);
  free(weights);
  free(sizes);
  free(locked);
  free(ti_hyperedges_p);
  free(ti_idx_hyperedges_p);
  free(limit_weight);
  free(sorted_by_criticality);

  return 0;

}

/**
 * @brief Compute a subhypergraph (h2) of (h1) using a subset of 
 * vertices clusters.
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
compute_subhypergraph(Hypergraph * h1, 
                      Hypergraph * h2, 
                      INT        * cluster, 
                      INT          nvp) 
{

  INT nep    = 0;
  INT nrp    = 0;
  INT npinsp = 0;
  INT nv     = h1->i_vertices;    /* number of vertices                        */
  INT ne     = h1->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h1->i_reds;        /* number of red vertices                    */
  INT np     = h1->i_pins;        /* number of pins                            */
  INT nwv    = h1->i_weights;     /* number of vertices weight                 */
  INT * reds = h1->ti_reds;       /* arrays of red vertices                    */
    
  INT * h2_ti_hyperedges = (INT *)malloc(sizeof(INT) * (2 * ne + np));
  MEM_ERROR(h2_ti_hyperedges);

  INT * h2_ti_idx_hyperedges = (INT *)malloc(sizeof(INT) * ne);
  MEM_ERROR(h2_ti_idx_hyperedges);

  bool * is_in = (bool *)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_in);

  bool * is_red = (bool *)malloc(sizeof(bool) * nv);
  MEM_ERROR(is_red);

  bool * is_red2 = (bool *)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_red2);

  INT * toh2 = (INT *)malloc(sizeof(INT) * nvp);
  MEM_ERROR(toh2);

  bool * in_cluster = (bool *)malloc(sizeof(bool) * nv);
  MEM_ERROR(in_cluster);
  
  for (INT i = 0; i < nv; i++)
    {
      in_cluster[i] = false;
      is_red[i]     = false;
    }
  
  for (INT i = 0; i < nr; i++)
    {
      is_red[reds[i]]     = true;
    }
    
  for (INT i = 0; i < nvp; i++)
    {
      is_in[i]               = false;
      toh2[cluster[i]]       = i;
      in_cluster[cluster[i]] = true;
      is_red2[i]             = false;
    }
    
  for (INT i = 0; i < nv; i++) 
   {
     is_red2[cluster[i]] = (is_red[i] || is_red2[cluster[i]]);
   }

  for (INT i = 0; i < nvp; i++) 
    {
      if (is_red2[i]) 
        {
          nrp++;
        }
    }

    INT idx_j2 = 0;
    for (INT j = 0; j < ne; j++) 
      {
        INT idx_j    = h1->ti_idx_hyperedges[j];
        INT size_j   = h1->ti_hyperedges[idx_j + 1];

        INT source   = h1->ti_hyperedges[idx_j + 2];
        bool deleted = true;

        for (INT i = 0; i < size_j - 1 &&! deleted; i++) 
          {
            INT sink = h1->ti_hyperedges[idx_j + 3 + i];
            
            if (cluster[sink] != cluster[source]) 
              {
                deleted = false;
              }
          }
        
        if (!deleted) 
          {
            h2_ti_idx_hyperedges[nep++] = idx_j2;
            INT weight                  = h1->ti_hyperedges[idx_j];
            INT size_j2                 = 1;
            INT idx_size                = 0;
            is_in[cluster[source]]      = true;

            h2_ti_hyperedges[idx_j2++] = weight;
            idx_size                   = idx_j2;
            h2_ti_hyperedges[idx_j2++] = size_j2;
            h2_ti_hyperedges[idx_j2++] = cluster[source];

            for (INT i = 0; i < size_j - 1; i++) 
              {
                INT sink = h1->ti_hyperedges[idx_j + 3 + i];
                if (!is_in[cluster[sink]]) 
                  {
                    is_in[cluster[sink]]       = true;
                    h2_ti_hyperedges[idx_j2++] = cluster[sink];
                    size_j2++;
                  }
              }
            
            h2_ti_hyperedges[idx_size] = size_j2;
            npinsp += size_j2;
            
            for (INT i = 0; i < size_j2; i++) 
              {
                INT u    = h2_ti_hyperedges[idx_j2 - (size_j2) + i];
                is_in[u] = false;
              }
         }
    }

    h2->i_vertices       = nvp;
    h2->i_hyperedges     = nep;
    h2->i_reds           = nrp;
    h2->i_pins           = npinsp;
    h2->i_weights        = nwv;

    rbh_init(h2);

    h2->s_rbh_name = (char*)malloc(sizeof(char) * (strlen(h1->s_rbh_name) + 1));
    strcpy(h2->s_rbh_name, h1->s_rbh_name);

    // create subhypergraph
    for (INT j = 0; j < nep; j++) 
      {
        INT idx_j                    = h2_ti_idx_hyperedges[j];
        h2->ti_idx_hyperedges[j]     = idx_j;
        INT weight                   = h2_ti_hyperedges[idx_j];
        h2->ti_hyperedges[idx_j]     = weight;
        INT len_j                    = h2_ti_hyperedges[idx_j + 1];
        h2->ti_hyperedges[idx_j + 1] = len_j;

        for (INT i = 0; i < len_j; i++) 
          {
            INT u                              = h2_ti_hyperedges[idx_j + 2 + i];
            h2->ti_hyperedges[idx_j + 2 + i]   = u;
          }
      }
    
    INT idx_r = 0;
    
    for (INT i = 0; i < nvp; i++) 
      {
        if (is_red2[i]) 
          {
            h2->ti_reds[idx_r++] = i;
          }
      }  

    for (INT i = 0; i < nv; i++) 
      {
        h2->ti_delays[cluster[i]]              = MAX(h2->ti_delays[cluster[i]] , h1->ti_delays[i]);
        h2->ti_criticalities_right[cluster[i]] = MAX(h2->ti_criticalities_right[cluster[i]], h1->ti_criticalities_right[i]);
        h2->ti_criticalities_left[cluster[i]]  = MAX(h2->ti_criticalities_left[cluster[i]], h1->ti_criticalities_left[i]);

        for (INT wi = 0; wi < nwv; wi++) 
          {
            h2->ti_weights[cluster[i] * nwv + wi] += h1->ti_weights[i * nwv + wi];
          }
      }  

    /* free section */
    free(in_cluster);
    free(toh2);
    free(is_red2);
    free(is_red);
    free(is_in);
    free(h2_ti_idx_hyperedges);
    free(h2_ti_hyperedges);

    return 0;
}

/**
 * @brief Compute a reduced hypergraph (h2) according to clustering of 
 * an input red-black hypergraph (h1).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
compute_hypergraph_kclustering(Hypergraph * h1, 
                               Hypergraph * h2, 
                               List      ** out_neighbors, 
                               List      ** in_neighbors,
                              INT         * sort, 
                              INT         * cluster, 
                              INT           nvp) 
{
  RAISIN_UNUSED(out_neighbors);
  RAISIN_UNUSED(in_neighbors);
  RAISIN_UNUSED(sort);

  INT nep    = 0;
  INT nrp    = 0;
  INT npinsp = 0;
  INT nv     = h1->i_vertices;    /* number of vertices                        */
  INT ne     = h1->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h1->i_reds;        /* number of red vertices                    */
  INT np     = h1->i_pins;        /* number of pins                            */
  INT nwv    = h1->i_weights;     /* number of vertices weight                 */
  INT * reds = h1->ti_reds;       /* arrays of red vertices                    */

  INT * h2_ti_hyperedges = (INT *)malloc(sizeof(INT) * (2 * ne + np));
  MEM_ERROR(h2_ti_hyperedges);

  INT * h2_ti_idx_hyperedges = (INT *)malloc(sizeof(INT) * ne);
  MEM_ERROR(h2_ti_idx_hyperedges);

  bool * is_in = (bool *)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_in);

  bool * is_red = (bool *)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_red);

  for (INT i = 0; i < nvp; i++) 
    {
      is_in[i]  = false;
      is_red[i] = false;
    }

  for (INT i = 0; i < nr; i++) 
    {
      if (!is_red[cluster[reds[i]]]) 
        {
          is_red[cluster[reds[i]]] = true;
        }
    }

  for (INT i = 0; i < nvp; i++) 
    {
      if (is_red[i]) 
        {
          nrp++;
        }
    }

  // cluster [u] --> v, u \in h1, v \in h2
  // if cluster[sink] != cluster[source] => ok
  INT idx_j2 = 0;
  
  for (INT j = 0; j < ne; j++) 
    {
      INT idx_j    = h1->ti_idx_hyperedges[j];
      INT size_j   = h1->ti_hyperedges[idx_j + 1];
      INT source   = h1->ti_hyperedges[idx_j + 2];
      bool deleted = true;
        
      for (INT i = 0; i < size_j - 1; i++) 
        {
          INT sink = h1->ti_hyperedges[idx_j + 3 + i];
          
          if (cluster[source] != cluster[sink]) 
            {
              deleted = false;
            }
        }
        
        if (!deleted) 
          {
            h2_ti_idx_hyperedges[nep++] = idx_j2;
            INT weight                  = h1->ti_hyperedges[idx_j];
            INT size_j2                 = 1;
            INT idx_size                = 0;
            is_in[cluster[source]]      = true;

            h2_ti_hyperedges[idx_j2++] = weight;
            idx_size                   = idx_j2;
            h2_ti_hyperedges[idx_j2++] = size_j2;
            h2_ti_hyperedges[idx_j2++] = cluster[source];

            for (INT i = 0; i < size_j - 1; i++) 
              {
                INT sink = h1->ti_hyperedges[idx_j + 3 + i];
                
                if (!is_in[cluster[sink]]) 
                  {
                    is_in[cluster[sink]]       = true;
                    h2_ti_hyperedges[idx_j2++] = cluster[sink];
                    size_j2++;
                  }
              }
            
            h2_ti_hyperedges[idx_size] = size_j2;
            npinsp += size_j2;
            
            for (INT i = 0; i < size_j2; i++) 
              {
                INT u    = h2_ti_hyperedges[idx_j2 - (size_j2) + i];
                is_in[u] = false;
              }
          }
      }

    /* Compute cluster delays */
    INT * cluster_delays = (INT*)malloc(sizeof(INT)*nvp);
    MEM_ERROR(cluster_delays);

    bool * is_in_cluster = (bool*)malloc(sizeof(bool)*nv);
    MEM_ERROR(is_in_cluster);

    #if SPEED_DELAY_CLUSTERING
        ;
    #else
      for (INT i = 0; i < nvp; i++) 
        {
          for (INT u = 0; u < nv; u++) 
            {
              is_in_cluster[u] = false;
              
              if (cluster[u] == i) 
                {
                  is_in_cluster[u] = true;
                }
            }
          
          cluster_delays[i] = compute_subpmax(h1, neighbors, in_neighbors, sort, is_in_cluster);

        }
    #endif

    h2->i_vertices       = nvp;
    h2->i_hyperedges     = nep;
    h2->i_reds           = nrp;
    h2->i_pins           = npinsp;
    h2->i_weights        = nwv;

    rbh_init(h2);

    h2->s_rbh_name = (char *)malloc(sizeof(char) * (strlen(h1->s_rbh_name) + 1));
    strcpy(h2->s_rbh_name, h1->s_rbh_name);

    // create clustered hypergraph

    for (INT j = 0; j < nep; j++) 
      {
        INT idx_j                    = h2_ti_idx_hyperedges[j];
        h2->ti_idx_hyperedges[j]     = idx_j;
        INT weight                   = h2_ti_hyperedges[idx_j];
        h2->ti_hyperedges[idx_j]     = weight;
        INT len_j                    = h2_ti_hyperedges[idx_j + 1];
        h2->ti_hyperedges[idx_j + 1] = len_j;
        
        for (INT i = 0; i < len_j; i++) 
          {
            INT u                             = h2_ti_hyperedges[idx_j + 2 + i];
            h2->ti_hyperedges[idx_j + 2 + i]  = u;
          }
      }
    
    INT idx_r = 0;
    
    for (INT i = 0; i < nvp; i++) 
      {
        if (is_red[i]) 
          {
            h2->ti_reds[idx_r++] = i;
          }
      }

    for (INT i = 0; i < nv; i++) 
      {
        #if SPEED_DELAY_CLUSTERING
          h2->ti_delays[cluster[i]] = MAX(h2->ti_delays[cluster[i]], h1->ti_delays[i]);
        #else
          h2->ti_delays[cluster[i]] = cluster_delays[cluster[i]];
        #endif
          h2->ti_criticalities_right[cluster[i]] = MAX(h2->ti_criticalities_right[cluster[i]], h1->ti_criticalities_right[i]);
          h2->ti_criticalities_left[cluster[i]]  = MAX(h2->ti_criticalities_left[cluster[i]], h1->ti_criticalities_left[i]);
        
          for (INT wi = 0; wi < nwv; wi++) 
            {
              h2->ti_weights[cluster[i] * nwv + wi] += h1->ti_weights[i * nwv + wi];
            }
      }

    /* free section */

    free(is_in_cluster);
    free(cluster_delays);
    free(is_red);
    free(is_in);
    free(h2_ti_idx_hyperedges);
    free(h2_ti_hyperedges);

    return 0;
}

/**
 * @brief Compute a reduced hypergraph (h2) according to matching of 
 * an input red-black hypergraph (h1).
 * 
 * @attention This algorithm only work for matching (cluster of size 2 or 1)  
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
compute_hypergraph_clustering(Hypergraph * h1, 
                              Hypergraph * h2, 
                              INT        * cluster, 
                              INT          nvp) 
{
    
  INT nep    = 0;
  INT nrp    = 0;
  INT npinsp = 0;
  INT nv     = h1->i_vertices;    /* number of vertices                        */
  INT ne     = h1->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h1->i_reds;        /* number of red vertices                    */
  INT np     = h1->i_pins;        /* number of pins                            */
  INT nwv    = h1->i_weights;     /* number of vertices weight                 */
  INT * reds = h1->ti_reds;       /* arrays of red vertices                    */

  INT * h2_ti_hyperedges = (INT *)malloc(sizeof(INT) * (2 * ne + np));
  MEM_ERROR(h2_ti_hyperedges);

  INT * h2_ti_idx_hyperedges = (INT*)malloc(sizeof(INT) * ne);
  MEM_ERROR(h2_ti_idx_hyperedges);

  bool * is_in = (bool*)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_in);

  bool * is_red = (bool*)malloc(sizeof(bool) * nvp);
  MEM_ERROR(is_red);

  for (INT i = 0; i < nvp; i++) 
    {
      is_in[i]  = false;
      is_red[i] = false;
    }

  for (INT i = 0; i < nr; i++) 
    {
      is_red[cluster[reds[i]]] = true;
    }

  for (INT i = 0; i < nvp; i++) 
    {
      if (is_red[i]) 
        {
          nrp++;
        }
    }

  // cluster [u] --> v, u \in h1, v \in h2
  // if cluster[sink] != cluster[source] => ok
  INT idx_j2 = 0;
  for (INT j = 0; j < ne; j++) 
    {
      INT idx_j    = h1->ti_idx_hyperedges[j];
      INT size_j   = h1->ti_hyperedges[idx_j + 1];
      INT source   = h1->ti_hyperedges[idx_j + 2];
      bool deleted = true;
      
      for (INT i = 0; i < size_j - 1; i++) 
        {
          INT sink = h1->ti_hyperedges[idx_j + 3 + i];
          
          if (cluster[source] != cluster[sink]) 
            {
              deleted = false;
            }
        }
      
      if (!deleted) 
        {
          h2_ti_idx_hyperedges[nep++] = idx_j2;
          INT weight                  = h1->ti_hyperedges[idx_j];
          INT size_j2                 = 1;
          INT idx_size                = 0;
          is_in[cluster[source]]      = true;
          h2_ti_hyperedges[idx_j2++]  = weight;
          idx_size                    = idx_j2;
          h2_ti_hyperedges[idx_j2++]  = size_j2;
          h2_ti_hyperedges[idx_j2++]  = cluster[source];
          
          for (INT i = 0; i < size_j - 1; i++) 
            {
              INT sink = h1->ti_hyperedges[idx_j + 3 + i];
              
              if (!is_in[cluster[sink]]) 
                {
                  is_in[cluster[sink]]       = true;
                  h2_ti_hyperedges[idx_j2++] = cluster[sink];
                  size_j2++;
                }
             }
            
            h2_ti_hyperedges[idx_size] = size_j2;
            npinsp                    += size_j2;

            for (INT i = 0; i < size_j2; i++) 
              {
                INT u    = h2_ti_hyperedges[idx_j2 - (size_j2) + i];
                is_in[u] = false;
              }
          }
    }

    h2->i_vertices       = nvp;
    h2->i_hyperedges     = nep;
    h2->i_reds           = nrp;
    h2->i_pins           = npinsp;
    h2->i_weights        = nwv;

    rbh_init(h2);

    h2->s_rbh_name = (char *)malloc(sizeof(char) * (strlen(h1->s_rbh_name) + 1));
    strcpy(h2->s_rbh_name, h1->s_rbh_name);

    // create subhypergraph
    for (INT j = 0; j < nep; j++) 
      {
        INT idx_j                    = h2_ti_idx_hyperedges[j];
        h2->ti_idx_hyperedges[j]     = idx_j;
        INT weight                   = h2_ti_hyperedges[idx_j];
        h2->ti_hyperedges[idx_j]     = weight;
        INT len_j                    = h2_ti_hyperedges[idx_j + 1];
        h2->ti_hyperedges[idx_j + 1]   = len_j;
        for (INT i = 0; i < len_j; i++) 
          {
            INT u                              = h2_ti_hyperedges[idx_j + 2 + i];
            h2->ti_hyperedges[idx_j + 2 + i]   = u;
          }
      }
    
    INT idx_r = 0;

    for (INT i = 0; i < nvp; i++) 
      {
        if (is_red[i]) 
          {
            h2->ti_reds[idx_r++] = i;
          }
      }

    for (INT i = 0; i < nv; i++) 
      {
        #if SPEED_DELAY_CLUSTERING
          h2->ti_delays[cluster[i]]            = MAX(h2->ti_delays[cluster[i]], h1->ti_delays[i]);
        #else
          h2->ti_delays[cluster[i]]           += h1->ti_delays[i];
        #endif

        h2->ti_criticalities_right[cluster[i]] = MAX(h2->ti_criticalities_right[cluster[i]], h1->ti_criticalities_right[i]);
        h2->ti_criticalities_left[cluster[i]]  = MAX(h2->ti_criticalities_left[cluster[i]], h1->ti_criticalities_left[i]);
        
        for (INT wi = 0; wi < nwv; wi++) 
          {
            h2->ti_weights[cluster[i] * nwv + wi] += h1->ti_weights[i * nwv + wi];
          }
      }

    /* free section */
    free(is_red);
    free(is_in);
    free(h2_ti_idx_hyperedges);
    free(h2_ti_hyperedges);

    return 0;
}

/**
 * @brief Compute a reduced hypergraph (h2) according to clustering of 
 * an input red-black hypergraph (h1).
 *        
 *
 * @param h1               Input hypergraph.
 * @param h2               Output hypergraph.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param cluster          Clusters.
 * @param nvp              Number of clusters.
 *
 * @return Return a flag.
 */
int 
best_phi_clustering(Hypergraph  *  h1, 
                    Hypergraph  *  h2,
                    Arch        *  a,
                    List       ** neighbors,
                    List       ** in_neighbors, 
                    INT        *  sort, 
                    INT           umap[],
                    INT           cluster_size, 
                    INT           epsilon, 
                    INT           relax) 
{
  RAISIN_UNUSED(sort);
  RAISIN_UNUSED(epsilon);
  
  bool verbose = false;

  INT nv     = h1->i_vertices;    /* number of vertices                        */
  INT ne     = h1->i_hyperedges;  /* number of hyperarcs                       */
  INT nr     = h1->i_reds;        /* number of red vertices                    */
  INT nwv    = h1->i_weights;     /* number of vertices weight                 */
  INT * reds = h1->ti_reds;       /* arrays of red vertices                    */

  bool * is_red = (bool*)malloc(sizeof(bool)*nv);
  MEM_ERROR(is_red);

  int crit_max = 0;

  for (INT i = 0; i < nv; i++) 
    {
      if (h1->ti_criticalities_right[i] > crit_max) 
        {
          crit_max = h1->ti_criticalities_right[i];
        }
    
      is_red[i] = false;
    }

  for (INT i = 0; i < nr; i++) 
    {
      is_red[reds[i]] = true;
    }

  INT D = 1;
  
  for (INT i = 0; i < a->i_m; i++) 
    {
      for (INT j = i; j < a->i_n; j++)
        {
          if (a->ti_delay[i * a->i_n + j] > D)
            {
        	    if (a->ti_delay[i * a->i_n + j] > D) 
                {
        		      D = a->ti_delay[i * a->i_n + j];
        	      }
            }
        }
    }

  INT Dp = (INT)( (float)D * ((float)relax / 100.0));
  D = D - (D / relax);
  D = D - Dp;

  INT start = 0;
  INT end   = ne * D + crit_max;
  
  if (verbose)
    printf("start : 0, end : %d\n", end);
  
  INT phi;
  INT best_phi = ne * D;
  RAISIN_UNUSED(best_phi);

  INT * map = (INT *)malloc(sizeof(INT) * nv);
  MEM_ERROR(map);
  
  INT * best_map = (INT *)malloc(sizeof(INT) * nv);
  MEM_ERROR(best_map);

  if (verbose) {
    printf("nv(%d), ne(%d)\n", nv, ne);
  }

  bool * locked = (bool *)malloc(sizeof(bool) * nv);
  MEM_ERROR(locked);

  bool * flag   = (bool *)malloc(sizeof(bool) * nv);
  MEM_ERROR(flag);

  bool * is_in = (bool *)malloc(sizeof(bool) * nv);
  MEM_ERROR(is_in);

  INT * subqueue = (INT *)malloc(sizeof(INT) * (nv));
  MEM_ERROR(subqueue);

  INT * queue    = (INT *)malloc(sizeof(INT) * nv);
  MEM_ERROR(queue);

  INT limit[nwv];

  INT * weights = (INT *)malloc(sizeof(INT) * (nv + 1) * nwv);
  MEM_ERROR(weights);

  INT * sizes = (INT *)malloc(sizeof(INT) * nv * nwv);
  MEM_ERROR(sizes);

  for (INT wi = 0; wi < nwv; wi++) 
    {
      limit[wi] = cluster_size;
    }
  
  INT dmax_ = 0;
  
  for (INT i = 0; i < nv; i++) 
    {
      best_map[i] = i;
      flag[i] = false;
    
      for (INT wi = 0; wi < nwv; wi++) 
        {
          sizes[i * nwv + wi] = h1->ti_weights[i * nwv + wi];
        }

      if (neighbors[i]->size + in_neighbors[i]->size > dmax_) 
        {
          dmax_ = neighbors[i]->size + in_neighbors[i]->size;
        }
    } 

  bool infeaseable;
  INT wmax[nwv];
  INT target_pmax = 0;

  INT *in_idx_bpc=(INT*)malloc(sizeof(INT)*(nv+1)); MEM_ERROR(in_idx_bpc);
  in_idx_bpc[0]=0;
  for(INT u=0;u<nv;u++) in_idx_bpc[u+1]=in_idx_bpc[u]+in_neighbors[u]->size;
  INT *rs_bpc=(INT*)calloc(in_idx_bpc[nv]?in_idx_bpc[nv]:1,sizeof(INT)); MEM_ERROR(rs_bpc);
  List *cell_in_neighbors;
  
  for (INT u = 0; u < nv; u++) 
    {
      INT w_hat = 0;
      // take the in neighbour
      cell_in_neighbors = in_neighbors[u];
      for (INT x = 0; x < in_neighbors[u]->size; x++) 
        {
          //selection of vertex v
          INT v = cell_in_neighbors->i;
          if (h1->ti_criticalities_left[v] > w_hat) 
            {
              w_hat = h1->ti_criticalities_left[v];
            }
            cell_in_neighbors = cell_in_neighbors->next;
        }

        cell_in_neighbors = in_neighbors[u];
        for (INT x = 0; x < in_neighbors[u]->size; x++) {
            INT v   = cell_in_neighbors->i;
            INT pos = in_idx_bpc[u] + x;
            rs_bpc[pos] = MAX(h1->ti_criticalities_right[v]-(w_hat-h1->ti_criticalities_left[v]), rs_bpc[pos]);
            if (is_red[v] || is_red[u])
                rs_bpc[pos] = MAX(h1->ti_criticalities_right[v], w_hat);
            cell_in_neighbors = cell_in_neighbors->next;
        }
    }

  while (end >= start) 
    {

      wmax[0]     = 0;
      target_pmax = 0;
      
      // initialization
      for (INT i = 0; i < nv; i++) 
        {
          map[i]=i;
          locked[i] = false;

          for (INT wi = 0; wi < nwv; wi++) 
            {
              weights[i * nwv + wi] = (h1->ti_weights[i * nwv + wi]);
            }
     }

    infeaseable = false;
    phi         = (start+end) / 2;
    
    if (verbose) {
        printf("pmax tentative : %d\n", phi);
    }

    for (INT u = 0; u < nv && !infeaseable; u++) 
      {
        List * cell        = in_neighbors[u];
        /* r_star_flat via rs_bpc/in_idx_bpc for u */
        for (INT o = 0; o < in_neighbors[u]->size; o++) 
          {
            //selection of vertex v
            INT v   = cell->i;
            INT r_v = rs_bpc[in_idx_bpc[u]+o];
            if (r_v + D >= phi && map[u] != map[v]) 
              {
                if (target_pmax < h1->ti_criticalities_right[v] + D) 
                  {
                    target_pmax = h1->ti_criticalities_right[v] + D;
                  }
                
                //fusion u-v
                if (locked[u] && locked[v]) 
                  {
                    if (verbose)
                        printf("map[%d](%d) and map[%d](%d) merged\n", u, v, map[u], map[v]);

                    for (INT wi = 0; wi < nwv; wi++) 
                      {
                        weights[map[u] * nwv + wi] += weights[map[v] * nwv + wi];
                        if (weights[map[u] * nwv + wi] > limit[wi]) 
                          {
                            infeaseable = true;
                            wmax[wi]    = weights[map[u] * nwv + wi];
                          }
                      }
                    
                    //merge u and v mapping
                    //map all vertices mapped to v in u
                    for (INT x = 0; x < nv; x++) 
                      {
                        if (map[x] == map[v] && v != x) 
                          {
                            map[x] = map[u];
                          }
                      }
                    map[v] = map[u];
                  }
                else if (locked[v]) 
                  {
                    if (verbose)
                        printf("%d to map[%d](%d)\n", u, v, map[v]);
                    
                    //if only v is locked, map u to v
                    for (INT wi = 0; wi < nwv; wi++) 
                      {
                        weights[map[v] * nwv + wi] += weights[map[u] * nwv + wi];
                        if (weights[map[v] * nwv + wi] > limit[wi]) 
                          {
                            infeaseable = true;
                            wmax[wi]    = weights[map[v] * nwv + wi];
                          }
                      }
                    map[u]    = map[v];
                    locked[u] = true;
                }else 
                  {
                    if (verbose)
                        printf("%d to %d\n", u, v);
                        
                    //if only u is locked, map v to u or
                    //if u and v doesn't locked, we map v to u
                    for (INT wi = 0; wi < nwv; wi++) 
                      {
                        weights[map[u] * nwv + wi] += weights[map[v] * nwv + wi];
                        if (weights[map[u] * nwv + wi] > limit[wi]) 
                          {
                            infeaseable = true;
                            wmax[wi]    = weights[map[u] * nwv + wi];
                          } 
                      }
                    map[v]    = map[u];
                    locked[v] = true;
                    locked[u] = true;
                 }
              }
            cell        = cell->next;
            
          }
      }
    if (infeaseable)
      {
        if (verbose)
            printf("phi(%d) infeaseable with max size %d (limit=%d)\n", phi, wmax[0], limit[0]);
        start = phi + 1;
      }
    else
      {
        if (verbose)
            printf("phi(%d) feasable with max size %d (limit=%d)\n", phi, wmax[0], limit[0]);
        
        end      = phi - 1;
        best_phi = phi;
        for (INT i = 0; i < nv; i++) 
          {
            best_map[i] = map[i];
            for (INT wi = 0; wi < nwv; wi++) 
              {
                sizes[i * nwv + wi] = weights[i * nwv + wi];
              }
            flag[i] = locked[i];
          }
      }
    }

  for (INT i = 0; i < nv; i++) 
    {
      locked[i] = flag[i];
    }
  
  free(map);
  free(weights);

  // merge cluster procedure
  start = 0;
  end   = 0;
  INT   s = 0, e = 0;

  for (INT i = 0; i < nv; i++) 
    {
      queue[i] = -1;
    }

  start = 0;
  end   = 0;
  INT u,v, r_v;

   // sort by criticality
   INT * sorted_by_criticality = (INT *)malloc(sizeof(INT) * h1->i_vertices);
   
   for (INT i = 0; i < nv; i++) 
     {
       sorted_by_criticality[i] = i;
     }

   quick_sort(h1->ti_criticalities_right, sorted_by_criticality, 0, nv-1);

   for (INT i = nv - 1; i >= 0; i--) 
    {
      u = sorted_by_criticality[i];
      if (!flag[u])
        queue[end++] = u;
    }

  while (start != end && start < end) 
    {
      u = queue[start++];
      flag[u]=true;
      e             = 0;
      s             = 0;
      subqueue[e++] = u;
      infeaseable   = false;
    
      while((s != e && s < e ) && (!infeaseable) && (neighbors[u]->size > 0 || in_neighbors[u]->size > 0)) 
        {
          u   = subqueue[s++];
          r_v = -1;
          v   = -1;
          
          List * cell = in_neighbors[u];
          for (INT o = 0; o < in_neighbors[u]->size; o++) {
              INT rs_val = rs_bpc[in_idx_bpc[u]+o];
              if (rs_val > r_v) { v = cell->i; r_v = rs_val; }
              cell = cell->next;
          }

          if (v < 0)
            continue;

          if (flag[u] && flag[v]) 
            {
              if (verbose)
                printf("map[%d](%d) and map[%d](%d) merged\n", u, v, best_map[u], best_map[v]);
        
              for (INT wi = 0; wi < nwv; wi++) 
                {
                  if (sizes[best_map[u] * nwv + wi] + sizes[best_map[v] * nwv + wi] <= limit[wi]) 
                    {
                      sizes[best_map[u] * nwv + wi] += sizes[best_map[v] * nwv + wi];
                      
                      for (INT x = 0; x < nv; x++) 
                        {
                          if (best_map[x] == best_map[v] && v != x) 
                            {
                              best_map[x] = best_map[u];
                            }
                        }
                      best_map[v] = best_map[u];
                    }
                  else
                    {
                      infeaseable = true;
                    }
                }
              //merge u and v mapping
              //map all vertices mapped to v in u
            }
          else if (flag[v]) 
            {
              if (verbose)
                printf("%d to map[%d](%d)\n", u, v, best_map[v]);
        
              //if only v is locked, map u to v
              for (INT wi = 0; wi < nwv; wi++) 
                {
                  if (sizes[best_map[u] * nwv + wi] + sizes[best_map[v] * nwv + wi] <= limit[wi]) 
                    {
                      sizes[best_map[u] * nwv + wi] += sizes[best_map[v] * nwv + wi];
                      best_map[u]                    = best_map[v];
                      flag[u]                        = true;
                    }
                  else
                    {
                      infeaseable = true;
                    }
                }
            }
          else 
            {
              if (verbose)
                printf("%d to %d\n", u, v);
        
              //if only u is locked, map v to u or
              //if u and v doesn't locked, we map v to u
              for (INT wi = 0; wi < nwv; wi++) 
                {
                  if (sizes[best_map[u] * nwv + wi] + sizes[best_map[v] * nwv + wi] <= limit[wi]) 
                    {
                      sizes[best_map[u] * nwv + wi] += sizes[best_map[v] * nwv + wi];
                      best_map[v]                    = best_map[u];
                      flag[v]                        = true;
                      flag[u]                        = true;
                    }
                  else
                    {
                      infeaseable = true;
                    }
                }
            }
          cell = in_neighbors[u];
          for (INT x = 0; x < in_neighbors[u]->size; x++) {
              if (rs_bpc[in_idx_bpc[u]+x] > r_v) { v = cell->i; subqueue[e++] = v; }
              cell = cell->next;
          }
    
          if (infeaseable) 
            {
              while(s != e) 
                {
                  u = subqueue[s++];
                  if (end < nv - 1)
                    {
                      queue[end++] = u;
                    }
                }
            }
        }
    }

  //refinement
  for (INT u = 0; u < nv; u++) 
    {
      bool fusionable;
      List * cell        = neighbors[u];
    
      for (INT o = 0; o < neighbors[u]->size; o++) 
        {
          //selection of vertex x
          INT x      = cell->i;
          fusionable = true;
        
          for (INT wi = 0; wi < nwv; wi++) 
            {
              if (x < nv && best_map[x] < nv &&
                  sizes[best_map[u] * nwv + wi] + sizes[best_map[x] * nwv + wi] > limit[wi]) 
                {
                  fusionable = false;
                }
            }
        
          if (fusionable) 
            {
              for (INT v = 0; v < nv; v++) 
                {
                  if ((v != x) && (best_map[x] == best_map[v])) 
                    {
                      best_map[v] = best_map[u];
                    }
                }
            
              for (INT wi = 0; wi < nwv; wi++) 
                {
                  sizes[best_map[u] * nwv + wi] += sizes[best_map[x] * nwv + wi];
                }
              best_map[x] = best_map[u];
          }
        cell       = cell->next;
      }
    cell        = in_neighbors[u];
    
    for (INT o = 0; o < in_neighbors[u]->size; o++) 
      { 
        //selection of vertex v
        INT x      = cell->i;
        fusionable = true;
        
        for (INT wi = 0; wi < nwv; wi++) 
          {
            if (x < nv && best_map[x] < nv &&
                sizes[best_map[u] * nwv + wi] + sizes[best_map[x] * nwv + wi] > limit[wi]) 
              {
                fusionable = false;
              }
          }
          
          if (fusionable) 
            {
              for (INT v = 0; v < nv; v++) 
                {
                  if ((v != x) && (best_map[x] == best_map[v])) 
                    {
                      best_map[v] = best_map[u];
                    }
                }
                
                for (INT wi = 0; wi < nwv; wi++) 
                  {
                    sizes[best_map[u] * nwv + wi] += sizes[best_map[x] * nwv + wi];
                  }
                best_map[x] = best_map[u];
              }
            cell       = cell->next;
         }

    }

  start = 0;
  end   = 0;

  unsigned int cluster_number = 0;

  for (INT i = 0; i < nv; i++)
    {
      is_in[i]=false;
    }

  for (INT i=0; i < nv; i++)
    {
      if (!is_in[best_map[i]])
       {
         is_in[best_map[i]] = true;
         INT clus = best_map[i];
      
         for (INT u = 0; u < nv; u++)
           {
             if (best_map[u] == clus)
               {
                 umap[u]=cluster_number;
               }
           }
        cluster_number++;
       }
    }

  compute_hypergraph_clustering(h1, h2, umap, cluster_number);

  for (INT i = 0; i < nv; i++)
    {
      assert(best_map[i] < nv);
      for (INT wi=0; wi < nwv; wi++)
        {
          if (sizes[best_map[i] * nwv + wi] > limit[wi])
            {
              printf("[ERROR] sizes[%d][%d] = %d \n", best_map[i], wi,sizes[best_map[i] * nwv + wi]);
            }
        }
    }

  /* free section */
  free(sorted_by_criticality);
  free(is_red);
  free(best_map);
  free(locked);
  free(flag);
  free(is_in);
  free(subqueue);
  free(queue);

  free(rs_bpc);
  free(in_idx_bpc);
  free(sizes);

  return 0;
}
