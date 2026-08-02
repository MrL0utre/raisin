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
/**   NAME       : fm.c                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are definitions for the     **/
/**                FM base refinement functions.           **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "fm.h"
#include "vth.h"
#include "bqueue.h"

/**
 * @brief Function implementing the DKFM algorithm defined and 
 * explained in the referenced paper available following this link : 
 *            https://doi.org/10.1007/978-3-031-36024-4_50 
 * The algorithm apply vertiex moves accross a partition in order
 * to reduce the maximum path degradation of the initial partition.
 * This algorithm is an adaptaion of existing FM strategy oriented 
 * for multi-partitionning. 
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
 *
 * @return Integer.
 */
INT 
dkfm(Hypergraph * h,  
     Arch       * a, 
     List      ** out_neighbors, 
     List      ** in_neighbors, 
     INT        * sort, 
     PART * partition, 
     INT          perform, 
     INT          tolerance,
     INT          k) 
{
    INT nv  = h->i_vertices;
    INT nwv = h->i_weights;

    bool  *is_in_halo = (bool *)calloc(nv, sizeof(bool));  MEM_ERROR(is_in_halo);
    bool  *is_locked  = (bool *)calloc(nv, sizeof(bool));  MEM_ERROR(is_locked);
    INT   *part_sizes = (INT  *)calloc(k * nwv, sizeof(INT)); MEM_ERROR(part_sizes);
    INT   *capacity_cst = (INT*)calloc(nwv, sizeof(INT));   MEM_ERROR(capacity_cst);
    float *cost       = (float*)malloc(sizeof(float) * nv); MEM_ERROR(cost);

    for (INT u = 0; u < nv; u++) {
        for (INT wi = 0; wi < nwv; wi++) {
            part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
            capacity_cst[wi]               += h->ti_weights[u*nwv+wi];
        }
        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            if (partition[u] != partition[cell->i]) is_in_halo[u] = true;
            cell = cell->next;
        }
        List *inc = in_neighbors[u];
        for (INT x = 0; x < in_neighbors[u]->size; x++) {
            if (partition[u] != partition[inc->i]) is_in_halo[u] = true;
            inc = inc->next;
        }
    }
    for (INT wi = 0; wi < nwv; wi++)
        capacity_cst[wi] = (INT)ceil((float)capacity_cst[wi] / (float)k * 1.05f);

    /* Pre-allocated criticality buffers — reused for every evaluation */
    INT  *_delays = (INT*)malloc(sizeof(INT)*nv); MEM_ERROR(_delays);
    bool *_is_red = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(_is_red);
    bool *_flag   = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(_flag);

    PART * partitionp = (PART*)calloc(nv, sizeof(PART)); MEM_ERROR(partitionp);
    /* min_crit: all-same partition (trivial lower bound) */
    INT min_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    for (INT u = 0; u < nv; u++) partitionp[u] = (PART)(u % k);
    INT max_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    INT cur_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition,  _delays, _is_red, _flag);

    INT min_gain = cur_crit_path - max_crit_path;
    INT max_gain = cur_crit_path - min_crit_path;
    if (max_gain <= min_gain) max_gain = min_gain + 1;  /* avoid zero-range */
    /* Clamp range to avoid huge BQueue allocation for large pmax graphs */
    if ((long long)max_gain - min_gain > 100000LL) {
        min_gain = cur_crit_path - 100000;
        max_gain = cur_crit_path;
    }

    for (INT u = 0; u < nv; u++) {
        partitionp[u] = partition[u];
        INT D = 1;
        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            INT v = cell->i;
            INT d = a->ti_delay[partition[u]*a->i_m + partition[v]];
            if (d > D) D = d;
            cell = cell->next;
        }
        cost[u] = ((float)tolerance / 100.0f) *
                  ((float)(max_crit_path - h->ti_criticalities_right[u]) / (float)D);
    }

    BQueue **bq = (BQueue**)malloc(sizeof(BQueue*) * k);
    MEM_ERROR(bq);
    for (INT p = 0; p < k; p++) {
        bq[p] = bq_alloc(nv, min_gain, max_gain);
        MEM_ERROR(bq[p]);
    }

    bool *is_computed = (bool*)malloc(sizeof(bool) * k); MEM_ERROR(is_computed);
    INT  anc_part, gain;

    for (INT i = 0; i < nv; i++) {
        if (!is_in_halo[i]) continue;
        for (INT p = 0; p < k; p++) is_computed[p] = false;
        INT ncomp = 0;

        List *cell = out_neighbors[i];
        for (INT x = 0; x < out_neighbors[i]->size && ncomp < k; x++) {
            INT v = cell->i;
            if (partition[i] != partition[v] && !is_computed[partition[v]]) {
                anc_part     = partition[i];
                partition[i] = partition[v];
                gain = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                bq[partition[i]]->gval[i] = gain;
                bq_insert(bq[partition[i]], i, gain);
                partition[i] = anc_part;
                is_computed[partition[v]] = true;
                ncomp++;
            }
            cell = cell->next;
        }
        List *inc = in_neighbors[i];
        for (INT x = 0; x < in_neighbors[i]->size && ncomp < k; x++) {
            INT v = inc->i;
            if (partition[i] != partition[v] && !is_computed[partition[v]]) {
                anc_part     = partition[i];
                partition[i] = partition[v];
                gain = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                bq[partition[v]]->gval[i] = gain;
                bq_insert(bq[partition[v]], i, gain);
                partition[i] = anc_part;
                is_computed[partition[v]] = true;
                ncomp++;
            }
            inc = inc->next;
        }
    }

    bool *to_balance = (bool*)calloc(k, sizeof(bool)); MEM_ERROR(to_balance);
    for (INT p = 0; p < k; p++)
        for (INT wi = 0; wi < nwv; wi++)
            if (part_sizes[p*nwv+wi] >= capacity_cst[wi]) to_balance[p] = true;

    INT base_cur   = cur_crit_path;
    INT best_crit  = cur_crit_path;
    PART * ex_part   = (PART*)malloc(sizeof(PART)*nv); MEM_ERROR(ex_part);
    memcpy(ex_part, partition, sizeof(PART)*nv);

    INT N_MOVES = (INT)ceil((float)nv*(1.0f/(float)perform));

    for (INT i = 0; i < N_MOVES; i++) {
        INT selected_part = -1;
        for (INT p = 0; p < k; p++)
            if (!to_balance[p] && !bq_empty(bq[p])) selected_part = p;
        if (selected_part == -1) break;

        INT u = bq_dequeue_max(bq[selected_part]);
        if (u == BQ_NONE) break;

        is_locked[u]   = true;
        cur_crit_path  = base_cur - bq[selected_part]->gval[u];

        for (INT wi = 0; wi < nwv; wi++) {
            part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
            part_sizes[partition[u]*nwv+wi]  -= h->ti_weights[u*nwv+wi];
            if (part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]) to_balance[selected_part] = true;
        }
        for (INT wi = 0; wi < nwv; wi++)
            to_balance[partition[u]] = (part_sizes[partition[u]*nwv+wi] > capacity_cst[wi]);

        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            INT v = cell->i;
            if (!is_locked[v]) {
                /* Check if v is still connected to old partition of u */
                bool still_connected = false;
                List *cp = out_neighbors[v];
                for (INT o = 0; o < out_neighbors[v]->size; o++) {
                    if (partition[cp->i] == partition[u]) { still_connected = true; break; }
                    cp = cp->next;
                }
                if (!still_connected && bq[partition[u]]->gval[v] != min_gain - 1) {
                    /* O(1) update-key via bucket queue */
                    bq_update(bq[partition[u]], v, min_gain);
                }
                if (partition[v] != selected_part) {
                    INT ancp = partition[v];
                    partition[u] = selected_part;
                    partition[v] = selected_part;
                    gain = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                    if (bq[selected_part]->gval[v] != min_gain - 1)
                        bq_update(bq[selected_part], v, gain);
                    else
                        bq_insert(bq[selected_part], v, gain);
                    partition[v] = ancp;
                }
            }
            cell = cell->next;
        }
        partition[u] = selected_part;
        cur_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);

        if (cur_crit_path < best_crit) {
            best_crit = cur_crit_path;
            memcpy(ex_part, partition, sizeof(PART)*nv);
        }
    }

    INT crit_p  = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    INT crit_ex = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, ex_part, _delays, _is_red, _flag);
    memcpy(partition, (crit_p < crit_ex) ? partitionp : ex_part, sizeof(PART)*nv);

    /* free */
    free(_delays); free(_is_red); free(_flag);
    free(is_in_halo); free(is_locked); free(partitionp); free(ex_part);
    for (INT p = 0; p < k; p++) bq_free(bq[p]);
    free(bq);
    free(to_balance); free(is_computed); free(cost);
    free(part_sizes); free(capacity_cst);

    return 0;   
}


/**
 * @brief Function implementing the fast version of the DKFM 
 * algorithm defined in the referenced document following 
 * this link : 
 *            https://theses.hal.science/tel-04731886/ 
 *
 * @param h                Input hypergraph.
 * @param a                Architecture.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Partition array.
 * @param perform          Percentage of moves. 
 * @param tolerance        Tolerance for cost degradation. 
 * @param k                Number of part. 
 *
 * @return INT.
 */
INT 
dkfm_fast(Hypergraph * h,  
         Arch       * a, 
         List      ** out_neighbors, 
         List      ** in_neighbors, 
         INT        * sort, 
         PART * partition, 
         INT          perform, 
         INT          tolerance,
         INT          k) 
{
    INT nv  = h->i_vertices;
    INT nwv = h->i_weights;

    bool  *is_in_halo  = (bool *)calloc(nv, sizeof(bool));  MEM_ERROR(is_in_halo);
    bool  *is_locked   = (bool *)calloc(nv, sizeof(bool));  MEM_ERROR(is_locked);
    INT   *part_sizes  = (INT  *)calloc(k * nwv, sizeof(INT)); MEM_ERROR(part_sizes);
    INT   *capacity_cst = (INT *)calloc(nwv, sizeof(INT));  MEM_ERROR(capacity_cst);
    float *cost        = (float*)malloc(sizeof(float) * nv); MEM_ERROR(cost);

    /* Build halo + part_sizes */
    for (INT u = 0; u < nv; u++) {
        for (INT wi = 0; wi < nwv; wi++) {
            part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
            capacity_cst[wi]               += h->ti_weights[u*nwv+wi];
        }
        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            if (partition[u] != partition[cell->i]) is_in_halo[u] = true;
            cell = cell->next;
        }
        List *inc = in_neighbors[u];
        for (INT x = 0; x < in_neighbors[u]->size; x++) {
            if (partition[u] != partition[inc->i]) is_in_halo[u] = true;
            inc = inc->next;
        }
    }
    for (INT wi = 0; wi < nwv; wi++)
        capacity_cst[wi] = (INT)ceil((float)capacity_cst[wi] / (float)k * 1.05f);

    /* Pre-allocated criticality buffers */
    INT  *_delays = (INT*)malloc(sizeof(INT)*nv); MEM_ERROR(_delays);
    bool *_is_red = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(_is_red);
    bool *_flag   = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(_flag);

    PART * partitionp = (PART*)calloc(nv, sizeof(PART)); MEM_ERROR(partitionp);
    INT min_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    for (INT u = 0; u < nv; u++) partitionp[u] = (PART)(u % k);
    INT max_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    INT cur_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition,  _delays, _is_red, _flag);

    INT min_gain = cur_crit_path - max_crit_path;
    INT max_gain = cur_crit_path - min_crit_path;
    if (max_gain <= min_gain) max_gain = min_gain + 1;
    /* Clamp range to avoid huge BQueue allocation for large pmax graphs */
    if ((long long)max_gain - min_gain > 100000LL) {
        min_gain = cur_crit_path - 100000;
        max_gain = cur_crit_path;
    }

    for (INT u = 0; u < nv; u++) {
        partitionp[u] = partition[u];
        INT D = 1;
        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            INT v = cell->i;
            INT d = a->ti_delay[partition[u]*a->i_m + partition[v]];
            if (d > D) D = d;
            cell = cell->next;
        }
        cost[u] = ((float)tolerance / 100.0f) *
                  ((float)(max_crit_path - h->ti_criticalities_right[u]) / (float)D);
    }

    BQueue **bq = (BQueue**)malloc(sizeof(BQueue*) * k); MEM_ERROR(bq);
    for (INT p = 0; p < k; p++) {
        bq[p] = bq_alloc(nv, min_gain, max_gain);
        MEM_ERROR(bq[p]);
    }

    /* Build halo array for random-selection variant */
    INT *halo = (INT*)malloc(sizeof(INT)*nv); MEM_ERROR(halo);
    INT  halo_size = 0;
    for (INT u = 0; u < nv; u++) if (is_in_halo[u]) halo[halo_size++] = u;

    bool *is_computed = (bool*)malloc(sizeof(bool) * k); MEM_ERROR(is_computed);

    /* Initial gain computation — randomly sample halo vertices */
    INT N_INIT = (INT)ceil((float)halo_size * (1.0f / (float)perform));
    for (INT ii = 0; ii < N_INIT && halo_size > 0; ii++) {
        INT idx = rand() % halo_size;
        INT i   = halo[idx];
        /* Remove from halo (swap-with-last) */
        halo[idx] = halo[--halo_size];

        for (INT p = 0; p < k; p++) is_computed[p] = false;
        INT ncomp = 0;

        List *cell = out_neighbors[i];
        for (INT x = 0; x < out_neighbors[i]->size && ncomp < k; x++) {
            INT v = cell->i;
            if (partition[i] != partition[v] && !is_computed[partition[v]]) {
                INT anc = partition[i];
                partition[i] = partition[v];
                INT gain = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                bq[partition[i]]->gval[i] = gain;
                bq_insert(bq[partition[i]], i, gain);
                partition[i] = anc;
                is_computed[partition[v]] = true;
                ncomp++;
            }
            cell = cell->next;
        }
        List *inc = in_neighbors[i];
        for (INT x = 0; x < in_neighbors[i]->size && ncomp < k; x++) {
            INT v = inc->i;
            if (partition[i] != partition[v] && !is_computed[partition[v]]) {
                INT anc = partition[i];
                partition[i] = partition[v];
                INT gain = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                bq[partition[v]]->gval[i] = gain;
                bq_insert(bq[partition[v]], i, gain);
                partition[i] = anc;
                is_computed[partition[v]] = true;
                ncomp++;
            }
            inc = inc->next;
        }
    }

    bool *to_balance = (bool*)calloc(k, sizeof(bool)); MEM_ERROR(to_balance);
    for (INT p = 0; p < k; p++)
        for (INT wi = 0; wi < nwv; wi++)
            if (part_sizes[p*nwv+wi] >= capacity_cst[wi]) to_balance[p] = true;

    INT base_cur  = cur_crit_path;
    INT best_crit = cur_crit_path;
    PART * ex_part  = (PART*)malloc(sizeof(PART)*nv); MEM_ERROR(ex_part);
    memcpy(ex_part, partition, sizeof(PART)*nv);

    INT N_MOVES = (INT)ceil((float)nv*(1.0f/(float)perform));

    for (INT i = 0; i < N_MOVES; i++) {
        INT selected_part = -1;
        for (INT p = 0; p < k; p++)
            if (!to_balance[p] && !bq_empty(bq[p])) selected_part = p;
        if (selected_part == -1) break;

        INT u = bq_dequeue_max(bq[selected_part]);
        if (u == BQ_NONE) break;

        is_locked[u]  = true;
        INT gain_u    = bq[selected_part]->gval[u];
        cur_crit_path = base_cur - gain_u;

        for (INT wi = 0; wi < nwv; wi++) {
            part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
            part_sizes[partition[u]*nwv+wi]  -= h->ti_weights[u*nwv+wi];
            if (part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]) to_balance[selected_part] = true;
        }
        for (INT wi = 0; wi < nwv; wi++)
            to_balance[partition[u]] = (part_sizes[partition[u]*nwv+wi] > capacity_cst[wi]);

        List *cell = out_neighbors[u];
        for (INT x = 0; x < out_neighbors[u]->size; x++) {
            INT v = cell->i;
            if (!is_locked[v]) {
                bool still_connected = false;
                List *cp = out_neighbors[v];
                for (INT o = 0; o < out_neighbors[v]->size; o++) {
                    if (partition[cp->i] == partition[u]) { still_connected = true; break; }
                    cp = cp->next;
                }
                if (!still_connected && bq[partition[u]]->gval[v] != min_gain - 1) {
                    bq_update(bq[partition[u]], v, min_gain);
                }
                if (partition[v] != selected_part) {
                    INT ancp = partition[v];
                    partition[u] = selected_part;
                    partition[v] = selected_part;
                    INT g = cur_crit_path - compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
                    if (bq[selected_part]->gval[v] != min_gain - 1)
                        bq_update(bq[selected_part], v, g);
                    else
                        bq_insert(bq[selected_part], v, g);
                    partition[v] = ancp;
                }
            }
            cell = cell->next;
        }
        partition[u] = selected_part;

        if ( (cost[u] <= 0) || (i % 100 == 99) ) {
            cur_crit_path = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partition, _delays, _is_red, _flag);
        }

        if (cur_crit_path < best_crit) {
            best_crit = cur_crit_path;
            memcpy(ex_part, partition, sizeof(PART)*nv);
        }
    }

    INT crit_p  = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, partitionp, _delays, _is_red, _flag);
    INT crit_ex = compute_partition_criticality_with_buf(h, a, out_neighbors, in_neighbors, sort, ex_part, _delays, _is_red, _flag);
    memcpy(partition, (crit_p < crit_ex) ? partitionp : ex_part, sizeof(PART)*nv);

    /* free */
    free(_delays); free(_is_red); free(_flag);
    free(is_in_halo); free(is_locked); free(partitionp); free(ex_part);
    for (INT p = 0; p < k; p++) bq_free(bq[p]);
    free(bq);
    free(to_balance); free(is_computed); free(cost); free(halo);
    free(part_sizes); free(capacity_cst);

    return 0;   
}


/**
 * @brief compute the connectivity cut cost according 
 *        to a partition and a target topology (architecture).
 *
 * @param h                Input hypergraph.
 * @param a                Target topology.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 *
 * @return The connectivity cut cost.
 */
INT
compute_partition_criticality_with_buf(Hypergraph *h,
                                       Arch       *a,
                                       List      **out_neighbors,
                                       List      **in_neighbors,
                                       INT        *sort,
                                       PART       *partition,
                                       INT        *delays,
                                       bool       *is_red,
                                       bool       *flag)
{
    RAISIN_UNUSED(in_neighbors);
    INT nv    = h->i_vertices;
    bool *isr = h->is_red;   /* precomputed static mask — no malloc */
    (void)is_red; (void)flag;

    for (INT i = 0; i < nv; i++)
        delays[i] = h->ti_delays[i];

    for (INT i = 0; i < nv; i++) {
        INT v  = sort[i];
        INT dv = isr[v] ? h->ti_delays[v] : delays[v];

        List *cell = out_neighbors[v];
        for (INT x = 0; x < out_neighbors[v]->size; x++) {
            INT u    = cell->i;
            INT base = dv + h->ti_delays[u];
            if (partition[u] != partition[v])
                base += a->ti_delay[partition[u] * a->i_m + partition[v]];
            if (base > delays[u]) delays[u] = base;
            cell = cell->next;
        }
    }

    INT max_delays = 0;
    for (INT v = 0; v < nv; v++)
        if (delays[v] > max_delays) max_delays = delays[v];

    return max_delays;
}

INT
compute_partition_criticality(Hypergraph * h,
                              Arch       * a,
                              List      ** out_neighbors,
                              List      ** in_neighbors,
                              INT        * sort,
                              PART       * partition)
{
    INT   nv     = h->i_vertices;
    INT  *delays = (INT  *)malloc(sizeof(INT)  * nv); MEM_ERROR(delays);
    bool *is_red = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(is_red);
    bool *flag   = (bool *)malloc(sizeof(bool) * nv); MEM_ERROR(flag);

    INT result = compute_partition_criticality_with_buf(
        h, a, out_neighbors, in_neighbors, sort, partition,
        delays, is_red, flag);

    free(delays);
    free(is_red);
    free(flag);
    return result;
}


INT 
compute_partition_cut(Hypergraph * h, 
                      Arch       * a, 
                      List      ** out_neighbors, 
                      List      ** in_neighbors, 
                      INT        * sort, 
                      PART * partition,
                      INT        * lambda,
                      INT          k)
{
  RAISIN_UNUSED(a);
  RAISIN_UNUSED(out_neighbors);
  RAISIN_UNUSED(in_neighbors);
  RAISIN_UNUSED(sort);

  INT ne     = h->i_hyperedges;  /* number of hyperarcs                       */

  
  INT * hyperedges     = h->ti_hyperedges;     /* array of hyperedges         */
  INT * idx_hyperedges = h->ti_idx_hyperedges; /* index of hyperedges         */

  INT cut_cost = 0;

  bool * is_in_part = (bool*)malloc(sizeof(bool) * k);
  MEM_ERROR(is_in_part);

  for (INT p = 0; p < k; p++) 
    {
      is_in_part[p] = false;
    }

  INT lambda_cost;

  for (INT j = 0; j < ne; j++) 
    {
      INT idx_j  = idx_hyperedges[j];
      INT wj     = hyperedges[idx_j];
      INT size_j = hyperedges[idx_j + 1];
      lambda[j]  = 0;
    
      for (INT i = 0; i < size_j; i++) 
        {
          INT e = hyperedges[idx_j + 2 + i];
          INT p = partition[e];
          is_in_part[p] = true;
        }

      lambda_cost = 0;
    
      for (INT p = 0; p < k; p++) 
        {
          if(is_in_part[p]) 
            {
              lambda_cost += 1;
              is_in_part[p] = false;
            }
        }

      lambda[j] = (lambda_cost - 1) * MAX(1,wj);
      cut_cost += (lambda_cost - 1) * MAX(1,wj);

    }
  
  free(is_in_part);

  return cut_cost;

}

/**
 * @brief Function implementing the K-FM algorithm: 
 * The algorithm apply vertiex moves accross a partition in order
 * to reduce the connectivity cut cost of the initial partition.
 * This algorithm is an adaptaion of existing FM strategy oriented 
 * for multi-partitionning.
 *
 * @param h                Input hypergraph.
 * @param a                Target topology.
 * @param out_neighbors    Outgoing adjacency lists.
 * @param in_neighbors     Incoming adjacency lists.
 * @param sort             Topological order of vertices.
 * @param partition        Cluster partition array.
 *
 * @return The connectivity cut cost.
 */
INT 
kfm(Hypergraph * h,  
    Arch       * a, 
    List      ** out_neighbors, 
    List      ** in_neighbors, 
    INT        * sort, 
    PART * partition, 
    INT          perform, 
    INT          tolerance,
    INT          k) 
{
  RAISIN_UNUSED(tolerance);
  INT nv  = h->i_vertices;
  INT ne  = h->i_hyperedges;
  INT nwv = h->i_weights;

  INT * hyperedges     = h->ti_hyperedges;
  INT * idx_hyperedges = h->ti_idx_hyperedges;

  VtH *vth = h->vth;

  INT * cnt = (INT*)calloc(ne * k, sizeof(INT));
  MEM_ERROR(cnt);

  /* Build initial cnt[] */
  for (INT j = 0; j < ne; j++) {
    INT base = idx_hyperedges[j];
    INT sz   = hyperedges[base + 1];
    for (INT i = 0; i < sz; i++)
      cnt[j * k + partition[hyperedges[base + 2 + i]]]++;
  }

  INT * neighbours_parts = (INT*)calloc(nv * k, sizeof(INT));
  MEM_ERROR(neighbours_parts);

  bool * is_in_halo = (bool*)calloc(nv, sizeof(bool));
  MEM_ERROR(is_in_halo);
  bool * is_locked  = (bool*)calloc(nv, sizeof(bool));
  MEM_ERROR(is_locked);

  INT * part_sizes  = (INT*)calloc(k * nwv, sizeof(INT));
  MEM_ERROR(part_sizes);
  INT * capacity_cst = (INT*)calloc(nwv, sizeof(INT));
  MEM_ERROR(capacity_cst);

  /* Compute lambda[j] = connectivity cost of arc j = (nb distinct parts - 1) * weight */
  INT * lambda = (INT*)calloc(ne, sizeof(INT));
  MEM_ERROR(lambda);
  for (INT j = 0; j < ne; j++) {
    INT base = idx_hyperedges[j];
    INT wj   = hyperedges[base];
    INT nparts = 0;
    for (INT p = 0; p < k; p++) if (cnt[j*k+p] > 0) nparts++;
    lambda[j] = (nparts - 1) * wj;
  }

  for (INT u = 0; u < nv; u++) {
    for (INT wi = 0; wi < nwv; wi++) {
      part_sizes[partition[u]*nwv+wi] += h->ti_weights[u*nwv+wi];
      capacity_cst[wi]               += h->ti_weights[u*nwv+wi];
    }
    List *cell = out_neighbors[u];
    for (INT x = 0; x < out_neighbors[u]->size; x++) {
      INT v = cell->i;
      if (partition[u] != partition[v]) is_in_halo[u] = true;
      neighbours_parts[u*k + partition[v]]++;
      cell = cell->next;
    }
    List *inc = in_neighbors[u];
    for (INT x = 0; x < in_neighbors[u]->size; x++) {
      INT v = inc->i;
      if (partition[u] != partition[v]) is_in_halo[u] = true;
      neighbours_parts[u*k + partition[v]]++;
      inc = inc->next;
    }
  }
  for (INT wi = 0; wi < nwv; wi++)
    capacity_cst[wi] = (INT)ceil((float)capacity_cst[wi] / (float)k * 1.05f);

  INT current_cut = compute_partition_cut(h, a, out_neighbors, in_neighbors, sort, partition, lambda, k);

  /* --- gain data structures: one pqueue per partition --- */
  INT ** partition_moves = (INT**)malloc(sizeof(INT*) * k);
  MEM_ERROR(partition_moves);
  INT ** partition_gains = (INT**)malloc(sizeof(INT*) * k);
  MEM_ERROR(partition_gains);
  INT * partition_gain_sizes = (INT*)calloc(k, sizeof(INT));
  MEM_ERROR(partition_gain_sizes);

  for (INT p = 0; p < k; p++) {
    partition_moves[p] = (INT*)calloc(nv, sizeof(INT));  MEM_ERROR(partition_moves[p]);
    partition_gains[p] = (INT*)malloc(sizeof(INT)*nv);   MEM_ERROR(partition_gains[p]);
    for (INT i = 0; i < nv; i++) partition_gains[p][i] = -ne;
  }

  bool *is_in_part = (bool*)calloc(k, sizeof(bool));
  MEM_ERROR(is_in_part);

  /*
   * Compute initial gain for each halo vertex v moving to part p:
   *
   *   gain(v, p) = current_cut - cost_if_v_moves_to_p
   *
   */
  for (INT i = 0; i < nv; i++) {
    if (!is_in_halo[i]) continue;
    INT p_cur = partition[i];
    for (INT p = 0; p < k; p++) {
      if (p == p_cur) continue;
      if (neighbours_parts[i*k + p] == 0) continue;

      INT delta = 0;
      const INT *je = vth_begin(vth, i);
      const INT *jend = vth_end(vth, i);
      for (; je != jend; je++) {
        INT j   = *je;
        INT wj  = hyperedges[idx_hyperedges[j]];
        /* nparts before move */
        INT nb_before = 0;
        for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp] > 0) nb_before++;
        /* simulate move */
        cnt[j*k + p_cur]--;
        cnt[j*k + p]++;
        INT nb_after = 0;
        for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp] > 0) nb_after++;
        delta += (nb_after - nb_before) * wj;
        /* undo */
        cnt[j*k + p_cur]++;
        cnt[j*k + p]--;
      }
      INT gain = -delta;   /* gain = cut_reduction */
      partition_gains[p][i] = gain;
      pqueue_add_element(partition_moves[p], partition_gains[p], i, partition_gain_sizes[p]);
      partition_gain_sizes[p]++;
    }
  }

  INT N_MOVES_MAX = 0;
  for (INT p = 0; p < k; p++) N_MOVES_MAX += partition_gain_sizes[p];

  bool *to_balance  = (bool*)calloc(k, sizeof(bool));  MEM_ERROR(to_balance);
  bool *locked_part = (bool*)calloc(k, sizeof(bool));  MEM_ERROR(locked_part);
  for (INT p = 0; p < k; p++) {
    for (INT wi = 0; wi < nwv; wi++) {
      if (part_sizes[p] >= capacity_cst[wi]) {
        locked_part[p] = true;
        to_balance[p]  = true;
      }
    }
  }

  PART * partitionp      = (PART*)malloc(nv*sizeof(PART)); MEM_ERROR(partitionp);
  memcpy(partitionp, partition, sizeof(PART)*nv);
  INT  best_cut_cost   = current_cut;
  INT  base_cur_cut    = current_cut;
  INT  N_MOVES = (INT)ceil((float)nv * (1.0f / (float)perform));

  for (INT i = 0; i < N_MOVES && i < N_MOVES_MAX; i++) {

    INT selected_part = -1;
    for (INT p = 0; p < k; p++)
      if (!to_balance[p] && partition_gain_sizes[p] > 0) selected_part = p;
    if (selected_part == -1) break;

    INT u = pqueue_dequeue(partition_moves[selected_part], partition_gains[selected_part], partition_gain_sizes[selected_part]);
    partition_gain_sizes[selected_part]--;

    if (partition[u] == selected_part) continue;

    is_locked[u]       = true;
    INT gain_u         = partition_gains[selected_part][u];
    current_cut        = base_cur_cut - gain_u;
    INT p_old          = partition[u];

    /* Apply move: update cnt[] for all arcs incident to u */
    {
      const INT *je = vth_begin(vth, u);
      const INT *jend = vth_end(vth, u);
      for (; je != jend; je++) {
        INT j  = *je;
        INT wj = hyperedges[idx_hyperedges[j]];
        INT nb_before = 0;
        for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp]>0) nb_before++;
        cnt[j*k + p_old]--;
        cnt[j*k + selected_part]++;
        INT nb_after = 0;
        for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp]>0) nb_after++;
        lambda[j] = ((nb_after)-1) * hyperedges[idx_hyperedges[j]];  /* keep lambda consistent */
        (void)wj; (void)nb_before;
      }
    }

    neighbours_parts[u*k + p_old]--;
    neighbours_parts[u*k + selected_part]++;
    partition[u] = selected_part;

    for (INT wi = 0; wi < nwv; wi++) {
      part_sizes[selected_part*nwv+wi] += h->ti_weights[u*nwv+wi];
      part_sizes[p_old*nwv+wi]        -= h->ti_weights[u*nwv+wi];
      if (part_sizes[selected_part*nwv+wi] >= capacity_cst[wi]) to_balance[selected_part] = true;
    }
    if (to_balance[p_old]) to_balance[p_old] = false;

    /* Update gains of unlocked neighbours */
    List *cell = out_neighbors[u];
    for (INT o = 0; o < out_neighbors[u]->size; o++) {
      INT v = cell->i;
      if (!is_locked[v]) {
        INT p_v = partition[v];
        for (INT p = 0; p < k; p++) {
          if (p == p_v) continue;
          if (neighbours_parts[v*k+p] == 0) continue;
          INT delta = 0;
          const INT *je = vth_begin(vth, v);
          const INT *jend = vth_end(vth, v);
          for (; je != jend; je++) {
            INT j   = *je;
            INT wj  = hyperedges[idx_hyperedges[j]];
            INT nb_before = 0;
            for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp]>0) nb_before++;
            cnt[j*k + p_v]--;
            cnt[j*k + p]++;
            INT nb_after = 0;
            for (INT pp = 0; pp < k; pp++) if (cnt[j*k+pp]>0) nb_after++;
            delta += (nb_after - nb_before) * wj;
            cnt[j*k + p_v]++;
            cnt[j*k + p]--;
          }
          INT new_gain = -delta;
          partition_gains[p][v] = new_gain;
          pqueue_sift_up(partition_moves[p], partition_gains[p], partition_gain_sizes[p]-1, partition_gain_sizes[p]);
        }
      }
      cell = cell->next;
    }

    if (current_cut < best_cut_cost) {
      best_cut_cost = current_cut;
      memcpy(partitionp, partition, sizeof(PART)*nv);
    }
  }

  memcpy(partition, partitionp, sizeof(PART)*nv);

  /* free section */
  free(cnt);
  free(is_in_halo);
  free(is_locked);
  free(partitionp);
  free(neighbours_parts);
  for (INT p = 0; p < k; p++) { free(partition_moves[p]); free(partition_gains[p]); }
  free(partition_moves);
  free(partition_gains);
  free(partition_gain_sizes);
  free(to_balance);
  free(is_in_part);
  free(lambda);
  free(locked_part);
  free(part_sizes);
  free(capacity_cst);

  return 0;
}




