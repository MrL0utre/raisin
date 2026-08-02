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
/**   NAME       : ipart.h                                 **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   :  O(1) bucket queue for FM-style         **/
/**                 partitioners.                          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


/*
 * @file bqueue.h
 * @brief O(1) bucket queue for FM-style partitioners.
 *
 * Gains in FM are bounded integers in [min_gain, max_gain].  A bucket
 * queue exploits this by keeping one doubly-linked list per gain level,
 * giving true O(1) insert / delete / update-key instead of O(log n)
 * for a binary heap.  The "max nonempty" cursor scans downward at most
 * O(range) times across the entire lifetime of one FM pass, so
 * delete-max is O(1) amortised.
 *
 * Memory layout
 * -------------
 *   head[g]   – index of first vertex with gain g, or BQ_NONE
 *   nxt[v]    – next vertex in the same gain bucket
 *   prv[v]    – previous vertex in the same gain bucket (for O(1) delete)
 *   gval[v]   – current gain of vertex v
 *
 * All arrays are pre-allocated by the caller via bq_init() to avoid
 * repeated malloc/free inside FM loops.
 */

#ifndef BQUEUE_H
#define BQUEUE_H

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "commons.h"

#define BQ_NONE  (-1)

typedef struct {
    INT *head;       /* head[g + offset] = first vertex at gain g  */
    INT *nxt;        /* nxt[v]  = next vertex in bucket            */
    INT *prv;        /* prv[v]  = prev vertex in bucket            */
    INT *gval;       /* gval[v] = current gain of v                */
    INT  offset;     /* shift so that head[0] = min_gain           */
    INT  range;      /* total number of buckets = max-min+1        */
    INT  top;        /* index of highest non-empty bucket          */
    INT  size;       /* number of elements currently in queue      */
} BQueue;

/* ------------------------------------------------------------------ */
/* Lifecycle                                                           */
/* ------------------------------------------------------------------ */

/*
 * bq_alloc – allocate all internal arrays.
 *   nv        : max number of vertices that will ever be inserted
 *   min_gain  : lower bound on any gain value
 *   max_gain  : upper bound on any gain value
 */
static inline BQueue *
bq_alloc(INT nv, INT min_gain, INT max_gain)
{
    BQueue *bq = (BQueue *)malloc(sizeof(BQueue));
    if (!bq) return NULL;

    bq->offset = -min_gain;
    bq->range  = max_gain - min_gain + 1;
    bq->top    = -1;
    bq->size   = 0;

    bq->head = (INT *)malloc(sizeof(INT) * bq->range);
    bq->nxt  = (INT *)malloc(sizeof(INT) * nv);
    bq->prv  = (INT *)malloc(sizeof(INT) * nv);
    bq->gval = (INT *)malloc(sizeof(INT) * nv);

    if (!bq->head || !bq->nxt || !bq->prv || !bq->gval) {
        free(bq->head); free(bq->nxt); free(bq->prv); free(bq->gval);
        free(bq); return NULL;
    }

    for (INT i = 0; i < bq->range; i++) bq->head[i] = BQ_NONE;
    for (INT v = 0; v < nv; v++) { bq->nxt[v] = BQ_NONE; bq->prv[v] = BQ_NONE; bq->gval[v] = min_gain - 1; }

    return bq;
}

/* bq_reset – O(range + nv): clear all buckets, reset gains */
static inline void
bq_reset(BQueue *bq, INT nv, INT min_gain)
{
    for (INT i = 0; i < bq->range; i++) bq->head[i] = BQ_NONE;
    for (INT v = 0; v < nv; v++) {
        bq->nxt[v] = BQ_NONE;
        bq->prv[v] = BQ_NONE;
        bq->gval[v] = min_gain - 1;
    }
    bq->top  = -1;
    bq->size = 0;
}

static inline void
bq_free(BQueue *bq)
{
    if (!bq) return;
    free(bq->head);
    free(bq->nxt);
    free(bq->prv);
    free(bq->gval);
    free(bq);
}

/* ------------------------------------------------------------------ */
/* Core operations – all O(1)                                         */
/* ------------------------------------------------------------------ */

/* Insert vertex v with gain g.  v must NOT already be in the queue. */
static inline void
bq_insert(BQueue *bq, INT v, INT g)
{
    INT b = g + bq->offset;
    bq->gval[v]  = g;
    bq->prv[v]   = BQ_NONE;
    bq->nxt[v]   = bq->head[b];
    if (bq->head[b] != BQ_NONE)
        bq->prv[bq->head[b]] = v;
    bq->head[b]  = v;
    if (b > bq->top) bq->top = b;
    bq->size++;
}

/* Remove vertex v from its current bucket.  v must be in the queue. */
static inline void
bq_remove(BQueue *bq, INT v)
{
    INT b = bq->gval[v] + bq->offset;
    if (bq->prv[v] != BQ_NONE)
        bq->nxt[bq->prv[v]] = bq->nxt[v];
    else
        bq->head[b] = bq->nxt[v];

    if (bq->nxt[v] != BQ_NONE)
        bq->prv[bq->nxt[v]] = bq->prv[v];

    /* Lazily update top: the dequeue-max loop will find the real top. */
    bq->nxt[v] = BQ_NONE;
    bq->prv[v] = BQ_NONE;
    bq->size--;
}

/* Update gain of v from old value to new_g.  O(1). */
static inline void
bq_update(BQueue *bq, INT v, INT new_g)
{
    bq_remove(bq, v);
    bq_insert(bq, v, new_g);
}

/*
 * Dequeue the vertex with maximum gain.  O(1) amortised.
 * Returns BQ_NONE if queue is empty.
 */
static inline INT
bq_dequeue_max(BQueue *bq)
{
    if (bq->size == 0) return BQ_NONE;

    /* Scan down to find the true top (lazy deletion of empty buckets). */
    while (bq->top >= 0 && bq->head[bq->top] == BQ_NONE)
        bq->top--;

    if (bq->top < 0) return BQ_NONE;

    INT v = bq->head[bq->top];
    bq_remove(bq, v);
    return v;
}

/* Peek at max gain without removing. */
static inline INT
bq_max_gain(BQueue *bq)
{
    while (bq->top >= 0 && bq->head[bq->top] == BQ_NONE)
        bq->top--;
    if (bq->top < 0) return INT_MIN;
    return bq->top - bq->offset;
}

static inline INT bq_empty(const BQueue *bq) { return bq->size == 0; }

#endif /* BQUEUE_H */
