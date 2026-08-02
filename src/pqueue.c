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
/**   NAME       : pqueue.h                                **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are priority queue          **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "pqueue.h"

/* -------------------------------------------------------------------------- */
/* Internal helpers                                                            */
/* -------------------------------------------------------------------------- */

static inline void
_swap(INT *queue, INT i, INT j)
{
    INT tmp   = queue[i];
    queue[i]  = queue[j];
    queue[j]  = tmp;
}

/* -------------------------------------------------------------------------- */
/* Public API                                                                  */
/* -------------------------------------------------------------------------- */

/**
 * pqueue_sift_up – iterative (was recursive).
 * Restores the heap property upward from position i.
 */
void
pqueue_sift_up(INT *queue, INT *queue_w, INT i, INT n)
{
    (void)n; /* n unused in sift_up but kept for API compatibility */
    while (i > 0) {
        INT parent = (i - 1) / 2;
        if (queue_w[queue[i]] > queue_w[queue[parent]]) {
            _swap(queue, i, parent);
            i = parent;
        } else {
            break;
        }
    }
}

/**
 * pqueue_sift_down – iterative (was recursive).
 * Restores the heap property downward from position i.
 */
void
pqueue_sift_down(INT *queue, INT *queue_w, INT i, INT n)
{
    while (1) {
        INT left  = 2 * i + 1;
        INT right = 2 * i + 2;
        INT largest = i;

        if (left  < n && queue_w[queue[left]]  > queue_w[queue[largest]])
            largest = left;
        if (right < n && queue_w[queue[right]] > queue_w[queue[largest]])
            largest = right;

        if (largest == i)
            break;

        _swap(queue, i, largest);
        i = largest;
    }
}

/**
 * pqueue_heapify – O(n) bottom-up heap construction (was empty stub).
 */
void
pqueue_heapify(INT *queue, INT *queue_w, INT n)
{
    /* Start from the last internal node and sift down each one */
    for (INT i = n / 2 - 1; i >= 0; i--)
        pqueue_sift_down(queue, queue_w, i, n);
}

/**
 * pqueue_dequeue – remove and return the max-priority element.
 *
 */
INT
pqueue_dequeue(INT *queue, INT *queue_w, INT n)
{
    INT top    = queue[0];
    queue[0]   = queue[n - 1];  /* move last element to root */
    /* caller must do: size-- */
    pqueue_sift_down(queue, queue_w, 0, n - 1);
    return top;
}

/**
 * pqueue_add_element – append element and sift up.
 */
void
pqueue_add_element(INT *queue, INT *queue_w, INT i, INT n)
{
    queue[n] = i;
    pqueue_sift_up(queue, queue_w, n, n + 1);
    /* caller must do: size++ */
}
