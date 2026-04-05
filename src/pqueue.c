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
 * @file pqueue.c
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  These lines are priority queue    
 *          functions definitions.
 */

#include "pqueue.h"

/**
 * @brief Function organizing the heap implementing the
 * priority queue. This function is utilized when an element
 * is added to the priority queue.    
 *
 * @param queue         Priority queue.
 * @param queue_w       Array of weight (priority).  
 * @param i             Integer (node).
 * @param n             Integer (size)
 *
 * @return Void.
 */
void 
pqueue_sift_up(INT * queue, 
              INT * queue_w, 
              INT   i, 
              INT   n)
{    
    INT i_father = (i - 1) / 2;
    INT tmp;

    if(queue_w[queue[i]] > queue_w[queue[i_father]])
      {
        tmp             = queue[i];
        queue[i]        = queue[i_father];
        queue[i_father] = tmp;

        if(i_father != 0)
            pqueue_sift_up(queue, queue_w, i_father, n);
      }
}

/**
 * @brief Function   
 *
 * @param queue         Priority queue.
 * @param queue_w       Array of weight (priority).  
 * @param i             Integer (node).
 * @param n             Integer (size)
 *
 * @return Void.
 */
void 
pqueue_sift_down(INT * queue, 
                INT * queue_w, 
                INT   i, 
                INT   n)
{    
    INT index_child_left  = i * 2 + 1;
    INT index_child_right = i * 2 + 2;
    INT m = queue_w[queue[i]];
    
    if(index_child_right < n)
      {
        m = MAX3(queue_w[queue[i]], queue_w[queue[index_child_left]], queue_w[queue[index_child_right]]);
      } 
    else if(index_child_left < n)
      {
        m = MAX(queue_w[queue[i]], queue_w[queue[index_child_left]]);
      }

    INT i_max;

    if(m > queue_w[queue[i]]) 
      {
        i_max = index_child_left;
        if(index_child_right < n && queue_w[queue[i_max]] < queue_w[queue[index_child_right]])
          {
            i_max = index_child_right;
          }
        
        INT tmp      = queue[i];
        queue[i]     = queue[i_max];
        queue[i_max] = tmp;
        
        if(2 * i_max + 1 < n)
          {
            pqueue_sift_down(queue, queue_w, i_max, n);
          }
    }
}

/**
 * @brief Function   
 *
 * @param queue         Priority queue.
 * @param queue_w       Array of weight (priority).  
 * @param i             Integer (node).
 * @param n             Integer (size)
 *
 * @return Void.
 */
void 
pqueue_heapify(INT * queue, 
               INT * queue_w, 
               INT   n)
{
    ;
}

/**
 * @brief Function deleting the first element and returning it.
 *
 * @param queue         Priority queue.
 * @param queue_w       Array of weight (priority).  
 * @param n             Integer (size)
 *
 * @return Integer (the dequeued element).
 */
INT 
pqueue_dequeue(INT * queue, 
               INT * queue_w, 
               INT n){

    INT i = queue[0];
    queue[0] = queue[n-1];
    n--;
    
    pqueue_sift_down(queue, queue_w, 0, n);
    return i;
}

/**
 * @brief Function adding an element according to its 
 * priority by calling pqueue_sift_up function.  
 *
 * @param queue         Priority queue.
 * @param queue_w       Array of weight (priority).  
 * @param i             Integer (node).
 * @param n             Integer (size)
 *
 * @return Void.
 */
void 
pqueue_add_element(INT * queue, 
                   INT * queue_w, 
                   INT   i, 
                   INT   n)
{     
    queue[n++] = i;

    pqueue_sift_up(queue, queue_w, n - 1, n);
}

