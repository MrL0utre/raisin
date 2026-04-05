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
 * @file pqueue.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  These lines are priority queue 
 *         functions declarations.      
 */
#ifndef PQUEUE_H

#include <stdlib.h>
#include "commons.h"

#define PQUEUE_H

void 
pqueue_sift_up(INT * queue, 
               INT * queue_w, 
               INT   i, 
               INT   n);

void 
pqueue_sift_down(INT * queue, 
                 INT * queue_w, 
                 INT   i, 
                 INT   n);

void 
pqueue_heapify(INT * queue, 
               INT * queue_w, 
               INT   n);

INT 
pqueue_dequeue(INT * queue, 
               INT * queue_w, 
               INT   n);

void 
pqueue_add_element(INT * queue, 
                   INT * queue_w, 
                   INT   i, 
                   INT   n);

#endif
