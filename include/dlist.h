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
 * @file dlist.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  This file is part of the project.    
 *        It contains declarations of double linked list functions
 *        for red-black hypergraphs partitioning software (raisin).
 *      
 */
#ifndef DLIST_H

#include <stdlib.h>
#include "commons.h"

#define DLIST_H

struct dlist{
  INT    i;
  INT    size;
  struct dlist *  next;
  struct dlist *  pred;
};
typedef struct dlist Dlist;

void new_dlist(Dlist * l);

void add_dlist(Dlist * l, Dlist * lb);

void dlist_add_element(Dlist * la, INT i);

void delete_dlist(Dlist * l);

#endif
