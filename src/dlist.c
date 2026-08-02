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
/**   NAME       : list.c                                  **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are double linked list      **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                 to   : xx xxx xxxx     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#include "dlist.h"

void new_dlist(Dlist * l) {

    l->i    = -1;
    l->size = 0;
    l->next = NULL;
    l->pred = NULL;
    return;
}

void add_dlist(Dlist * la, Dlist * lb){
    if (la != NULL && lb != NULL) {
        la->next = lb;
        la->size += lb->size;
        lb->pred = la;
    }
    return;
}

void dlist_add_element(Dlist * la, INT i){
    if (la != NULL) {

        if (la->size>0) {
            Dlist * lb = (Dlist*)malloc(sizeof(Dlist));
            MEM_ERROR(lb);
            new_dlist(lb);
            lb->i = i;
            lb->size=1;
            Dlist * cell = la;
            while(cell->next != NULL) {
                cell = cell->next;
            }
            cell->next = lb;
            la->size +=1;
            lb->pred = cell;
        } else {
            la->i = i;
            la->size+=1;
        }
    }
    return;
}

void delete_dlist(Dlist * l) {

    Dlist * tmp;
    while (l != NULL) {

        tmp = l->next;
        free(l);
        l   = tmp;
    }
    return;
}
