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
