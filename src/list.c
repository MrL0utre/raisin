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
/**   FUNCTION   : These lines are list                    **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/


#include "list.h"

void new_list(List * l) {

    l->i    = -1;
    l->size = 0;
    l->next = NULL;
    return;
}

void add_list(List * la, List * lb){
    if (la != NULL && lb != NULL) {
        la->next = lb;
        la->size += lb->size;
    }
    return;
}

void list_add_element(List * la, INT i){
    if (la != NULL) {

        if (la->size>0) {
            List * lb = (List*)malloc(sizeof(List));
            MEM_ERROR(lb);
            new_list(lb);
            lb->i = i;
            lb->size=1;
            List * cell = la;
            while(cell->next != NULL) {
                cell = cell->next;
            }
            cell->next = lb;
            la->size +=1;
        } else {
            la->i = i;
            la->size+=1;
        }
    }
    return;
}

void list_del_next_element(List * la){
    List * lb = la->next;
    if (la != NULL && lb != NULL) {

        if (lb->next != NULL) {
            List * lc = lb->next;
            la->next = lc;
            free(lb);
        } else {
            free(lb);
            la->next = NULL;
        }
    }
    return;
}

void delete_list(List * l) {

    List * tmp;
    while (l != NULL) {

        tmp = l->next;
        free(l);
        l   = tmp;
    }
    return;
}
