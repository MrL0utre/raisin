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
 * @file list.c
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 13 oct 2023
 * @brief  These lines are list functions definitions. 
 */
#include "list.h"

/**
 * @brief Function initializing a list of size 0 and 
 * value -1. 
 *
 * @param l                List.
 *
 * @return Void.
 */
void 
new_list(List * l) 
{

    l->i    = -1;
    l->size = 0;
    l->next = NULL;
    return;
}

/**
 * @brief The function add_list concatenate the lists la
 * and lb. The list lb was pushed after la. 
 *
 * @param la                List.
 * @param la                List.
 *
 * @return Void.
 */
void
add_list(List * la, 
         List * lb)
{
    if (la != NULL && lb != NULL) 
      {
        la->next = lb;
        la->size += lb->size;
      }
    return;
}

/**
 * @brief This function add an element (i) at the 
 * end of the list la and increase the size of la
 * plus one. 
 *
 * @param la                List.
 * @param i                 Integer.
 *
 * @return Void.
 */
void 
list_add_element(List * la, 
                 INT    i)
{
    if (la != NULL) 
      {

        if (la->size > 0) 
          {
            List * lb = (List*)malloc(sizeof(List));
            MEM_ERROR(lb);

            new_list(lb);

            lb->i = i;
            lb->size = 1;
            List * cell = la;
            
            while(cell->next != NULL) 
              {
                cell = cell->next;
              }
            
            cell->next = lb;
            la->size += 1;
          } 
        else 
          {
            la->i = i;
            la->size += 1;
          }
      }
    return;
}

/**
 * @brief This function delete the next element of 
 * a list and decrease its size. 
 *
 * @param la                List.
 * @param i                 Integer.
 *
 * @return Void.
 */
void 
list_del_next_element(List *la)
{
    List * lb = la->next;
    if (la != NULL && lb != NULL) 
      {
        if (lb->next != NULL) 
          {
            List * lc = lb->next;
            la->next = lc;
            la->size -= 1;
            free(lb);
          } 
        else 
          {
            free(lb);
            la->next = NULL;
            la->size -= 1;
          }
      }
    return;
}

/**
 * @brief This function delete a list. 
 *
 * @param l                List.
 *
 * @return Void.
 */
void 
delete_list(List *l) 
{
    List * tmp;
    while (l != NULL) 
      {
        tmp = l->next;
        free(l);
        l   = tmp;
      }
    return;
}
