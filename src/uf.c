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
/**   NAME       : uf.c                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are union-find              **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "uf.h"

Union_find * new_set(INT element) {
    Union_find * set = (Union_find*)malloc(sizeof(Union_find));
    MEM_ERROR(set);
    set->value      = element;
    set->parent     = NULL;
    set->depth      = 1;
    return set;
}

Union_find * find    (Union_find * element) {
    if(!element) {
        return NULL;
    }
    Union_find * base_element = element;
    while(element->parent != NULL) {
        element = element->parent;
    }
    if (element != base_element->parent && element != base_element) {
        base_element->parent = element;
    }
    return element;
}

Union_find * merge   (Union_find * a, Union_find * b) {

    Union_find * root_a   = find(a);
    Union_find * root_b   = find(b);
    Union_find * new_root = root_a;

    if(root_a != root_b){
        if(root_a->depth > root_b->depth){
            root_b->parent = root_a;
            root_a->depth++;
        }else{
            root_a->parent = root_b;
            root_b->depth++;
            new_root = root_b;
        }
    }
    return new_root;
}