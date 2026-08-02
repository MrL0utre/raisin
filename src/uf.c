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
/**   NAME       : uf.c                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are union-find              **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include "uf.h"

Union_find *
new_set(INT element)
{
    Union_find *set = (Union_find *)malloc(sizeof(Union_find));
    MEM_ERROR(set);
    set->value  = element;
    set->parent = NULL;
    set->depth  = 1;
    return set;
}

/*
 * find() with full path compression (two-pass).
 *
 * Pass 1: walk to the root.
 * Pass 2: point every node on the path directly at the root.
 *
 */
Union_find *
find(Union_find *element)
{
    if (!element)
        return NULL;

    /* Pass 1 - locate root */
    Union_find *root = element;
    while (root->parent != NULL)
        root = root->parent;

    /* Pass 2 - flatten all intermediate pointers */
    while (element->parent != NULL) {
        Union_find *next = element->parent;
        element->parent  = root;
        element          = next;
    }

    return root;
}

Union_find *
merge(Union_find *a,
      Union_find *b)
{
    Union_find *root_a   = find(a);
    Union_find *root_b   = find(b);
    Union_find *new_root = root_a;

    if (root_a == root_b)
        return root_a;

    if (root_a->depth >= root_b->depth) {
        root_b->parent = root_a;
        if (root_a->depth == root_b->depth)
            root_a->depth++;
    } else {
        root_a->parent = root_b;
        root_b->depth++;
        new_root = root_b;
    }

    return new_root;
}
