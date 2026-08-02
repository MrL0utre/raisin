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


#include "vector.h"

void
new_vector(Vector *v, INT size)
{
    v->v    = (INT *)malloc(sizeof(INT) * size);
    MEM_ERROR(v->v);
    v->size = size;
    v->i    = 0;
}

void
push_back(Vector *v, INT val)
{
    if (v->i == v->size) {
        INT new_size = v->size * 2;
        void *tmp    = realloc(v->v, sizeof(INT) * new_size);
        /* Fix: check realloc result and update the pointer */
        if (tmp == NULL) {
            printf("vector.c: unable to realloc buffer\n");
            exit(1);
        }
        v->v    = (INT *)tmp;
        v->size = new_size;
    }
    v->v[v->i++] = val;
}

/*
 * delete_vector_data() frees only the inner data buffer.
 * Use this for Vectors embedded in a Matrix
 */
void
delete_vector_data(Vector *v)
{
    if (v && v->v) {
        free(v->v);
        v->v    = NULL;
        v->size = 0;
        v->i    = 0;
    }
}

/*
 * delete_vector() frees both the data buffer and the Vector struct itself.
 * Only call this when the Vector was individually heap-allocated via malloc.
 */
void
delete_vector(Vector *v)
{
    if (!v) return;
    free(v->v);
    free(v);
}
