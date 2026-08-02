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
