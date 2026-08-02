#include "rbh.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static List **
allocate_lists(INT count)
{
    List **lists = (List **)calloc((size_t)count, sizeof(List *));
    assert(lists != NULL);
    for (INT i = 0; i < count; i++) {
        lists[i] = (List *)malloc(sizeof(List));
        assert(lists[i] != NULL);
        new_list(lists[i]);
    }
    return lists;
}

static void
free_lists(List **lists, INT count)
{
    for (INT i = 0; i < count; i++)
        delete_list(lists[i]);
    free(lists);
}

int
main(int argc, char **argv)
{
    Hypergraph h;
    List **out_neighbors;
    List **in_neighbors;
    INT *sort;
    bool *red_mask;
    INT red_count = 0;

    if (argc != 2) {
        fprintf(stderr, "usage: %s HYPERGRAPH\n", argv[0]);
        return EXIT_FAILURE;
    }

    memset(&h, 0, sizeof(h));
    assert(rbh_load_base(&h, argv[1], 0, false) == 0);
    assert(rbh_validate(&h) == 0);
    assert(h.is_red != NULL);
    assert(h.vth != NULL);

    red_mask = (bool *)malloc((size_t)h.i_vertices * sizeof(bool));
    assert(red_mask != NULL);
    memcpy(red_mask, h.is_red, (size_t)h.i_vertices * sizeof(bool));

    for (INT v = 0; v < h.i_vertices; v++)
        red_count += h.is_red[v] ? 1 : 0;
    assert(red_count == h.i_reds);

    out_neighbors = allocate_lists(h.i_vertices);
    in_neighbors = allocate_lists(h.i_vertices);
    sort = (INT *)malloc((size_t)h.i_vertices * sizeof(INT));
    assert(sort != NULL);

    assert(compute_list_neighbors(&h, out_neighbors) == 0);
    assert(compute_list_in_neighbors(&h, out_neighbors, in_neighbors) == 0);
    assert(topological_sort(&h, out_neighbors, in_neighbors, sort) == 0);
    assert(compute_criticality(&h, out_neighbors, in_neighbors, sort) == 0);
    assert(memcmp(red_mask, h.is_red,
                  (size_t)h.i_vertices * sizeof(bool)) == 0);

    for (INT v = 0; v < h.i_vertices; v++) {
        assert(h.ti_criticalities_left[v] > 0);
        assert(h.ti_criticalities_right[v] > 0);
    }

    free(sort);
    free_lists(in_neighbors, h.i_vertices);
    free_lists(out_neighbors, h.i_vertices);
    free(red_mask);
    rbh_free(&h);
    return EXIT_SUCCESS;
}
