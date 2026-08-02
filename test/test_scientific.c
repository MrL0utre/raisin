/* SPDX-License-Identifier: GPL-3.0-only */

#include "capacity.h"
#include "fm.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

static List **
new_lists(INT count)
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
    Hypergraph hypergraph;
    Arch *arch;
    List **out_neighbors;
    List **in_neighbors;
    INT sort[3];
    INT lambda[2];
    PART local[] = {0, 0, 0};
    PART routed[] = {0, 2, 0};
    CommunicationUsage usage;

    assert(argc == 3);
    memset(&hypergraph, 0, sizeof(hypergraph));
    assert(rbh_load_base(&hypergraph, argv[1], 0, false) == 0);
    assert(rbh_validate(&hypergraph) == 0);

    arch = (Arch *)calloc(1, sizeof(Arch));
    assert(arch != NULL);
    assert(arch_load(arch, argv[2], false) == 0);

    out_neighbors = new_lists(hypergraph.i_vertices);
    in_neighbors = new_lists(hypergraph.i_vertices);
    assert(compute_list_neighbors(&hypergraph, out_neighbors) == 0);
    assert(compute_list_in_neighbors(&hypergraph, out_neighbors,
                                     in_neighbors) == 0);
    assert(topological_sort(&hypergraph, out_neighbors, in_neighbors,
                            sort) == 0);

    assert(compute_partition_criticality(
               &hypergraph, arch, out_neighbors, in_neighbors,
               sort, local) == 35);
    assert(compute_partition_cut(
               &hypergraph, arch, out_neighbors, in_neighbors,
               sort, local, lambda, 3) == 0);

    assert(compute_partition_criticality(
               &hypergraph, arch, out_neighbors, in_neighbors,
               sort, routed) == 43);
    assert(compute_partition_cut(
               &hypergraph, arch, out_neighbors, in_neighbors,
               sort, routed, lambda, 3) == 2);
    assert(lambda[0] == 1);
    assert(lambda[1] == 1);

    assert(communication_usage_compute(
               &usage, &hypergraph, arch, routed, 3) == 0);
    assert(communication_usage_link_load(&usage, 0, 1) == 2);
    assert(communication_usage_link_load(&usage, 1, 2) == 2);
    assert(usage.i_routed_signal_hops == 4);
    assert(usage.i_overloaded_links == 2);
    communication_usage_free(&usage);

    free_lists(in_neighbors, hypergraph.i_vertices);
    free_lists(out_neighbors, hypergraph.i_vertices);
    arch_free(arch);
    rbh_free(&hypergraph);
    return EXIT_SUCCESS;
}
