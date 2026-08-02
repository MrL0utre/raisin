/* SPDX-License-Identifier: GPL-3.0-only */

#include "capacity.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

int
main(int argc, char **argv)
{
    Arch *arch;
    Hypergraph hypergraph;
    CommunicationUsage usage;
    PART partition[] = {0, 2, 2, 0, 2};
    INT indices[] = {0, 5};
    INT hyperedges[] = {
        1, 3, 0, 1, 2,
        1, 2, 3, 4
    };

    assert(argc == 2);
    arch = (Arch *)calloc(1, sizeof(Arch));
    assert(arch != NULL);
    assert(arch_load(arch, argv[1], false) == 0);

    memset(&hypergraph, 0, sizeof(hypergraph));
    hypergraph.i_vertices = 5;
    hypergraph.i_hyperedges = 2;
    hypergraph.ti_idx_hyperedges = indices;
    hypergraph.ti_hyperedges = hyperedges;

    assert(communication_usage_compute(&usage, &hypergraph, arch,
                                       partition, 3) == 0);
    assert(communication_usage_link_load(&usage, 0, 1) == 2);
    assert(communication_usage_link_load(&usage, 1, 2) == 2);
    assert(communication_usage_link_load(&usage, 0, 2) == 0);
    assert(usage.i_routed_signal_hops == 4);
    assert(usage.i_overloaded_links == 2);
    assert(usage.i_max_overload == 1);
    assert(!communication_usage_is_feasible(&usage));

    communication_usage_free(&usage);
    arch_free(arch);
    return EXIT_SUCCESS;
}
