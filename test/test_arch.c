/* SPDX-License-Identifier: GPL-3.0-only */

#include "a.h"

#include <assert.h>
#include <stdlib.h>

static Arch *
new_arch(void)
{
    Arch *arch = (Arch *)calloc(1, sizeof(Arch));
    assert(arch != NULL);
    return arch;
}

int
main(int argc, char **argv)
{
    Arch *arch;

    assert(argc == 4);

    arch = new_arch();
    assert(arch_load(arch, argv[1], false) == 0);
    assert(arch->i_n == 3);
    assert(arch_has_link(arch, 0, 1));
    assert(arch->ti_link_delay[0 * 3 + 1] == 0);
    assert(!arch_has_link(arch, 0, 2));
    assert(arch->ti_delay[0 * 3 + 2] == 4);
    assert(arch_next_hop(arch, 0, 2) == 1);
    assert(arch_next_hop(arch, 2, 0) == 1);
    assert(arch->ti_capacity[0 * 3 + 2] == 0);
    arch_free(arch);

    arch = new_arch();
    assert(arch_load(arch, argv[2], false) == 7);
    assert(arch->ti_delay == NULL);
    arch_free(arch);

    arch = new_arch();
    assert(arch_load(arch, argv[3], false) == 6);
    assert(arch->ti_delay == NULL);
    arch_free(arch);

    return EXIT_SUCCESS;
}
