/* SPDX-License-Identifier: GPL-3.0-only */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "rng.h"

static uint32_t raisin_rng_state = 1U;

void
raisin_rng_seed(uint32_t seed)
{
    raisin_rng_state = seed;
}

uint32_t
raisin_rng_next(void)
{
    raisin_rng_state = raisin_rng_state * UINT32_C(214013) +
                       UINT32_C(2531011);
    return (raisin_rng_state >> 16U) & UINT32_C(0x7fff);
}

int
raisin_rng_bounded(int upper_bound)
{
    if (upper_bound <= 0) {
        fprintf(stderr, "raisin_rng_bounded: upper bound must be positive\n");
        abort();
    }

    return (int)(raisin_rng_next() % (uint32_t)upper_bound);
}
