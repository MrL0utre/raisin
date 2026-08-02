/* SPDX-License-Identifier: GPL-3.0-only */

#ifndef RAISIN_RNG_H
#define RAISIN_RNG_H

#include <stdint.h>

void     raisin_rng_seed(uint32_t seed);
uint32_t raisin_rng_next(void);
int      raisin_rng_bounded(int upper_bound);

#endif
