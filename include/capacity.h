/* SPDX-License-Identifier: GPL-3.0-only */

#ifndef RAISIN_CAPACITY_H
#define RAISIN_CAPACITY_H

#include "a.h"
#include "rbh.h"

typedef struct communication_usage {
  INT *ti_link_load;
  INT  i_parts;
  INT  i_routed_signal_hops;
  INT  i_overloaded_links;
  INT  i_max_overload;
} CommunicationUsage;

int  communication_usage_compute(CommunicationUsage *usage,
                                 const Hypergraph *hypergraph,
                                 const Arch *arch,
                                 const PART *partition,
                                 INT part_count);
void communication_usage_free(CommunicationUsage *usage);
bool communication_usage_is_feasible(const CommunicationUsage *usage);
INT  communication_usage_link_load(const CommunicationUsage *usage,
                                   INT u,
                                   INT v);

#endif
