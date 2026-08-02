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

typedef struct resource_usage {
  INT *ti_part_load;
  INT  i_parts;
  INT  i_resources;
  INT  i_overloaded_dimensions;
  INT  i_overloaded_parts;
  INT  i_max_overload;
} ResourceUsage;

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

int  resource_usage_compute(ResourceUsage *usage,
                            const Hypergraph *hypergraph,
                            const Arch *arch,
                            const PART *partition,
                            INT part_count);
void resource_usage_free(ResourceUsage *usage);
bool resource_usage_is_feasible(const ResourceUsage *usage);
INT  resource_usage_part_load(const ResourceUsage *usage,
                              INT part,
                              INT resource);
int  resource_partition_repair(const Hypergraph *hypergraph,
                               const Arch *arch,
                               PART *partition,
                               INT part_count,
                               INT *move_count);
int  communication_partition_repair(const Hypergraph *hypergraph,
                                    const Arch *arch,
                                    PART *partition,
                                    INT part_count,
                                    INT *move_count);

#endif
