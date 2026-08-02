/* SPDX-License-Identifier: GPL-3.0-only */

#include "capacity.h"

#include <limits.h>
#include <string.h>

static void
communication_usage_reset(CommunicationUsage *usage)
{
  usage->ti_link_load = NULL;
  usage->i_parts = 0;
  usage->i_routed_signal_hops = 0;
  usage->i_overloaded_links = 0;
  usage->i_max_overload = 0;
}

int
communication_usage_compute(CommunicationUsage *usage,
                            const Hypergraph *hypergraph,
                            const Arch *arch,
                            const PART *partition,
                            INT part_count)
{
  bool *used_links;
  size_t matrix_size;

  if (usage == NULL || hypergraph == NULL || arch == NULL ||
      partition == NULL || part_count <= 0 || part_count > arch->i_n ||
      arch->ti_next_hop == NULL || arch->ti_capacity == NULL)
    return 1;

  communication_usage_reset(usage);
  matrix_size = (size_t)arch->i_n * (size_t)arch->i_n;
  usage->ti_link_load = (INT *)calloc(matrix_size, sizeof(INT));
  MEM_ERROR(usage->ti_link_load);
  used_links = (bool *)calloc(matrix_size, sizeof(bool));
  MEM_ERROR(used_links);
  usage->i_parts = arch->i_n;

  for (INT vertex = 0; vertex < hypergraph->i_vertices; vertex++) {
    if ((INT)partition[vertex] >= part_count) {
      free(used_links);
      communication_usage_free(usage);
      return 2;
    }
  }

  for (INT edge = 0; edge < hypergraph->i_hyperedges; edge++) {
    INT base = hypergraph->ti_idx_hyperedges[edge];
    INT pin_count = hypergraph->ti_hyperedges[base + 1];
    INT source_vertex = hypergraph->ti_hyperedges[base + 2];
    INT source_part = partition[source_vertex];

    memset(used_links, 0, matrix_size * sizeof(bool));
    for (INT pin = 1; pin < pin_count; pin++) {
      INT sink_vertex = hypergraph->ti_hyperedges[base + 2 + pin];
      INT sink_part = partition[sink_vertex];
      INT current = source_part;
      INT hop_count = 0;

      while (current != sink_part) {
        INT next = arch_next_hop(arch, current, sink_part);
        INT low;
        INT high;
        size_t link;

        if (next < 0 || next >= arch->i_n || next == current ||
            !arch_has_link(arch, current, next) || hop_count++ >= arch->i_n) {
          free(used_links);
          communication_usage_free(usage);
          return 3;
        }
        low = MIN(current, next);
        high = MAX(current, next);
        link = (size_t)low * (size_t)arch->i_n + (size_t)high;
        used_links[link] = true;
        current = next;
      }
    }

    for (INT u = 0; u < arch->i_n; u++) {
      for (INT v = u + 1; v < arch->i_n; v++) {
        size_t link = (size_t)u * (size_t)arch->i_n + (size_t)v;
        if (!used_links[link])
          continue;
        if (usage->ti_link_load[link] == INT_MAX ||
            usage->i_routed_signal_hops == INT_MAX) {
          free(used_links);
          communication_usage_free(usage);
          return 4;
        }
        usage->ti_link_load[link]++;
        usage->ti_link_load[(size_t)v * (size_t)arch->i_n + (size_t)u]++;
        usage->i_routed_signal_hops++;
      }
    }
  }

  for (INT u = 0; u < arch->i_n; u++) {
    for (INT v = u + 1; v < arch->i_n; v++) {
      size_t link = (size_t)u * (size_t)arch->i_n + (size_t)v;
      if (!arch_has_link(arch, u, v))
        continue;
      INT overload = usage->ti_link_load[link] - arch->ti_capacity[link];
      if (overload > 0) {
        usage->i_overloaded_links++;
        usage->i_max_overload = MAX(usage->i_max_overload, overload);
      }
    }
  }

  free(used_links);
  return 0;
}

void
communication_usage_free(CommunicationUsage *usage)
{
  if (usage == NULL)
    return;
  free(usage->ti_link_load);
  communication_usage_reset(usage);
}

bool
communication_usage_is_feasible(const CommunicationUsage *usage)
{
  return usage != NULL && usage->ti_link_load != NULL &&
         usage->i_overloaded_links == 0;
}

INT
communication_usage_link_load(const CommunicationUsage *usage, INT u, INT v)
{
  if (usage == NULL || usage->ti_link_load == NULL ||
      u < 0 || v < 0 || u >= usage->i_parts || v >= usage->i_parts)
    return -1;
  return usage->ti_link_load[u * usage->i_parts + v];
}

static void
resource_usage_reset(ResourceUsage *usage)
{
  usage->ti_part_load = NULL;
  usage->i_parts = 0;
  usage->i_resources = 0;
  usage->i_overloaded_dimensions = 0;
  usage->i_overloaded_parts = 0;
  usage->i_max_overload = 0;
}

int
resource_usage_compute(ResourceUsage *usage,
                       const Hypergraph *hypergraph,
                       const Arch *arch,
                       const PART *partition,
                       INT part_count)
{
  size_t load_count;

  if (usage == NULL || hypergraph == NULL || arch == NULL ||
      partition == NULL || part_count <= 0 || part_count > arch->i_n ||
      !arch_has_part_capacities(arch, hypergraph->i_weights))
    return 1;

  resource_usage_reset(usage);
  load_count = (size_t)arch->i_n * (size_t)arch->i_resources;
  usage->ti_part_load = (INT *)calloc(load_count, sizeof(INT));
  MEM_ERROR(usage->ti_part_load);
  usage->i_parts = arch->i_n;
  usage->i_resources = arch->i_resources;

  for (INT vertex = 0; vertex < hypergraph->i_vertices; vertex++) {
    INT part = partition[vertex];
    if (part >= part_count) {
      resource_usage_free(usage);
      return 2;
    }
    for (INT resource = 0; resource < hypergraph->i_weights; resource++) {
      size_t index = (size_t)part * (size_t)usage->i_resources +
                     (size_t)resource;
      INT weight = hypergraph->ti_weights[
          vertex * hypergraph->i_weights + resource];
      if (weight > INT_MAX - usage->ti_part_load[index]) {
        resource_usage_free(usage);
        return 3;
      }
      usage->ti_part_load[index] += weight;
    }
  }

  for (INT part = 0; part < part_count; part++) {
    bool part_overloaded = false;
    for (INT resource = 0; resource < usage->i_resources; resource++) {
      size_t index = (size_t)part * (size_t)usage->i_resources +
                     (size_t)resource;
      INT overload = usage->ti_part_load[index] -
                     arch_part_capacity(arch, part, resource);
      if (overload > 0) {
        usage->i_overloaded_dimensions++;
        usage->i_max_overload = MAX(usage->i_max_overload, overload);
        part_overloaded = true;
      }
    }
    if (part_overloaded)
      usage->i_overloaded_parts++;
  }
  return 0;
}

void
resource_usage_free(ResourceUsage *usage)
{
  if (usage == NULL)
    return;
  free(usage->ti_part_load);
  resource_usage_reset(usage);
}

bool
resource_usage_is_feasible(const ResourceUsage *usage)
{
  return usage != NULL && usage->ti_part_load != NULL &&
         usage->i_overloaded_dimensions == 0;
}

INT
resource_usage_part_load(const ResourceUsage *usage, INT part, INT resource)
{
  if (usage == NULL || usage->ti_part_load == NULL ||
      part < 0 || resource < 0 || part >= usage->i_parts ||
      resource >= usage->i_resources)
    return -1;
  return usage->ti_part_load[part * usage->i_resources + resource];
}
