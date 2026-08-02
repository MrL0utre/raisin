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


#include "a.h"
#include <limits.h>

static void
arch_release_arrays(Arch *this)
{
  free(this->ti_capacity);
  free(this->ti_link_delay);
  free(this->ti_delay);
  free(this->ti_next_hop);
  this->ti_capacity = NULL;
  this->ti_link_delay = NULL;
  this->ti_delay = NULL;
  this->ti_next_hop = NULL;
  this->i_m = 0;
  this->i_n = 0;
}

static int
arch_compute_routes(Arch *this, const char *path)
{
  INT n = this->i_n;

  for (INT i = 0; i < n; i++) {
    for (INT j = 0; j < n; j++) {
      size_t index = (size_t)i * (size_t)n + (size_t)j;
      if (i == j) {
        this->ti_delay[index] = 0;
        this->ti_next_hop[index] = i;
      } else if (this->ti_link_delay[index] >= 0) {
        this->ti_delay[index] = this->ti_link_delay[index];
        this->ti_next_hop[index] = j;
      } else {
        this->ti_delay[index] = INT_MAX;
        this->ti_next_hop[index] = -1;
      }
    }
  }

  for (INT via = 0; via < n; via++) {
    for (INT from = 0; from < n; from++) {
      INT from_via = this->ti_delay[from * n + via];
      if (from_via == INT_MAX)
        continue;
      for (INT to = 0; to < n; to++) {
        INT via_to = this->ti_delay[via * n + to];
        if (via_to == INT_MAX || from_via > INT_MAX - via_to)
          continue;
        INT candidate = from_via + via_to;
        if (candidate < this->ti_delay[from * n + to]) {
          this->ti_delay[from * n + to] = candidate;
          this->ti_next_hop[from * n + to] =
              this->ti_next_hop[from * n + via];
        }
      }
    }
  }

  for (INT from = 0; from < n; from++) {
    for (INT to = 0; to < n; to++) {
      if (this->ti_delay[from * n + to] == INT_MAX) {
        fprintf(stderr,
                "Architecture %s is disconnected: no route from %d to %d\n",
                path, from, to);
        return 1;
      }
    }
  }
  return 0;
}

int arch_init(Arch * this, INT m, INT n) {

  this->ti_capacity      = (INT *)calloc((size_t)m * (size_t)n, sizeof(INT));
  MEM_ERROR(this->ti_capacity);

  this->ti_link_delay    = (INT *)malloc((size_t)m * (size_t)n * sizeof(INT));
  MEM_ERROR(this->ti_link_delay);

  this->ti_delay         = (INT *)calloc((size_t)m * (size_t)n, sizeof(INT));
  MEM_ERROR(this->ti_delay);

  this->ti_next_hop      = (INT *)malloc((size_t)m * (size_t)n * sizeof(INT));
  MEM_ERROR(this->ti_next_hop);

  for (INT i = 0; i < m * n; i++) {
    this->ti_link_delay[i] = -1;
    this->ti_next_hop[i] = -1;
  }

 
  this->i_m = m;
  this->i_n = n;

  return (0);
}

int arch_free(Arch * this) {

  arch_release_arrays(this);
  if (this->s_arch_name != NULL) {
    free(this->s_arch_name);
  }
  free(this);
  return (0);
}

int arch_load(Arch * this, const char * s_path, bool verbose)
{
  char * buffer = (char*)malloc(BUFSIZE);
  MEM_ERROR(buffer);

  this->s_arch_name = NULL;
  this->ti_capacity = NULL;
  this->ti_link_delay = NULL;
  this->ti_delay = NULL;
  this->ti_next_hop = NULL;
  this->i_m = 0;
  this->i_n = 0;

  FILE * in = NULL;
  in = fopen(s_path, "r");
  if (in == NULL) {
    fprintf(stderr, "Cannot open file %s\n", s_path);
    free(buffer);
    return 1;
  }

  INT i_n, i_m, i_ncon;
  char extra;

  if (read_line(in, buffer, BUFSIZE) <= 0) {
    fprintf(stderr, "Architecture file %s has no header\n", s_path);
    free(buffer);
    fclose(in);
    return 2;
  }
  
  int n = sscanf(buffer, "%d %d %c", &i_m, &i_ncon, &extra);
  
  i_n = i_m;
  if (verbose) {
    printf("%d %d\n", i_m, i_n); 
  }

  if (n != 2 || i_m <= 0 || i_ncon < 0) {
    fprintf(stderr, "Invalid architecture header in %s\n", s_path);
    free(buffer);
    fclose(in);
    return 3;
  }

  arch_init(this, i_m, i_n);

  INT       capacity   = 0;
  INT       delay      = 0;
  INT       u;
  INT       v;
  for (INT j = 0; j < i_ncon; j++) {

    if (read_line(in, buffer, BUFSIZE) <= 0) {
      fprintf(stderr, "Architecture file %s ends before connection %d\n", s_path, j);
      arch_release_arrays(this);
      free(buffer);
      fclose(in);
      return 4;
    }
    if (verbose) {
      printf("j(%d)\n%s\n", j, buffer);
    }
    n = sscanf(buffer, "%d %d %d %d %c", &capacity, &delay, &u, &v, &extra);
    if (n != 4 || capacity < 0 || delay < 0 ||
        u < 0 || u >= i_n || v < 0 || v >= i_n || u == v) {
      fprintf(stderr, "Invalid architecture connection %d in %s\n", j, s_path);
      arch_release_arrays(this);
      free(buffer);
      fclose(in);
      return 5;
    }

    if (this->ti_link_delay[u * i_n + v] >= 0) {
      fprintf(stderr, "Duplicate architecture connection %d-%d in %s\n",
              u, v, s_path);
      arch_release_arrays(this);
      free(buffer);
      fclose(in);
      return 6;
    }

    this->ti_capacity[u*i_n+v] = capacity;
    this->ti_link_delay[u*i_n+v] = delay;

    this->ti_capacity[v*i_n+u] = capacity;
    this->ti_link_delay[v*i_n+u] = delay;

  }
  
  if (arch_compute_routes(this, s_path) != 0) {
    arch_release_arrays(this);
    free(buffer);
    fclose(in);
    return 7;
  }

  free(buffer);
  fclose(in);

  return (0);
}

bool
arch_has_link(const Arch *this, INT u, INT v)
{
  if (this == NULL || this->ti_link_delay == NULL ||
      u < 0 || v < 0 || u >= this->i_m || v >= this->i_n || u == v)
    return false;
  return this->ti_link_delay[u * this->i_n + v] >= 0;
}

INT
arch_next_hop(const Arch *this, INT u, INT v)
{
  if (this == NULL || this->ti_next_hop == NULL ||
      u < 0 || v < 0 || u >= this->i_m || v >= this->i_n)
    return -1;
  return this->ti_next_hop[u * this->i_n + v];
}
