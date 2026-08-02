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

int arch_init(Arch * this, INT m, INT n) {

  this->ti_capacity      = (INT *)calloc((size_t)m * (size_t)n, sizeof(INT));
  MEM_ERROR(this->ti_capacity);

  this->ti_delay         = (INT *)calloc((size_t)m * (size_t)n, sizeof(INT));
  MEM_ERROR(this->ti_delay);

 
  this->i_m = m;
  this->i_n = n;

  return (0);
}

int arch_free(Arch * this) {

  free(this->ti_capacity);
  free(this->ti_delay);
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
  this->ti_delay = NULL;
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
      free(this->ti_capacity);
      free(this->ti_delay);
      this->ti_capacity = NULL;
      this->ti_delay = NULL;
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
      free(this->ti_capacity);
      free(this->ti_delay);
      this->ti_capacity = NULL;
      this->ti_delay = NULL;
      free(buffer);
      fclose(in);
      return 5;
    }

    this->ti_capacity[u*i_n+v] = capacity;
    this->   ti_delay[u*i_n+v] = delay;

    this->ti_capacity[v*i_n+u] = capacity;
    this->   ti_delay[v*i_n+u] = delay;

  }
  
  free(buffer);
  fclose(in);

  return (0);
}
