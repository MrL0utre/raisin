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

int archInit(Arch * this, INT m, INT n) {

  this->ti_capacity      = (INT *)calloc(m*n, sizeof(INT));
  MEM_ERROR(this->ti_capacity);

  this->ti_delay         = (INT *)calloc(m*n, sizeof(INT));
  MEM_ERROR(this->ti_delay);

 
  this->i_m = m;
  this->i_n = n;

  return (0);
}

int archFree(Arch * this) {

  free(this->ti_capacity);
  free(this->ti_delay);
  if (this->s_arch_name != NULL) {
    free(this->s_arch_name);
  }
  free(this);
  return (0);
}

int archLoad(Arch * this, const char * s_path, bool verbose)
{
  char * buffer = (char*)malloc(BUFSIZE);
  MEM_ERROR(buffer);

  FILE * in = NULL;
  in = fopen(s_path, "r");
  if (in == NULL) {
    fprintf(stderr, "Cannot open file %s\n", s_path);
    exit(0);
  }

  char ** delimiter = (char**)malloc(sizeof(char*)*1);
  *(delimiter) = (char*)malloc(sizeof(char)*2);
  strcpy(*(delimiter), ".");

  INT i_n, i_m, i_ncon;

  readLine(in, buffer, BUFSIZE);
  
  int n = sscanf(buffer,"%d %d",&i_m,&i_ncon);
  
  i_n = i_m;
  if (n<0) {
    printf("error\n");
  }
  if (verbose) {
    printf("%d %d\n", i_m, i_n); 
  }

  FATAL(n!=2,"Invalid arch file format");

  archInit(this, i_m, i_n);

  INT       capacity   = 0;
  INT       delay      = 0;
  INT       u;
  INT       v;
  char *    raw;
  char      token[2]   = " ";

  for (INT j = 0; j < i_ncon; j++) {

    readLine(in, buffer, BUFSIZE);
    if (verbose) {
      printf("j(%d)\n%s\n", j, buffer);
    }
    raw = strtok(buffer, token);

    capacity = atoi(raw);
    raw      = strtok(NULL, token);
    delay    = atoi(raw);
    raw      = strtok(NULL, token);
    u        = atoi(raw);
    raw      = strtok(NULL, token);
    v        = atoi(raw);
    raw      = strtok(NULL, token);

    this->ti_capacity[u*i_n+v] = capacity;
    this->   ti_delay[u*i_n+v] = delay;

    this->ti_capacity[v*i_n+u] = capacity;
    this->   ti_delay[v*i_n+u] = delay;

  }
  
  free(buffer);
  free(raw);
  fclose(in);

  return (0);
}
