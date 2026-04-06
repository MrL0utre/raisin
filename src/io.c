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


/************************************************************/
/**                                                        **/
/**   NAME       : io.c                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are input and output        **/
/**                functions definitions.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10 jun 2022     **/
/**                                                        **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include"io.h"

/**
 * @brief Read a line and removing trailing '\n'.
 *
 * @param file             Input file.
 * @param buffer           String.
 * @param size             Size of string. 
 *
 * @return The length of the readed line.
 */
int 
read_line(FILE * file, 
          char * buffer, 
          int    size) 
{
  int len = 0;

  if(fgets(buffer, size, file)) 
    {
      len = strlen(buffer);
      FATAL(len == 0, "length is zero");
      FATAL(buffer[len-1] != '\n', "buffer size too small");
      buffer[--len] = '\0';
    }
  return len;
}

/**
 * @brief Read a file describing a partition.
 *
 * @param i_vertices       Number of vertices.
 * @param partition        Partition assignation.
 * @param s_path           Path of file. 
 *
 * @return Integer.
 */
int 
load_partition(INT          i_vertices, 
               PART       * partition,
               const char * s_path){

  char * buffer = malloc(BUFSIZE);
  MEM_ERROR(buffer);

  FILE * in = NULL;
  in = fopen(s_path, "r");

  if (in == NULL) 
    {
      fprintf(stderr, "Cannot open file %s\n", s_path);
      exit(0);
    }

  for(INT i = 0; i < i_vertices; i++) 
    {
      read_line(in, buffer, BUFSIZE);
      partition[i] = atoi(buffer);
    }

  free(buffer);
  fclose(in);
  
  return (0);
}

/**
 * @brief Write a file describing a partition.
 *
 * @param i_vertices       Number of vertices.
 * @param partition        Partition assignation.
 * @param s_path           Path of file. 
 *
 * @return Integer.
 */
int 
write_partition(INT          i_vertices, 
                PART       * partition, 
                const char * file_path)
{
  char * buffer = malloc(BUFSIZE);
  MEM_ERROR(buffer);

  char file_name[strlen(file_path) + 5];
  strcpy(file_name, file_path);
  strcat(file_name,".sol");

  FILE * out = NULL;
  out = fopen(file_name, "w+");
  
  if (out == NULL) 
    {
      fprintf(stderr, "Cannot open file %s\n", file_name);
      exit(0);
    }

  for(INT i = 0; i < i_vertices; i++) 
    {
      fprintf(out,"%d\n",partition[i]);
    }

  fclose(out);
  free(buffer);
  
  return (0);
}
