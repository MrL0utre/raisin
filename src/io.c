/*
 * Copyright (C) 2025 Commissariat à l'énergie atomique et aux énergies
 * alternatives (CEA) and Institut national de recherche en sciences et
 * technologies du numérique (INRIA)
 * Contributor(s): Julien Rodriguez <julien.ro34@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
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
#include <errno.h>

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
  int len;

  if (fgets(buffer, size, file) == NULL)
    return -1;

  len = (int)strlen(buffer);
  if (len > 0 && buffer[len - 1] == '\n')
    buffer[--len] = '\0';
  else
    FATAL(!feof(file), "buffer size too small");

  if (len > 0 && buffer[len - 1] == '\r')
    buffer[--len] = '\0';
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
      free(buffer);
      return 1;
    }

  for(INT i = 0; i < i_vertices; i++) 
    {
      char *end = NULL;
      long value;

      if (read_line(in, buffer, BUFSIZE) <= 0)
        {
          fprintf(stderr, "Partition file %s ends before vertex %d\n", s_path, i);
          free(buffer);
          fclose(in);
          return 2;
        }

      errno = 0;
      value = strtol(buffer, &end, 10);
      if (errno == ERANGE || end == buffer || *end != '\0' ||
          value < 0 || value > RAISIN_PART_MAX)
        {
          fprintf(stderr, "Invalid partition value for vertex %d: %s\n", i, buffer);
          free(buffer);
          fclose(in);
          return 3;
        }
      partition[i] = (PART)value;
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
  static const char suffix[] = ".sol";
  size_t file_name_size = strlen(file_path) + sizeof(suffix);
  char *file_name = (char *)malloc(file_name_size);
  MEM_ERROR(file_name);
  snprintf(file_name, file_name_size, "%s%s", file_path, suffix);

  FILE * out = NULL;
  out = fopen(file_name, "w+");
  
  if (out == NULL) 
    {
      fprintf(stderr, "Cannot open file %s\n", file_name);
      free(file_name);
      return 1;
    }

  for(INT i = 0; i < i_vertices; i++) 
    {
      fprintf(out,"%d\n",partition[i]);
    }

  fclose(out);
  free(file_name);
  
  return (0);
}
