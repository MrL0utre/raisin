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
/**   NAME       : io.h                                    **/
/**                                                        **/
/**   AUTHOR     : Julien RODRIGUEZ                        **/
/**                                                        **/
/**   FUNCTION   : These lines are read and write function **/
/**                declarations.                           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 04 feb 2022     **/
/**                                 to   : 05 apr 2025     **/
/**                                                        **/
/**                                                        **/
/************************************************************/


/** 
 * @file io.h
 * @author Julien Rodriguez
 * @date 10 jun 2022 – 05 apr 2025
 * @brief  These lines are read and write function  
 *        declarations.
 *      
 */
#ifndef IO_H

#define IO_H
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"

int 
read_line(FILE * file, 
          char * buffer, 
          int size);

int 
load_partition(INT i_vertices, 
               PART * partition, 
               const char * s_path);

int 
write_partition(INT i_vertices,
                PART * partition, 
                const char * file_path);

int
write_int_partition(INT i_vertices,
                    const INT *partition,
                    const char *file_path);

#endif
