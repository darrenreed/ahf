#ifndef IO_NCHILADA_HEADER_DEF_H 
#define IO_NCHILADA_HEADER_DEF_H

/**
 * \file io_nchilada_header_def.h
 *
 * Provides the structure definition for the NChilada header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdint.h>
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/*
 * The size (in bytes) reserved at the beginning of a file for
 * the header 
 */
#define NCHILADA_HEADER_SIZE 28

/*
 * We actually store additional information in the header structure 
 */
#define NCHILADA_HEADER_EXTRA 48

/*
 * The header structure itself
 */
struct io_nchilada_header_struct {
  uint32_t    magic;
  double time;
  uint32_t    iHighWord;
  uint32_t    nbodies;
  uint32_t    ndim;
  uint32_t    code;
  
  /* this is the extra information not found in the actual header */
  double  omega0;
  double  lambda0;
  double  boxsize;
  double  vunit;
  double  munit;
  double  eunit;
};

/** Convenient typedef */
typedef struct io_nchilada_header_struct io_nchilada_header_struct_t;

/** Convenient typedef */
typedef io_nchilada_header_struct_t *io_nchilada_header_t;


#endif /* IO_NCHILADA_HEADER_DEF_H */
