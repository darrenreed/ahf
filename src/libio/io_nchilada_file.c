/**
 * \file io_nchilada_file.c
 *
 * Provides functions for reading and writing NChilada field files
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

///////#include "io_nchilada.h"
#include "io_nchilada_header.h"
#include "io_nchilada_file_def.h"
#include "io_nchilada_file.h"
/////////#include "io_util.h"

////////#include "../common.h"
////////////////// TODO:  should this be a separate file, or just a local function   local_io_nchilada_file_init


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/* io_nchilada_t is general NChilada IO struct 
   io_nchilada_file_t is specific struct for NChilada field data file ("pos","mass" etc.)
*/


/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_nchilada_open more readable.
 *
 * \param log   The logging object.
 * \param f     The NCHILADA file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_nchilada_file_t
local_openopen(io_logging_t log, io_nchilada_file_t f, io_file_mode_t mode);

 /**
  * \brief Helper funtion to set the swap-state and write some
  *        messages
  *
  * \param log      The logging object.
  * \param f        The NCHILADA file object sofar.
  * \param swapped  The swap state.
  *
  * \return Returns the file object or NULL in case of an error.
  */
inline static void
local_openswapped(io_logging_t log,
                  io_nchilada_file_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
  *
  * \return nothing
 */
inline static io_nchilada_file_t
local_opengetswap(io_logging_t log, io_nchilada_file_t f);

/**
 * \brief Check before writing particles and place file pointer at the
 *        right spot.
 *
 * \param log
 * \param f
 * \param pskip
 * \param bytes
 *
 * \return 
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_nchilada_file_t f,
                   uint64_t pskip,
                   int32_t bytes);




extern io_nchilada_file_t
io_nchilada_file_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader)
{

  io_nchilada_file_t f;
  /* Get memory for the structure */
  f = (io_nchilada_file_t)malloc(sizeof(io_nchilada_file_struct_t)); 
  if (f == NULL) {
    io_logging_memfatal(log,  "io_nchilada structure");
    return NULL;
  }

  /* Start filling the structure */

  /* Store the filename */
  f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
  if (f->fname == NULL) {
    io_logging_memfatal(log, "filename of NChilada File");
    free(f);
    return NULL;
  }
  strncpy(f->fname, fname, strlen(fname)+1);

  /* Okay, we are a NChilada file */
  f->ftype = IO_FILE_NCHILADA;


  // TODO does this MPI stuff go here?
	/* And we can just copy in the parallel information */ 
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	if (f->size >= reader) {
		/* TODO 
		 * THIS IS JUST A QUICK HACK TO PREVENT NCHILADA FILES TO
		 * AGAIN TRY TO SPLIT THE COMMUNICATOR, THAT IS ALREADY DONE
		 * IN mnchilada.c ? 
		 * TODO 
		 */
		MPI_Comm_split(MPI_COMM_WORLD, 1,
		               f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->mycomm = MPI_COMM_NULL;
		f->size_mycomm = -1;
		f->rank_mycomm = -1;
	}
#	endif

	fprintf(stderr,"rank %d of %d ready to open file %s %s \n",f->rank, f->size, f->fname, fname); fflush(stderr);////////////////// TODO debugging ////////////

        /* Try to open the file and set mode */
        if (local_openopen(log, f, mode) == NULL) {
                free(f->fname);
                free(f);
                return NULL;
        }
    	fprintf(stderr,"rank %d of %d opening file %s %s \n",f->rank, f->size, f->fname, fname); fflush(stderr);////////////////// TODO debugging ////////////

        /* Set swapping */
	///////////// TODO: get swap/endianness from magic number in header, or always use the AHF -DBYTESWAP flag?
        local_openswapped(log, f, swapped);
        if (    (f->mode == IO_FILE_READ)
             && (f->swapped == IO_FILE_UNKOWN_SWAPPING)) {
                if (local_opengetswap(log, f) != f) {
                        io_logging_fatal(log, "Cannot open this file.");
                        free(f->fname);
                        free(f);
                        return NULL;
                }
        }
	fprintf(stderr,"rank %d of %d opened file %s %s \n",f->rank, f->size, f->fname, fname); fflush(stderr);////////////////// TODO debugging ////////////

	        /* Identify NChilada format */
        local_openversion(log, f); // Is this needed?  Is there more than 1 nchilada version?

        /* Nothing for the header for now */
        f->header = NULL;

        /* Set some dummy values */
        f->no_part           = UINT64_C(0);
        f->no_part_with_mass = UINT64_C(0);
        f->multimass         = INT8_C(0);
        f->mmass             = 1e40;
        f->minweight         = 1e40;
        f->maxweight         = 0.0;
        f->sumweight         = 0.0;
        f->no_species        = INT32_C(0);
        f->posscale          = 1.0;
        f->weightscale       = 1.0;

	fprintf(stderr,"rank %d of %d returning from open file %s %s \n",f->rank, f->size, f->fname, fname); fflush(stderr);////////////////// TODO debugging ////////////

        return f;

}

extern void
io_nchilada_file_close(io_logging_t log,
                io_nchilada_file_t *f)
{


    if ((*f)->file == NULL) {fprintf(stderr,"WARNING file is null but has not been closed %s\n",(*f)->fname);}////////////////// TODO debugging

        /* Catch NULLs */
        if (f == NULL || *f == NULL)
                return;

	        /* Put header to the file if necessary */
  /* XXX Only relevant for writing */  //////////// TODO: do we want to support NCHILADA writing ??
        /* Put header to the file if necessary */
/************************
	if (    ((*f)->mode == IO_FILE_WRITE)
             && ((*f)->header != NULL)) {
                io_nchilada_header_write(log, (*f)->header, *f);
        }
*************************/

	fprintf(stderr,"free nchilada header xxxxxxxxxxxx\n");fflush(stderr);////////////////// TODO debugging
        /* Close */
        if ((*f)->header != NULL)
                io_nchilada_header_del(log, &((*f)->header));
	fprintf(stderr,"free file fname xxxxxxxxxxxx\n");fflush(stderr);////////////////// TODO debugging
        if ((*f)->fname != NULL)
                free((*f)->fname);
	fprintf(stderr,"free file mycomm xxxxxxxxxxxx\n");fflush(stderr);////////////////// TODO debugging
#       ifdef WITH_MPI
        if ((*f)->mycomm != MPI_COMM_NULL)
                MPI_Comm_free(&((*f)->mycomm));
#       endif

    	fprintf(stderr,"actually close the file xxxxxxxxxxxx\n");fflush(stderr);////////////////// TODO debugging

        /* Actually close the file */
        if ((*f)->file != NULL)
                fclose((*f)->file);
        fprintf(stderr,"done actually closed the file dddddddxxxxxxxxxxxx\n");fflush(stderr);////////////////// TODO debugging


        /* Cleaning */
        free(*f);
        *f = NULL;

        return;
}


extern void
io_nchilada_file_init(io_logging_t log,
               io_nchilada_file_t f)
{
        if (f == NULL)
                return;

        if (f->header != NULL) {
                io_logging_warn(log, INT32_C(1),
                                "Already have the header information! Rereading.");
                io_nchilada_header_del(log, &(f->header));
        }

        if (f->mode != IO_FILE_READ) {
                io_logging_warn(log, INT32_C(1),
                                "%s is not opened for reading. "
                                "Will do nothing.",
                                f->fname);
                return;
        }

        io_logging_msg(log, INT32_C(5),
                       "Starting to initialize file object from %s",
                       f->fname);
        f->header = io_nchilada_header_get(log, f);
        io_logging_msg(log, INT32_C(5),
                       "Done with initializing file object from %s",
                       f->fname);

        /* sum up the particle count */
	f->no_part = (uint64_t)f->header->nbodies + ((uint64_t)f->header->iHighWord << 32); // assumes n is a single 64bit int. not 2 signed int32  /////// TODO: check  
        f->multimass = 1;
	f->no_part_with_mass = f->no_part;  // assuming
	f->sumweight = 0.;

////////    fprintf(stderr,"after header read myrank globalrank  name no_part no_part_with_mass no_part_manual iHighWord nbodies  %d %d %s %llu %llu %llu %lu %lu\n",f->rank_mycomm,f->rank,f->fname,f->no_part,f->no_part_with_mass, f->header->nbodies + (f->header->iHighWord << 32),f->header->iHighWord,f->header->nbodies); ////////// TODO debugging //////// TODO debugging verify number of particles on all ranks


        return;
}



/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/


////////// TODO:  delete local functions that are not used.


inline static io_nchilada_file_t
local_openopen(io_logging_t log, io_nchilada_file_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			return NULL;
		}
	}

	f->mode = mode;

	return f;    /////// TODO should return nothing?
}



inline static void
local_openswapped(io_logging_t log,
                  io_nchilada_file_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch(swapped) {
			case IO_FILE_ISNOT_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming unswapped file");
				break;
			case IO_FILE_IS_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming swapped file");
				break;
			case IO_FILE_UNKOWN_SWAPPING:
			default:
				io_logging_msg(log, INT32_C(3),
				               "Will try to find out swap status");
		}
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}


inline static io_nchilada_file_t
local_opengetswap(io_logging_t log, io_nchilada_file_t f)
{

  /* We simply let the user decide! */
  
#ifdef BYTESWAP
  f->swapped = IO_FILE_IS_SWAPPED;
#else
  f->swapped = IO_FILE_ISNOT_SWAPPED;
#endif
  
  
	return f;
}


local_openversion(io_logging_t log, io_nchilada_file_t f)
{
  f->ver = 1;
  return;
}






inline static int32_t
local_write_common(io_logging_t log,
                   io_nchilada_file_t f,
                   uint64_t pskip,
                   int32_t bytes)
{
	return 0;
}


#define CHECK_FLOATBYTES(bfile, bstore) {\
	if (bfile == sizeof(float)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses float for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else if (bfile == sizeof(double)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses double for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else { \
		io_logging_fatal(log, \
		                 "No clue what kind of floating point uses " \
		                 "%" PRIi32 " bytes. Aborting reading.", \
		                 bfile); \
		return UINT64_C(0); \
	}\
	if (bfile < bstore) { \
		io_logging_msg(log, INT32_C(1), \
		               "The floating point values in the file have " \
		               "less precision than the particle storage " \
		               "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		               "No problem, will upcast.", \
		                bfile, bstore); \
	} else if (bfile > bstore) { \
		io_logging_warn(log, INT32_C(1), \
		                "The floating point values in the file have " \
		                "a higher precision than the particle storage " \
		                "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		                "Will downcast, but precision might be lost.", \
		                bfile, bstore); \
	} \
}


