#ifndef IO_NCHILADA_FILE_H
#define IO_NCHILADA_FILE_H

/**
 * \file io_nchilada_file.h
 *
 * Provides functions for reading and writing NCHILADA files.
 *  at the individual file level (e.g. mass/pos, mass/vel, gas/density
 */
/******  TODO: change fname from single char to char array **/

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

//////////#include "io_nchilada_header_def.h"   //////// TODO: which .h files needed?
//////////#include "io_nchilada_def.h"
/////////#include "io_file.h"
////////////#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Just tries to open a NCHILADA file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename of the NCHILADA file.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according NCHILADA file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a NCHILADA file as the interal mechanism will
 *                  normally detect the right state.
 * \param mode      Tells if the file should be opened for reading or for
 *                  writing. If opened for writing, the value for swapped
 *                  will be ignored.
 * \param reader    Number of processes reading. Only important if in
 *                  MPI mode, otherwise it will be forced to 1.
 *
 * \return Returns a partially initialized file object, or NULL if the
 *         file could not be opened.
 */
extern io_nchilada_file_t
io_nchilada_file_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader);


/**
 * \brief This will close and finalize an NCHILADA file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the NCHILADA file object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_nchilada_file_close(io_logging_t log,
                io_nchilada_file_t *f);

/**
 * \brief Initializes an opened for reading NCHILADA file.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_nchilada_file_init(io_logging_t log,
               io_nchilada_file_t f);



#endif /* IO_NCHILADA_H */
