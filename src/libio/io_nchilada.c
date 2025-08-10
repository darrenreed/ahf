/**
 * \file io_nchilada.c
 *
 * Provides functions for reading and writing NCHILADA files.
 */
/******** TODO     
             --change the 32 bit nparticle headers into 64 bit	     ?
             --change io_nchilada_field to io_nchilada_file   because a nchilada field is always a file  ?
	     --check endianness/swap using nchilada header magic number ?
             --should we be compatible with pseudo-nchilada formats that have no "iord" files?
                    or the user should need to creat a dummy iord files (e.g. ordered 1-npart)?
***/
/*
 *
 *  Requires that sim base filename not contain "gas" "dark" or "star"
 *     or "pos" or "vel" or "id" or "temperature" or "iord" or "mass" or "GasDensity" or "soft"
 *
 *  ...anyways, please we very careful with the usage of the
 *  thermal energy of gas particles: AHF expect it to be in km/sec!
 *          (if in doubt, please get in touch...)
 *
 *
 */

/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_nchilada.h"
#include "io_nchilada_header.h"
#include "io_nchilada_file.h"
#include "io_util.h"

#include "../common.h"

/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/

#ifdef FSEEKO    // needed for input files more than 2**31 bytes 
#define _FILE_OFFSET_BITS 64
#define _LARGE_FILE_API
#endif

#define H0              100.  // not known as param.h is not included

#define MAX(A,B)        ((A)>(B)?(A):(B))

// internally used bytes for various particle types
#define NCHILADA_GASVAR    12
#define NCHILADA_DARKVAR    9
#define NCHILADA_STARVAR   11

// ptype values as used by AHF
// NOTE: this has to be identical to the definitions found in param.h
#define PGAS               ( 0.0)   /* identifier for gas particles; has to be exactly 0.0!!!! */
#define PDM                (-1.0)   /* identifier for dm particles; whatever negative value */
#define PSTAR              (-4.0)   /* identifier for star particles; whatever negative value */
#define PDMbndry           (-5.0)   /* identifier for star particles; whatever negative value */

/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/* io_nchilada_t is general NChilada IO struct 
   io_nchilada_file_t is specific struct for NChilada field data file ("pos","mass" etc.)
*/


/**
 *  \brief  identifier to help keep track of which particle field goes with which file
 *
 *   values correspond to the NChilada files that AHF needs to read
 *      same values read as for tipsy, except that gas softening is added
 *      iord not used
 *
 *   Should this enum be in this file?
 */
typedef enum {
    GASMASS,        // 0
    GASPOS,         // 1
    GASVEL,         // 2
    GASTEMPERATURE, // 3
    GASID,
#ifdef STORE_MORE
    GASDENSITY,  
    GASSOFT,
#endif    
    DARKMASS, 
    DARKPOS,  
    DARKVEL,
    DARKID,
#ifdef STORE_MORE    //  !DARK and STAR density are never set (nor in io_tipsy.c)
    DARKDENSITY,     //  !STORE_MORE should not be used?
    DARKSOFT,  
#endif
    STARMASS, 
    STARPOS,
    STARVEL,
#ifndef STORE_MORE   
    STARID
#endif
#ifdef STORE_MORE
    STARID,
    STARDENSITY,  
    STARSOFT
#endif    
} particle_field_enum;


/**
 * \brief Check whether the open file (normally description.xml) 
 *        contains the particle family name
 *
 * \param file          expects description.xml
 * \param family_name   (NChilada can have gas, dark, or star)
 *
 * \return Returns a 1 or 0 indicating whether family name is present in snapshot
 */
static int32_t
local_check_family(FILE *file,
                   const char *family);
		   
/**
 * \brief Check whether this file is gas particles from filename 
 *
 * \param filename
 *
 * \return Returns a 1 or 0  (true or false)
 */
static int32_t
local_i_am_gasfile(const char *filename);


/**
 * \brief Check whether this file is dark particles from filename 
 *
 * \param filename
 *
 * \return Returns a 1 or 0  (true or false)
 */
static int32_t
local_i_am_darkfile(const char *filename);
  
/**
 * \brief Check whether this file is star particles from filename 
 *
 * \param filename
 *
 * \return Returns a 1 or 0  (true or false)
 */
static int32_t
local_i_am_starfile(const char *filename);


// identify species from filename.  assumes NChilada naming conventions ///////// IS THIS NEEDED???
  /**
 * \brief Check whether this file is star particles from filename 
 *
 * \param filename
 *
 * \return Returns a -4 -1 or 0     (-4:star  -1:dark  0:gas)
 */
static int32_t
local_get_species_tag_from_filename(const char *filename);
  

// identify species from filename.  assumes NChilada naming conventions ///////// IS THIS NEEDED??? /////////// TODO: delete this function  (instead use local_get_species_tag_from_filename)
static int32_t
local_get_species_from_filename(const char *filename);

  
///////int32_t local_get_species_from_filename() TODO

/**
 * \brief get names for all NChilada snapshot files to open 
 *        
 * \param log
 * \param path
 * \param fnames 
 *
 * \return  Returns number of files
 */
static int32_t
local_findfiles(io_logging_t log,
		  const char *path,
                  char ***fnames);


/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_pos(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);

/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_vel(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);

/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_mass(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);


/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_u(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);


/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_density(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);


/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_softening(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);



/**
 * \brief read particle positions field from the open file
 *        
 * \param log
 * \param f
 * \param pskip
 * \param pread
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_id(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);


/**
 * \brief set particle field u to value tag to indicate dark or star
 *        
 * \param log
 * \param pread
 * \param tag
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_u_nongas_label(io_logging_t log,
                    uint64_t *pread,
		    double tag,
                    io_file_strg_struct_t strg);

/**
 * \brief set particle field density to value tag (instead of reading from input)
 *        
 * \param log
 * \param pread
 * \param tag
 * \param strg
 *
 * \return  Returns number of files
 */
static uint64_t
local_get_block_density_nongas_label(io_logging_t log,
                    uint64_t *pread,
		    double tag,
                    io_file_strg_struct_t strg);

//////local_io_nchilada_init    ////////// TODO?
  

  

/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_nchilada_t
io_nchilada_open(io_logging_t log,
	       char *fname,            /* AHF.input: ic_filename==path to description.xml */
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader)
{
        int32_t i;
        io_nchilada_t f;
	char **fnames;

#ifdef VERBOSE

#endif	
	
#ifndef GAS_PARTICLES
  fprintf(stderr,"you have not defined -DGAS_PARTICLES for your NCHILADA simulation\n");
  fprintf(stderr,"            PLEASE RECOMPILE AHF!\n");
  return NULL;
#endif

#ifndef FSEEKO    // needed for input files more than 2**31 bytes 
  fprintf(stderr,"WARNING: You have not defined -DFSEEKO\n");
  fprintf(stderr,"    ahf will probably fail unless all files are smaller than ~2GB\n");
#endif
  
	/* Get memory for the structure */
	f = (io_nchilada_t)malloc(sizeof(io_nchilada_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_nchilada structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->path = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->path == NULL) {
		io_logging_memfatal(log, "filename of NCHILADAfile");
		free(f);
		return NULL;
	}
	strncpy(f->path, fname, strlen(fname)+1); 

	/* Okay, we are a NCHILADA file */
	f->ftype = IO_FILE_NCHILADA;


	  /* And we can just copy in the parallel information */
#       ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1, f->rank, &(f->mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
#       endif

	
	/* Get the filenames */
	////  #if (defined WITH_MPI && defined BCASTHEADER) // TODO: ? only master task reads headers, as in io_mgadget.c

	f->numfiles = local_findfiles(log, f->path, &fnames);
	fprintf(stderr,"found %d files in path %s\n", f->numfiles, f->path);////////////////// TODO debugging
	
        /* Glue the separate field files into the NChilada structure, similar to io_mgadget vs io_gadget */
        f->files = (io_nchilada_file_t *)malloc( sizeof(io_nchilada_file_t)*(f->numfiles)); 
	if (f->files == NULL) {
	  io_logging_memfatal(log,  "io_nchilada structure (2)");
                for (i=0; i<f->numfiles; i++)
                        free(fnames[i]);
                free(fnames);
                free(f->path);
                free(f);
                return NULL;
        }

	/* open all files on all mpi ranks even though it might not be needed */
	for (i=0; i<f->numfiles; i++) {
#ifdef DEBUG_NCHILADA
	  io_logging_msg(log,6,"\ntrying to open file %s (%d of %d) for reading ... ", fnames[i],i,f->numfiles);
#endif

#               ifdef WITH_MPI
      /* TODO
       * THIS   IS JUST A NASTY HACK TO PREVENT io_gadget.c FROM REDOING
       * THE MPI-SPLIT. CALLED FROM HERE IT IS ONLY SUPPOSED TO ACT AS
       * A DUMMY INTERFACE TO A GADGET FILE.
       * TODO
       */
      (f->files)[i] = io_nchilada_file_open(log, fnames[i], swapped, mode, f->size+1);
      
      ///////      fprintf(stderr,"fff i %d f_ver %d f_nopart_init %lld\n",i,((f->files)[i])->ver, ((f->files)[i])->no_part); fflush(stderr);////////////////// TODO debugging 
#               else
      (f->files)[i] = io_nchilada_file_open(log, fnames[i], swapped, mode, reader);
#               endif

      if ((f->files)[i] == NULL) {
        int32_t j;
        for (j=i; i<f->numfiles; j++)
          free(fnames[j]);
        free(fnames);
        while (i>0) {
          i--;
          io_nchilada_file_close(log, &((f->files)[i]));
        }
        free(f->path);
        free(f);
        return NULL;
      } // end if files[i] == NULL

      fprintf(stderr,"rank %d of %d opened files %d %s %s \n",f->rank,f->size,i,((f->files)[i])->fname, fnames[i]); fflush(stderr);////////////////// TODO debugging ////////////

      
    } // end loop over all f->numfiles 

    
    /* the whole char *fnames[] is not needed anymore */ //////////// TODO --- need to put fname into f->files?
    for (i=0; i<f->numfiles; i++) {
      free(fnames[i]);
    }
    free(fnames);

    ////////////////////////////////////////////////  TODO:  only need to read??    Are file names already known to all tasks?

      
	/* Set some dummy values */
	f->no_part = UINT64_C(0);
	f->no_part_with_mass = UINT64_C(0);
	f->multimass = INT8_C(0);
	f->mmass = 1e40;
	f->minweight = 1e40;
	f->maxweight = 0.0;
	f->sumweight = 0.0;
	f->no_species = INT32_C(0);
	f->posscale = 1.0;
	f->weightscale = 1.0;

	return f;
}

extern void
io_nchilada_close(io_logging_t log,
                io_nchilada_t *f)
{ 
  int32_t i;

  
  /* Catch NULLs */
  if (f == NULL || *f == NULL)
    return;

	/* Actually close the file */
	  /* free() all related structures and fclose() the files (in case the file pointer is not NULL) */
  for (i=0; i<(*f)->numfiles; i++) {
    io_nchilada_file_close(log, &(((*f)->files)[i]));
  }
  if((*f)->files) free((*f)->files);
  if((*f)->path)  free((*f)->path);  
  if(*f)          free(*f);
  *f = NULL;

  return;
}


///////////////////////// DONE?
extern void
io_nchilada_init(io_logging_t log,
               io_nchilada_t f)
{ 
  int32_t i;
  uint64_t ngas, ndark, nstar;
  
  if (f == NULL)
    return;

  if (f->files[0]->mode != IO_FILE_READ) {
    io_logging_warn(log, INT32_C(1),
		    "%s (first file of %" PRIi32 ") is not opened "
		    "for reading. Will do nothing.",
		    f->files[0]->fname, f->numfiles);
    return;
  }

#if (defined WITH_MPI && defined BCASTHEADER)
fprintf(stderr,"You have defined WITH_MPI and BCASTHEADER is defined rank %d of %d    myrank %d of %d  numfiles %d\n",f->rank_mycomm,f->size_mycomm, f->rank, f->size, f->numfiles); ////////// TODO debugging 
#endif

  
  ngas=0;
  ndark=0;
  nstar=0;
  f->no_part=0;
  
#if (defined WITH_MPI && defined BCASTHEADER)
  if(f->rank_mycomm == 0) {
#endif
    for (i=0; i<f->numfiles; i++) {      
      io_logging_msg(log, INT32_C(5),
		     "Starting to initialize file object from %s",
		     f->files[i]->fname);
      io_nchilada_file_init(log, (f->files)[i]); // reads header 
      io_logging_msg(log, INT32_C(5),
		     "Done with initializing file object from %s",
		     f->files[i]->fname);
      /////////      fprintf(stderr,"header myrank globalrank i name box omega lambda vunit munit eunit %d %d %d %s %g %g %g %g %g %g\n",f->rank_mycomm,f->rank,i,f->files[i]->fname,f->files[i]->header->boxsize,f->files[i]->header->omega0,f->files[i]->header->lambda0,f->files[i]->header->vunit,f->files[i]->header->munit,f->files[i]->header->eunit); ////////// TODO debugging /////// only rank 0 gets here.  why?????
      
      /* get total particle number by counting particles in files named "pos" (PROBABLY UNSAFE) */ 
      if (strstr(f->files[i]->fname, "pos") != NULL) { // e.g. gets tripped 3 times. once for gas,stars,dark
	f->no_part += f->files[i]->no_part;  // we do not yet need ngas, ndark, nstar -- just their sum
      }
    } // end i<numfiles

    //////////////// TODO debugging   NOW only task 0 has file[0]->header->boxsize   FIX because they all need it for scaling
    
    
    f->multimass = 1;
    f->no_part_with_mass = f->no_part; // assuming NO massless particles
#if (defined WITH_MPI && defined BCASTHEADER)
  }  // endif f->rank_mycomm == 0
  else {
    /* allocate memory on this MPI task to hold the header for each file */
    for (i=0; i<f->numfiles; i++) {
      // TODO debugging careful: "header" is more than just NChilada file format header
      ((f->files)[i])->header = (io_nchilada_header_t)malloc(sizeof(io_nchilada_header_struct_t));  // with mpi, rank 0 header malloc in io_nchilada_file_init
    }
  } 

    /* Broadcast the header of each file from rank==0 to all other MPI tasks */ 
  for (i=0; i<f->numfiles; i++) {
    MPI_Bcast(((f->files)[i])->header,sizeof(io_nchilada_header_struct_t),MPI_BYTE,0,f->mycomm); 

    /* io_nchilada_file_init() also set some more relevant values in f */   
    MPI_Bcast(&(((f->files)[i])->multimass),           sizeof(int8_t),   MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(((f->files)[i])->no_part),             sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);  // number particles in file
    MPI_Bcast(&(((f->files)[i])->no_part_with_mass),   sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(f->multimass),                         sizeof(int8_t),   MPI_BYTE, 0, f->mycomm);  // (this is unnecessarily broadcasted multiple times)
    MPI_Bcast(&(f->no_part),                           sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);  // number particles total
    ////////    fprintf(stderr,"after init Bcast header myrank globalrank i name box omega lambda vunit munit eunit %d %d %d %s %g %g %g %g %g %g\n",f->rank_mycomm,f->rank,i,f->files[i]->fname,f->files[i]->header->boxsize,f->files[i]->header->omega0,f->files[i]->header->lambda0,f->files[i]->header->vunit,f->files[i]->header->munit,f->files[i]->header->eunit); ////////// TODO debugging
    ////////    fprintf(stderr,"after init Bcast header myrank globalrank i name no_part iHighWord nbodies %d %d %d %s %llu %lu %lu\n",f->rank_mycomm,f->rank,i,f->files[i]->fname,f->files[i]->no_part,f->files[i]->header->iHighWord,f->files[i]->header->nbodies); ////////// TODO debugging //////// TODO debugging verify number of particles on all ranks
  } 
#endif



  
  return;
}

extern uint64_t
io_nchilada_readpart(io_logging_t log,
                   io_nchilada_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg)
{
	uint64_t particles_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

  /* 
	 * First read the particles unscaled. This will set important
	 * scaling information
	 */
	particles_read = io_nchilada_readpart_raw(log, f, pskip, pread, strg); 


  if (particles_read != pread) {
		return UINT64_C(0);
	}

  /* And do the scaling */
#ifdef WITH_MPI
	io_nchilada_scale_global(log, f->mycomm,  f->maxpos, f->minpos,
                        &(f->mmass), &(f->minweight), &(f->maxweight), &(f->sumweight));
#endif
  
  //f->mmass = 1.0;
	/* Assume that every NChilada has same sim data in header (should check) */
  tmp = io_nchilada_scale_particles(log,
				    f->maxpos, f->minpos,
				    f->files[0]->header->boxsize,
				    f->files[0]->header->time,
				    f->posscale,
				    f->files[0]->header->vunit,
				    f->files[0]->header->eunit,
				    f->mmass,
				    particles_read, strg);
	if (tmp != particles_read) {
		return tmp;
	}

  /* Wow, we are done! */
	return particles_read;
}


/* No we define a bunch of macros to make life easier */
                                  /* Most of these are not useful for NChilada */
// #define SKIP {io_util_readint32(f->file, &blocksize, f->swapped);} 
// #define SKIP2 {io_util_readint32(f->file, &blocksize2, f->swapped);}
/***
#define CHECK_BLOCK {				\
	if (blocksize != blocksize2) {\
		io_logging_fatal(log,\
		                 "The block boundaries (beginning: %" PRIi32\
		                 " end: %" PRIi32 ") are not identical. "\
		                 "Corrupt file?", blocksize, blocksize2);\
/		return UINT64_C(0);\
	} else {\
		io_logging_msg(log, INT32_C(5),\
		               "Block claimed correctly to be %f MB long.", \
		               (float)(blocksize/1024./1024.));\
	}\
}
#define DESCRIBE_BLOCK {\
	if (f->ver == 2) {\
		SKIP;\
		io_util_readstring(f->file, str, (size_t)4);\
		io_util_readint32(f->file, &nextblocksize, f->swapped);\
		io_logging_msg(log, INT32_C(1),\
		               "Arrived at block %s, size of it will be %" \
		               PRIi32, str, nextblocksize);\
		SKIP2;\
		CHECK_BLOCK;\
	}\
}
***/


  /* NChilada has a single file for each species and each field. 
        e.g. 1 file for gas positions. */
extern uint64_t
io_nchilada_readpart_raw(io_logging_t log,
                       io_nchilada_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg)
{
	uint64_t j, psum;
	int ptype;
	uint32_t bytes_file, bytes_int_file;
#ifdef FSEEKO
	int64_t blocksize, blocksize2, partsize;
#else
	long blocksize, blocksize2, partsize;
#endif
	int32_t nextblocksize;
	long skipsize;
	double fweight=0., oldfweight=0.;
	double fu;
  double ftemp;
  double frho, feps;
	uint64_t fid;
	float dummy;
	uint32_t dummy_int;
	char str[5];
	/** Used to figure out at which particle type we are */
	uint32_t curprtt;     ///////////// TODO: delete?
	uint64_t ndark, ngas, nstar;
	uint32_t ispecies;
	uint64_t pskip_file, pread_file, pread_todo, pread_done, pread_function_return; // 
	int32_t i;   // max i will be numfiles
	
	/* Check if we actually have to do something */
	if ( (f == NULL) || (f->files[0]->header == NULL) )
		return UINT64_C(0);

	/* Need to read gas, dark, and star particles with separate calls */
	ngas = 0;
	ndark = 0;
	nstar = 0;

	for (i=0; i<f->numfiles; i++) {    // 1st time, just get the number of particles from headers
	  /* */   
	  /******************************************************************* \
	   *  Skip HEADER so ready to read particles                           *
          \*******************************************************************/
	  skipsize = NCHILADA_HEADER_SIZE;
	  fseek(f->files[i]->file, skipsize, SEEK_SET);
	  
    /* get number of particles of each species/family (gas, dark, star) 
        by checking files named "pos" (PROBABLY UNSAFE) 
          already have total (no_part) from io_nchilada_init */	  
	  if (strstr(f->files[i]->fname, "pos") != NULL) { // e.g. gets tripped 3 times. once for gas,stars,dark
	    if (local_i_am_gasfile(f->files[i]->fname)) {  
	      ngas = f->files[i]->no_part;
	    }
	    else if (local_i_am_darkfile(f->files[i]->fname)) {
	      ndark = f->files[i]->no_part;
	    }
	    else if (local_i_am_starfile(f->files[i]->fname)) {
	      nstar = f->files[i]->no_part;
	    }       
	    else { //  particle family not found
	      io_logging_fatal(log,
			       "Particle family in filename '%s' not recognized.",
			       f->files[i]->fname);
	    }
	  } // close if pos
	  fprintf(stderr,"file ngas ndark nstar %s %llu %llu %llu\n", f->files[i]->fname,ngas,ndark,nstar); //////////////////////// TODO debugging
      } // close of numfiles loop


	/*********************************************************************** \
	* 1st read gas, then dark, then stars (pskip, pread, pdark)            *
        *                                                                      *
        * This is all quite a mess, but maybe it works.                        *
	/***********************************************************************/

	/* Set extreme position detectors */
	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;

   /* Three particle types, 3 sets of files to read: gas, dark, star */
   /* Read 1st dark, then gas, then star */
   /******************************************************************* \
    *  SKIP to the first particle to be read                          *
   \*******************************************************************/
   // first particle lies in GAS block -- we read from the gas files
  //==================================
	pread_done = 0;
	pread_todo = pread;
	pskip_file = 0;  // might be safer to have separate pskip_gas_file, pskip_dark_file, pskip_star_file etc.
	pread_file = 0;
	
	    
     if (pskip < ngas) { //   this task will read gas particles 
       // how many gas particles to read
       pskip_file = pskip;
       if ((pread+pskip) > ngas) { // will also read other particle types
	 pread_file = ngas - pskip_file;
       }
       else {
	 pread_file = pread;
       }
       ////////////////////////// TODO: if ngas,ndark or nstar = 0, verify: does this all work?
       /* This could be made cleaner */
       for (i=0; i<f->numfiles; i++) { // read the particles

	 fprintf(stderr,"gas i file %d %s rank %d\n",i,f->files[i]->fname,f->rank); fflush(stderr); ////////////// TODO debugging
	 
	 if (local_i_am_gasfile(f->files[i]->fname)) {  
	   if (strstr(f->files[i]->fname, "mass") != NULL) { //  read masses
	     pread_function_return = local_get_block_mass(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle masses, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     
	     // no accumulation of sumweight as we do an MPI_SUM reduction below. 
	     f->sumweight += f->files[i]->sumweight;  // Here is summed only local gas+dark+star
	     if (isless(f->files[i]->minweight, f->minweight))          f->minweight = f->files[i]->minweight;
	     if (isgreater(f->files[i]->maxweight, f->maxweight)) f->maxweight = f->files[i]->maxweight;	
	   }
	   
	   if (strstr(f->files[i]->fname, "pos") != NULL) { //
	     pread_function_return = local_get_block_pos(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle positions, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     if (isgreater(f->files[i]->maxpos[0], f->maxpos[0])) f->maxpos[0] = f->files[i]->maxpos[0];
	     if (isgreater(f->files[i]->maxpos[1], f->maxpos[1])) f->maxpos[1] = f->files[i]->maxpos[1];
	     if (isgreater(f->files[i]->maxpos[2], f->maxpos[2])) f->maxpos[2] = f->files[i]->maxpos[2];
	     if (isless(f->files[i]->minpos[0], f->minpos[0]))    f->minpos[0] = f->files[i]->minpos[0];
	     if (isless(f->files[i]->minpos[1], f->minpos[1]))    f->minpos[1] = f->files[i]->minpos[1];
	     if (isless(f->files[i]->minpos[2], f->minpos[2]))    f->minpos[2] = f->files[i]->minpos[2];	     
	   }
	   
	   if (strstr(f->files[i]->fname, "vel") != NULL) { // 		    
	     pread_function_return = local_get_block_vel(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle velocities, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
	   if (strstr(f->files[i]->fname, "temperature") != NULL) { //
	     pread_function_return = local_get_block_u(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" energies, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
	   if (strstr(f->files[i]->fname, "iord") != NULL) { //
	     pread_function_return = local_get_block_id(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" energies, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
		   
#ifdef STORE_MORE
	   if (strstr(f->files[i]->fname, "GasDensity") != NULL) { // 		    
	     pread_function_return = local_get_block_density(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" densities, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
	   if (strstr(f->files[i]->fname, "soft") != NULL) { // 		    
	     pread_function_return = local_get_block_softening(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" softenings, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
#endif
	   
	 } // end of gas_file read
       } // end of loop over all files       
       pread_done += pread_file;
       pread_todo -= pread_file;
	
      /* Move the particle pointers to match the pointer increments in local_get_block_ */
       
      strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride*pread_file);
      strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride*pread_file);
      strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride*pread_file);
      strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride*pread_file);
      strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride*pread_file);
      strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride*pread_file);
      if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride*pread_file);
      if (strg.id.val     != NULL) strg.id.val     = (void *)(((char *)strg.id.val)     + strg.id.stride*pread_file);
      if (strg.u.val      != NULL) strg.u.val      = (void *)(((char *)strg.u.val)      + strg.u.stride*pread_file);
#ifdef STORE_MORE
      strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride*pread_file);
      strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride*pread_file);
#endif

       
     } // end of if this task reads gas particles


     ////////////     fprintf(stderr,"rank %d of %d gas pskip pread pread_done pread_todo pread_file pskip_file %llu %llu %llu %llu %llu %llu\n", f->rank, f->size,pskip,pread,pread_done,pread_todo,pread_file,pskip_file); //////////////////////// TODO debugging
     pskip_file = 0; 
     pread_file = 0;
       
     /////////////////////////////////// TODO: debugging should not set nongas tags if no dark particles ////////
     
     // read dark matter
     //================================== // have dark particles to read 
     /* last particle is not gas.  and 1st particle is not star. (order is gas, dark, star)   */
     if ( ((pskip + pread) > ngas) && (pskip  < (ngas + ndark)) ) { // will read dark particles (last particle to read is not gas. 1st particle to read is not star)
       if (pskip > ngas) { // 1st particle to read is dark 
	 pskip_file = pskip - ngas;
       }
       else {
	 pskip_file = 0; 
       }
       /* get how many dark particles to read by this task */
       if ((pskip + pread) > (ngas + ndark)) { // last particle is star
	 pread_file = ndark - pskip_file;
       }
       else { // last particle is dark
	 pread_file = pread_todo;
       }


       for (i=0; i<f->numfiles; i++) { // scan all filenames

	 fprintf(stderr,"dark i file %d %s rank %d\n",i,f->files[i]->fname,f->rank); fflush(stderr); ////////////// TODO debugging

	 if (local_i_am_darkfile(f->files[i]->fname)) {  // dark files
	   if (strstr(f->files[i]->fname, "mass") != NULL) { //  read masses		    
	     pread_function_return = local_get_block_mass(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle masses, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     // no accumulation of sumweight as we do an MPI_SUM reduction below. 
	     f->sumweight += f->files[i]->sumweight;  // Here is summed only local gas+dark+star
	     /* mmass: lowest mass dark matter particle */
	     if (isless(f->files[i]->mmass, f->mmass))                  f->mmass     = f->files[i]->mmass; 
	     if (isless(f->files[i]->minweight, f->minweight))          f->minweight = f->files[i]->minweight;
	     if (isgreater(f->files[i]->maxweight, f->maxweight)) f->maxweight = f->files[i]->maxweight;
	   }
	   if (strstr(f->files[i]->fname, "pos") != NULL) { // 		    
	     pread_function_return = local_get_block_pos(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle positions, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     if (isgreater(f->files[i]->maxpos[0], f->maxpos[0])) f->maxpos[0] = f->files[i]->maxpos[0];
	     if (isgreater(f->files[i]->maxpos[1], f->maxpos[1])) f->maxpos[1] = f->files[i]->maxpos[1];
	     if (isgreater(f->files[i]->maxpos[2], f->maxpos[2])) f->maxpos[2] = f->files[i]->maxpos[2];
	     if (isless(f->files[i]->minpos[0], f->minpos[0]))    f->minpos[0] = f->files[i]->minpos[0];
	     if (isless(f->files[i]->minpos[1], f->minpos[1]))    f->minpos[1] = f->files[i]->minpos[1];
	     if (isless(f->files[i]->minpos[2], f->minpos[2]))    f->minpos[2] = f->files[i]->minpos[2];	     
	   }
	   
	   if (strstr(f->files[i]->fname, "vel") != NULL) { // 		    
	     pread_function_return = local_get_block_vel(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle velocities, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	    
	   if (strstr(f->files[i]->fname, "iord") != NULL) { //
	     pread_function_return = local_get_block_id(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" energies, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
#ifdef STORE_MORE
	   if (strstr(f->files[i]->fname, "soft") != NULL) { // 		    
	     pread_function_return = local_get_block_softening(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" softenings, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
#endif
	 } // end of dark files read
       } // end of loop over all files
       local_get_block_u_nongas_label(log, &pread_file, PDM, strg); // label particle type as dark in energy block
#ifdef STORE_MORE
       local_get_block_density_nongas_label(log, &pread_file, PDM, strg); // TODO: Init density to what?  Is this needed at all?
#endif
       pread_done += pread_file;
       pread_todo -= pread_file;


      /* Move the particle pointers */
      strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride*pread_file);
      strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride*pread_file);
      strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride*pread_file);
      strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride*pread_file);
      strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride*pread_file);
      strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride*pread_file);
      if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride*pread_file);
      if (strg.id.val     != NULL) strg.id.val     = (void *)(((char *)strg.id.val)     + strg.id.stride*pread_file);
      if (strg.u.val      != NULL) strg.u.val      = (void *)(((char *)strg.u.val)      + strg.u.stride*pread_file);
#ifdef STORE_MORE
      strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride*pread_file);
      strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride*pread_file);
#endif




     } // end of this task reads dark matter
     fprintf(stderr,"rank %d of %d dark pskip pread pread_done pread_todo pread_file pskip_file %llu %llu %llu %llu %llu %llu\n", f->rank, f->size,pskip,pread,pread_done,pread_todo,pread_file,pskip_file); //////////////////////// TODO debugging
     pskip_file = 0; 
     pread_file = 0;

     // read stars     
     //================================== // have star particles to read 
                  /* last particle to read is not gas or dark */
     if ( ((pskip + pread) > (ngas + ndark)) ) { // will read star particles
       if (pskip > (ngas + ndark)) {      // 1st particle is star
	 pskip_file = pskip - ngas - ndark;
       }
       else {
	 pskip_file = 0;
       }
       /* get how many star particles to read by this task */
       if ((pskip + pread) == nstar) { // last particle is last star
	 pread_file = nstar - pskip_file; // (pread_file == pread_todo)
       }
       else { 
	 pread_file = pread_todo - pskip_file;
       }
       



       for (i=0; i<f->numfiles; i++) { // scan all filenames

	 fprintf(stderr,"star i file %d %s rank %d\n",i,f->files[i]->fname,f->rank); fflush(stderr); ////////////// TODO debugging

	 if (local_i_am_starfile(f->files[i]->fname)) {  // star files
	   if (strstr(f->files[i]->fname, "mass") != NULL) { //  read masses		    
	     pread_function_return = local_get_block_mass(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle masses, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     // no accumulation of sumweight as we do an MPI_SUM reduction below. 
	     f->sumweight += f->files[i]->sumweight;  // Here is summed only local gas+dark+star
	     if (isless(f->files[i]->minweight, f->minweight))          f->minweight = f->files[i]->minweight;
	     if (isgreater(f->files[i]->maxweight, f->maxweight)) f->maxweight = f->files[i]->maxweight;
	   }
	   if (strstr(f->files[i]->fname, "pos") != NULL) { // 		    
	     pread_function_return = local_get_block_pos(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle positions, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	     if (isgreater(f->files[i]->maxpos[0], f->maxpos[0])) f->maxpos[0] = f->files[i]->maxpos[0];
	     if (isgreater(f->files[i]->maxpos[1], f->maxpos[1])) f->maxpos[1] = f->files[i]->maxpos[1];
	     if (isgreater(f->files[i]->maxpos[2], f->maxpos[2])) f->maxpos[2] = f->files[i]->maxpos[2];
	     if (isless(f->files[i]->minpos[0], f->minpos[0]))    f->minpos[0] = f->files[i]->minpos[0];
	     if (isless(f->files[i]->minpos[1], f->minpos[1]))    f->minpos[1] = f->files[i]->minpos[1];
	     if (isless(f->files[i]->minpos[2], f->minpos[2]))    f->minpos[2] = f->files[i]->minpos[2];
	   }
	   
	   if (strstr(f->files[i]->fname, "vel") != NULL) { // 		    
	     pread_function_return = local_get_block_vel(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" particle velocities, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
	   if (strstr(f->files[i]->fname, "iord") != NULL) { //
	     pread_function_return = local_get_block_id(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" energies, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
	   
#ifdef STORE_MORE
	   if (strstr(f->files[i]->fname, "soft") != NULL) { // 		    
	     pread_function_return = local_get_block_softening(log, f->files[i], &pskip_file, &pread_file, strg);
	     if (pread_function_return != pread_file) {
	       io_logging_fatal(log, 
				"Expected to read %"PRIu64
				" softenings, but only got %"PRIu64
				".  Aborting.", pread_file, pread_function_return);
	       return UINT64_C(0);
	     }
	   }
#endif
	 } // end of star_file read
       } // end of loop over all files
       local_get_block_u_nongas_label(log, &pread_file, PSTAR, strg); // label particle type as star in energy block
#ifdef STORE_MORE
       local_get_block_density_nongas_label(log, &pread_file, PSTAR, strg); // TODO: Init density to what?  Is this needed at all?
#endif
       
       pread_done += pread_file;
       pread_todo -= pread_file;



     /* Move the particle pointers */  // TODO: -- delete because we are done reading?
      strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride*pread_file);
      strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride*pread_file);
      strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride*pread_file);
      strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride*pread_file);
      strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride*pread_file);
      strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride*pread_file);
      if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride*pread_file);
      if (strg.id.val     != NULL) strg.id.val     = (void *)(((char *)strg.id.val)     + strg.id.stride*pread_file);
      if (strg.u.val      != NULL) strg.u.val      = (void *)(((char *)strg.u.val)      + strg.u.stride*pread_file);
#ifdef STORE_MORE
      strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride*pread_file);
      strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride*pread_file);
#endif


	  }  // end of this task reads stars 


     //// TODO verify pread_done == pread
     fprintf(stderr,"rank %d of %d star pskip pread pread_done pread_todo pread_file pskip_file  %llu %llu %llu %llu %llu %llu\n", f->rank, f->size,pskip,pread,pread_done,pread_todo,pread_file,pskip_file); //////////////////////// TODO debugging


#if (defined WITH_MPI)
  {
    /* Reduce all values */
    double buffer[3];
    
    MPI_Allreduce(&(f->mmass),     buffer, 1, MPI_DOUBLE, MPI_MIN, f->mycomm);
    f->mmass = buffer[0];
    
    MPI_Allreduce(&(f->minweight), buffer, 1, MPI_DOUBLE, MPI_MIN, f->mycomm);
    f->minweight = buffer[0];
    
    MPI_Allreduce(&(f->maxweight), buffer, 1, MPI_DOUBLE, MPI_MAX, f->mycomm);
    f->maxweight = buffer[0];
    
    MPI_Allreduce(&(f->sumweight), buffer, 1, MPI_DOUBLE, MPI_SUM, f->mycomm);
    f->sumweight = buffer[0];
    
    /* note, this reduction is in fact repeated by io_gadget_scale_global() but what the heck... */
    MPI_Allreduce(&(f->maxpos),    buffer, 3, MPI_DOUBLE, MPI_MAX, f->mycomm);
    f->maxpos[0] = buffer[0];
    f->maxpos[1] = buffer[1];
    f->maxpos[2] = buffer[2];
    
    MPI_Allreduce(&(f->minpos),    buffer, 3, MPI_DOUBLE, MPI_MIN, f->mycomm);
    f->minpos[0] = buffer[0];
    f->minpos[1] = buffer[1];
    f->minpos[2] = buffer[2];
    
  }
#endif

  /* indicate that all files are ready for closure */ ///////////// TODO: why set to NULL before closing????????? is this bug in io_mgizmo.c ??? 
  //////////    for (i=0; i<f->numfiles; i++) {
  //////////      (f->files[i])->file = NULL;
  ///////////} 


  	/* Return the number of particles read */
	return pread_done;
}

extern uint64_t
io_nchilada_writepart(io_logging_t log,
                    io_nchilada_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg)
{
	return 0;
}

extern uint64_t
io_nchilada_writepart_ord(io_logging_t log,
                        io_nchilada_t f,
                        uint64_t pskip,
                        uint64_t pwrite,
                        void *nxt_part,
                        io_file_strg_struct_t strg)
{
	return UINT64_C(0);
}

/* same as io_tipsy_get ? */
extern bool
io_nchilada_get(io_logging_t log,
              io_nchilada_t f,
              io_file_get_t what,
              void *res)
{
  if ( (f == NULL) || (f->files[0]->header == NULL) )  // will use 1st NChilada file header
    return false;
  
  switch (what) {
    case IO_FILE_GET_NOPART_IN_FILE:  // TODO: verify not needed?
    case IO_FILE_GET_NOPART:
      *((long *)res) = (long)(f->no_part);  // used by io_file_get (in startrun.c). f->no_part is already set by io_nchilada_init()
      break;
    case IO_FILE_GET_NOVPART:
      *((double *)res) = (double)f->sumweight/f->mmass;
      break;
    case IO_FILE_GET_NOSPECIES: ////////////// TODO: verify for nchilada ??????????
      *((int *)res) = (int)f->no_species;
      break;
    case IO_FILE_GET_BOXSIZE:
      *((double *)res) = f->files[0]->header->boxsize;
      break;
  case IO_FILE_GET_PMASS:  //  this defines ahf mass unit  (m_fac   = simu.pmass;)
      *((double *)res) = f->files[0]->header->munit*f->mmass;
      break;
    case IO_FILE_GET_ZINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "zinitial is not set in a NCHILADA file, "
                      "using current redshift");
    case IO_FILE_GET_Z:
      *((double *)res) = 1./f->files[0]->header->time-1.;
      break;
    case IO_FILE_GET_AINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "ainitial is not set in a NCHILADA file, "
                      "using current expansion.");
    case IO_FILE_GET_A:
      *((double *)res) = f->files[0]->header->time;
      break;
    case IO_FILE_GET_OMEGA0:
      *((double *)res) = f->files[0]->header->omega0;
      break;
    case IO_FILE_GET_OMEGAL:
      *((double *)res) = f->files[0]->header->lambda0;
      break;
    case IO_FILE_GET_H:
      io_logging_warn(log, INT32_C(1),
                      "NCHILADA files don't store the timestep. "
                      "Setting to 100.0");
      *((double *)res) = H0;
      break;
    case IO_FILE_GET_DOUBLE:
      io_logging_warn(log, INT32_C(1),
                      "NCHILADA files don't store the use of "
                      "double precision. Assuming it is not "
                      "double precision.");
      *((int *)res) = 0;
      break;
    case IO_FILE_GET_MMASS:
        *((int *)res) = 1;
      break;
    case IO_FILE_GET_NOTSTEP:
      io_logging_warn(log, INT32_C(1),
                      "NCHILADA files don't store the timestep. "
                      "Setting to 500");
      *((int32_t *)res) = 500; //TODO: is nstep important for NCHILADA files?
      break;
    case IO_FILE_GET_TSTEP:
      io_logging_warn(log, INT32_C(1),
                      "NCHILADA files don't store the timestep. "
                      "Setting to 0.0");
      *((double *)res) = 0.0;
      break;
    case IO_FILE_GET_HEADERSTR:
      io_logging_warn(log, INT32_C(1),
                      "NCHILADA files don't have a header string. "
                      "Using a dummy one.");
      *((char **)res) = "No header string.";
      break;
    case IO_FILE_GET_MINWEIGHT:
      *((double *)res) = f->minweight/f->mmass;
      break;
    case IO_FILE_GET_MAXWEIGHT:
      *((double *)res) = f->maxweight/f->mmass; 
      break;
    default:
      io_logging_fatal(log, "Requesting something unkown in %s.",
                       __func__);
      return false;
  }
  
  return true;
}

///////// TODO: do we need to have io_nchilada_file_log?????????
extern void
io_nchilada_log(io_logging_t log, io_nchilada_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:             %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:             %s",
	               f->path);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:                 %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:             %" PRIi8,
	               f->swapped);
	io_logging_msg(log, INT32_C(5),
	               "  Number of files:          %" PRIi8,
	               f->numfiles);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles:        %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles w/mass: %" PRIu64,
	               f->no_part_with_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:            %" PRIi8,
	               f->multimass);
	io_logging_msg(log, INT32_C(5),
	               "  MMass (Halo parts):   %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight:       %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight:       %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of all weights:   %g",
	               f->sumweight);
	io_logging_msg(log, INT32_C(5),
	               "  No. of species:       %" PRIi32,
	               f->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Position scale:       %g",
	               f->posscale);
	io_logging_msg(log, INT32_C(5),
	               "  Weight scale:         %g",
	               f->weightscale);
	//////////	io_nchilada_header_log(log, f->header);   //TODO: put this in io_nchilada_file_log??

	return;
}

extern void
io_nchilada_resetscale(io_logging_t log,
                     io_nchilada_t f,
                     double posscale,
                     double weightscale) {
	if (f == NULL)
		return;

	io_logging_msg(log, INT32_C(8),
	               "Old posscale: %g   New posscale: %g",
	               f->posscale, posscale);
	io_logging_msg(log, INT32_C(8),
	               "Old weightscale: %g   New weightscale: %g",
	               f->weightscale, weightscale);
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}

/* same as io_tipsy_scale_particles ? */
extern uint64_t
io_nchilada_scale_particles(io_logging_t log,
                          double maxpos[],
                          double minpos[],
                          double boxsize,
                          double expansion,
                          double posscale,
                          double vunit,
                          double eunit,
                          double mmass,
                          uint64_t particles_read,
                          io_file_strg_struct_t strg)
{
	double box[3], shift[3], B;
	double scale_pos, scale_mom, scale_weight, scale_e, scale_rho;
	uint64_t i;
  
	/* Now we can do the scaling */
  B      = 1.0; // the standard szenario is that we do not want to scale the positions


  fprintf(stderr,"sssssssscaling %g %g %g %g %g %g %g %g %g %g %g %g box a posscale vunit eunit mmass posminmax",boxsize,expansion,posscale,vunit,eunit,mmass, maxpos[0],maxpos[1],maxpos[2],minpos[0],minpos[1],minpos[2]); fflush(stderr);/////////////////////// TODO: debugging 
  
	box[0] = fabs(maxpos[0] - minpos[0]);
	box[1] = fabs(maxpos[1] - minpos[1]);
	box[2] = fabs(maxpos[2] - minpos[2]);
	if (isgreater(box[0], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "x-Separation of particles exceeds boxsize "
		                "(%g > %g). Setting scale_pos accordingly!",
		                box[0], boxsize);
    B=MAX(B,box[0]);
	}
	if (isgreater(box[1], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "y-Separation of particles exceeds boxsize "
		                "(%g > %g). Setting scale_pos accordingly!",
		                box[1], boxsize);
    B=MAX(B,box[1]);
	}
	if (isgreater(box[2], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "z-Separation of particles exceeds boxsize "
		                "(%g > %g). Setting scale_pos accordingly!",
		                box[2], boxsize);
    B=MAX(B,box[2]);
	}
	io_logging_msg(log, INT32_C(4),
	               "Extreme positions: xmin = %g  xmax = %g",
	               minpos[0], maxpos[0]);
	io_logging_msg(log, INT32_C(4),
	               "                   ymin = %g  ymax = %g",
	               minpos[1], maxpos[1]);
	io_logging_msg(log, INT32_C(4),
	               "                   zmin = %g  zmax = %g",
	               minpos[2], maxpos[2]);
	shift[0] = (isless(minpos[0], 0.0) ? -minpos[0] : 0.0);
	shift[1] = (isless(minpos[1], 0.0) ? -minpos[1] : 0.0);
	shift[2] = (isless(minpos[2], 0.0) ? -minpos[2] : 0.0);
	io_logging_msg(log, INT32_C(4),
	               "Applying shift: (%g, %g, %g)",
	               shift[0], shift[1], shift[2]);

  
  
  
  
  
  
  
	/*=================================================================*\ 
    *                    NCHILADA UNIT SCALING                           *
    *             we only need to scale VELOCTIES!                    *
   \*=================================================================*/
#ifdef Hydra
  scale_pos      = 1.0/boxsize;
#else
  scale_pos      = 1.0/B;
#endif
  scale_mom      = expansion*expansion * vunit/(H0*boxsize);
  scale_weight   = 1.0/mmass;   // TODO: why scale by mmass (minimum dark matter particle mass)?  gadget does the same arbitrary re-scaling? 
  scale_e        = eunit;
  scale_rho      = scale_weight/(scale_pos*scale_pos*scale_pos);
  
	io_logging_msg(log, INT32_C(3),
	               "Scaling by:  positions:  %g", scale_pos);
	io_logging_msg(log, INT32_C(3),
	               "             velocities: %g", scale_mom);
	io_logging_msg(log, INT32_C(3),
	               "             weights:    %g", scale_weight);
	io_logging_msg(log, INT32_C(3),
	               "             energies:   %g", scale_e);
  
  
  
  /* keep track of the applied shifts and scales */
  simu.pos_shift[0] = shift[0];
  simu.pos_shift[1] = shift[1];
  simu.pos_shift[2] = shift[2];
  simu.pos_scale    = scale_pos;
  

  
  
  
  
  

	/* Define the actual scaling calls type independent */
#ifdef NCHILADA_ZOOMDATA
#	define SCALE_CALL(type) {\
    *((type *)strg.posx.val) *= (type)scale_pos; \
    *((type *)strg.posy.val) *= (type)scale_pos; \
    *((type *)strg.posz.val) *= (type)scale_pos; \
    *((type *)strg.posx.val)  = fmod(*((type *)strg.posx.val) + 1.5, 1.); \
    *((type *)strg.posy.val)  = fmod(*((type *)strg.posy.val) + 1.5, 1.); \
    *((type *)strg.posz.val)  = fmod(*((type *)strg.posz.val) + 1.5, 1.); \
    *((type *)strg.momx.val) *= (type)(scale_mom); \
    *((type *)strg.momy.val) *= (type)(scale_mom); \
    *((type *)strg.momz.val) *= (type)(scale_mom); \
    if (strg.u.val != NULL) {if(*((type *)strg.u.val) >= 0.) *((type *)strg.u.val) *= (type)(scale_e);} \
    if (strg.weight.val != NULL) *((type *)strg.weight.val) *= (type)(scale_weight); \
    if (strg.rho.val != NULL) *((type *)strg.rho.val) *= (type)(scale_rho); \
    if (strg.eps.val != NULL) *((type *)strg.eps.val) *= (type)(scale_pos); \
    strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride); \
    strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride); \
    strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride); \
    strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride); \
    strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride); \
    strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride); \
    if (strg.u.val != NULL) strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride); \
    if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride); \
    if (strg.rho.val != NULL) strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride); \
    if (strg.eps.val != NULL) strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride); \
}
#else
#	define SCALE_CALL(type) {\
		*((type *)strg.posx.val) += (type)(shift[0]); \
		*((type *)strg.posx.val) *= (type)(scale_pos); \
		*((type *)strg.posy.val) += (type)(shift[1]); \
		*((type *)strg.posy.val) *= (type)(scale_pos); \
		*((type *)strg.posz.val) += (type)(shift[2]); \
		*((type *)strg.posz.val) *= (type)(scale_pos); \
		*((type *)strg.momx.val) *= (type)(scale_mom); \
		*((type *)strg.momy.val) *= (type)(scale_mom); \
		*((type *)strg.momz.val) *= (type)(scale_mom); \
    if (strg.u.val != NULL) {if(*((type *)strg.u.val) >= 0.) *((type *)strg.u.val) *= (type)(scale_e);} \
		if (strg.weight.val != NULL) *((type *)strg.weight.val) *= (type)(scale_weight); \
    if (strg.rho.val != NULL) *((type *)strg.rho.val) *= (type)(scale_rho); \
    if (strg.eps.val != NULL) *((type *)strg.eps.val) *= (type)(scale_pos); \
    strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride); \
		strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride); \
		strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride); \
		strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride); \
		strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride); \
		strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride); \
    if (strg.u.val != NULL) strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride); \
		if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride); \
		if (strg.rho.val != NULL) strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride); \
		if (strg.eps.val != NULL) strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride); \
 }
#endif /* NCHILADA_ZOOMDATA */

	/* Do the scaling depending on the storage type */
	if (strg.bytes_float == sizeof(float)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(float);
		}
	} else if (strg.bytes_float == sizeof(double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(double);
		}
	} else if (strg.bytes_float == sizeof(long double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(long double);
		}
	} else {
		io_logging_fatal(log,
		                 "Don't know which floating point types "
		                 "has %" PRIi32 " bytes. Aborting read.",
		                 strg.bytes_float);
		return UINT64_C(0);
	}

	/* Clean the macro away */
#	undef SCALE_CALL
  
	/* And we are done */
	return particles_read;
}

/* same as io_tipsy_scale_global ? */
#ifdef WITH_MPI
extern void
io_nchilada_scale_global(io_logging_t log,
                       MPI_Comm comm,
                       double *maxpos,
                       double *minpos,
                       double *mmass,
                       double *minweight,
                       double *maxweight,
                       double *sumweight)
{
  int size, rank;
  double buffer[3];
  
  io_logging_msg(log, INT32_C(5),
                 "Updating local scale values to global values.");

  io_logging_msg(log, INT32_C(5),
                 "local : maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  MPI_Allreduce((void *)maxpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MAX, comm);
  maxpos[0] = buffer[0];
  maxpos[1] = buffer[1];
  maxpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  

  io_logging_msg(log, INT32_C(5),
                 "local : minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  MPI_Allreduce((void *)minpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MIN, comm);
  minpos[0] = buffer[0];
  minpos[1] = buffer[1];
  minpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  

  io_logging_msg(log, INT32_C(5), "local : sumweight = %g", *sumweight);
  MPI_Allreduce((void *)sumweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_SUM, comm);
  *sumweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: sumweight = %g", *sumweight);

  io_logging_msg(log, INT32_C(5), "local : minweight = %g", *minweight);
  MPI_Allreduce((void *)minweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MIN, comm);
  *minweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: minweight = %g", *minweight);
  
  io_logging_msg(log, INT32_C(5), "local : maxweight = %g", *maxweight);
  MPI_Allreduce((void *)maxweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MAX, comm);
  *maxweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: maxweight = %g", *maxweight);
  
  io_logging_msg(log, INT32_C(5), "local : mmass = %g", *mmass);
  MPI_Allreduce((void *)mmass, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MIN, comm);
  *mmass = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: mmass = %g", *mmass);
  
  return;
}
#endif /* WITH_MPI */


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

////////// TODO:  delete local functions that are not used.



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


// check for presence of particle family name (in file description.xml)
static int32_t
local_check_family(FILE *file,
                  const char *family)
{
  char line[MAXSTRING];    
  char pattern[MAXSTRING];
  //////////  char family_name[MAXSTRING];

  // the search pattern
  snprintf(pattern, sizeof(pattern), "<family name=\"%s\">", family);
  //////////////  snprintf(pattern, sizeof(pattern), "%s", family);

  // Read the file line by line and check whether each line contains family name
  rewind(file);  
  while (fgets(line, sizeof(line), file)) {
    ///////    fprintf(stderr,"line %s pattern %s\n", line, pattern);////////////////// TODO debugging
    if (strstr(line, pattern) != NULL) {
      fprintf(stderr,"found family: line %s pattern %s\n", line, pattern);////////////////// TODO debugging
      return 1; // pattern found. family present.
        }
    }
  return 0; // family not present in snapshot
}


static int32_t
local_i_am_gasfile(const char *filename)
{
  char pattern[MAXSTRING];

  snprintf(pattern, sizeof(pattern), "gas");
  if (strstr(filename, pattern) != NULL) {
    return 1;   
  }
  return 0;
}

static int32_t
local_i_am_darkfile(const char *filename)
{
  char pattern[MAXSTRING];

  snprintf(pattern, sizeof(pattern), "dark");
  if (strstr(filename, pattern) != NULL) {
    return 1;   
  }
  return 0;
}

static int32_t
local_i_am_starfile(const char *filename)
{
  char pattern[MAXSTRING];

  snprintf(pattern, sizeof(pattern), "star");
  if (strstr(filename, pattern) != NULL) {
    return 1;   
  }
  return 0;
}

// identify species from filename.  assumes NChilada naming conventions ///////// IS THIS NEEDED???
static int32_t
local_get_species_tag_from_filename(const char *filename)
{
  if (local_i_am_gasfile(filename)) {
    return 0;   
  }
  if (local_i_am_darkfile(filename)) {
    return -1;   
  }
  if (local_i_am_starfile(filename)) {
    return -4;   
  }
    
    return 1; // //////////// TODO: should produce an error somewhere
}


// identify species from filename.  assumes NChilada naming conventions ///////// IS THIS NEEDED??? /////////// TODO: delete this function  (instead use local_get_species_tag_from_filename)
static int32_t
local_get_species_from_filename(const char *filename)
{
  char pattern[MAXSTRING];
  
  snprintf(pattern, sizeof(pattern), "gas");
    if (strstr(filename, pattern) != NULL) {
      return 0;   
        }

  snprintf(pattern, sizeof(pattern), "dark");
    if (strstr(filename, pattern) != NULL) {
      return -1;   
        }

  snprintf(pattern, sizeof(pattern), "star");
    if (strstr(filename, pattern) != NULL) {
      return -4;   
        }
    
    return 1; // //////////// TODO: should produce an error somewhere
}


static int32_t
local_findfiles(io_logging_t log,
		  const char *path,
                  char ***fnames)
{
        int32_t numfiles, i;
        int32_t len;
	int32_t have_gas, have_dark, have_star;
	char filename[MAXSTRING];
	char familyname[MAXSTRING];
	char fullpath[MAXSTRING];
	FILE *fdescription;
	particle_field_enum field;


  //  Check for each Family type (gas,dark,star) in description.xml
  have_gas = 0;
  have_dark = 0;
  have_star = 0;

  /////////////////   strcpy(filename,"./description.xml");  // description.xml describes which files are present in snapshot
  snprintf(filename, sizeof(fullpath), "%s/description.xml",path);
  io_logging_msg(log, INT32_C(5),
		 "Will read nchilada file %s",filename);
  
  fdescription = fopen(filename, "r");
  if (fdescription == NULL) {
    io_logging_fatal(log,
		     "Could not open '%s' for reading.",
		     filename);
    return NULL;
  }
  strcpy(familyname,"gas");
  have_gas = local_check_family(fdescription, familyname);
  strcpy(familyname,"dark");
  have_dark = local_check_family(fdescription, familyname);
  strcpy(familyname,"star");
  have_star = local_check_family(fdescription, familyname);
  fclose(fdescription);
  
  io_logging_msg(log, INT32_C(5),
		 "presence of: gas %" PRIi32 "   dark %" PRIi32
		 "  star %" PRIi32 "  ",have_gas,have_dark,have_star);

  
  // number of NChilada files to be read
  // gas: mass,pos,vel,temperature,iord,(GasDensity,soft)
  // dark: mass,pos,vel,iord,(soft)
  // star: mass,pos,vel,iord,(soft)
  // !!!This should be set in a better way in case fields change in future
  numfiles = 5*have_gas + 4*have_dark + 4*have_star;
#ifdef STORE_MORE
  numfiles = 7*have_gas + 5*have_dark + 5*have_star;
#endif

  io_logging_msg(log, INT32_C(5),
		 "need to read %" PRIi32 " NChilada data files",
		  numfiles);

  *fnames = (char **)malloc(sizeof(char *) * numfiles);
  if (*fnames == NULL) {
    return INT32_C(0);
  }

  fprintf(stderr, "families fnames malloc done\n"); ///////////////////////////////// 

  
  ///////////////// TODO:  remove the complicated enum/field stuff because it is not used 
  /////////////     TODO:  verify the file names are present from description.xml instead of assuming?
  

  /* Will use only the filename later to identify which NChilada files are dark masses, star positions etc. 
   *   Will read pos, vel, temperature, etc. with a separate local_get_block function 
   *      that also fills in the AHF data structure for each field.
   */
  
  // Allocate memory for each filename string and assign the values
  // unsafely assuming filenames for each family -- would be better to read them from description.xml
  if (have_gas == 1) {
    field=GASMASS;    
    snprintf(familyname, sizeof(familyname), "gas");

    snprintf(filename, sizeof(filename), "mass");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }
    strcpy((*fnames)[field], fullpath);	

    fprintf(stderr, "field %d fullpath %s fnames  %s is %s\n",field,fullpath,(*fnames)[field],(*fnames)[0]);   /////////////////////////
    
    field=GASPOS;
    fprintf(stderr, "field pos is  %d \n",field);   /////////////////////////
    snprintf(filename, sizeof(filename), "pos");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////

    field=GASVEL;    
    snprintf(filename, sizeof(filename), "vel");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////

    field=GASTEMPERATURE;    
    snprintf(filename, sizeof(filename), "temperature");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
        fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////

    field=GASID;    
    snprintf(filename, sizeof(filename), "iord");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
        fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////
	
#ifdef STORE_MORE
    field=GASDENSITY;    
    snprintf(filename, sizeof(filename), "GasDensity");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////

    field=GASSOFT;    
    snprintf(filename, sizeof(filename), "soft");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
        fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   /////////////////////////
#endif  
  } // if have_gas
  


 if (have_dark == 1) {
    field=DARKMASS;    
    snprintf(familyname, sizeof(familyname), "dark");

    snprintf(filename, sizeof(filename), "mass");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);	     

    field=DARKPOS;    
    snprintf(filename, sizeof(filename), "pos");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);

    field=DARKVEL;    
    snprintf(filename, sizeof(filename), "vel");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    
    field=DARKID;    
    snprintf(filename, sizeof(filename), "iord");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   	

#ifdef STORE_MORE
    field=DARKSOFT;    
    snprintf(filename, sizeof(filename), "soft");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
#endif  
  } // if have_dark
  


  
  
 if (have_star == 1) {
    field=STARMASS;    
    snprintf(familyname, sizeof(familyname), "star");

    snprintf(filename, sizeof(filename), "mass");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);	     

    field=STARPOS;    
    snprintf(filename, sizeof(filename), "pos");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);

    field=STARVEL;    
    snprintf(filename, sizeof(filename), "vel");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);

    field=STARID;    
    snprintf(filename, sizeof(filename), "iord");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
    fprintf(stderr, "fullpath %s fnames  %s \n",fullpath,(*fnames)[field]);   
	    
#ifdef STORE_MORE
    field=STARSOFT;    
    snprintf(filename, sizeof(filename), "soft");
    snprintf(fullpath, sizeof(fullpath), "%s/%s/%s",path,familyname,filename);
    (*fnames)[field] = (char *)malloc((strlen(fullpath) + 1));
    if ((*fnames)[field] == NULL) {
      fprintf(stderr, "fnames allocation failed\n");
      return 1;
    }  
    strcpy((*fnames)[field], fullpath);
#endif  
  } // if have_star
  
 fprintf(stderr, "findfiles done %d files fullpath %s\n",numfiles,fullpath);   ///////////////////////// TODO debugging /////////


 /////////////////////// TODO debugging
 for (i=0; i<numfiles; i++)
   {
     fprintf(stderr, "xxxxxxxx file i %d fnames  %s \n",i,(*fnames)[i]);   /////////////////////////
   }
 //////////////////////////// TODO end debugging
 
  return numfiles;
  

}




///////////// TODO: make sure byte endian/magic number is checked.
///////////// TODO:    simplest (not cleanest) way is to read each field, 1 particle at a time. -- below is a copy of the gadget reader
///////////////////////* These local_get_block  functions could be combined into a single function that has some if statements to take care of strg.momx vs strg.posx vs strg.u etc.  */
static uint64_t
local_get_block_pos(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
        double fposx, fposy, fposz;
        float dummy;
        uint64_t i;

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = 3*bytes_file;
	blocksize = (uint64_t)partsize * (*pread);   ///// TODO: this blocksize calculation is a remnant of gadget readers and should be removed
	
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR);  // assumes that we are already at end of header

        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle positions "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored, which should equal %f",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.),
		       (float)(blocksize/1024./1024.));
        /* Checks whether the number of particles to read make sense */  // TODO: should have a more robust check? and in all block readers
        io_logging_msg(log, INT32_C(3),
                       "Asked to read %" PRIu64 " and to skip %" PRIu64
                       " particles. Checking those numbers.",
                       *pread, *pskip);
        if (*pskip > f->no_part) {
                *pskip = f->no_part;
                io_logging_msg(log, INT32_C(3),
                               "Cannot skip more than there is, will now "
                               "only skip %" PRIu64 " particles.",
                               *pskip);
        }
        if ( f->no_part - *pskip < *pread ) {
                *pread = f->no_part - *pskip;
                io_logging_msg(log, INT32_C(3),
                               "Cannot read more than there is left after "
                               "skipping. Will only read %" PRIu64
                               "particles.", *pread);
        }

	
        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);  // fseeko for "large" files. fseek has 32bit skip limit.
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR);  ////////// TODO: assumes that header has been read and at start of data
#endif
        /* Set extreme position detectors */
        f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
        f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;

        /* Loop over the particle positions */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data from the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fposx = (double)dummy;
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fposy = (double)dummy;
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fposz = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &fposx, f->swapped);
                        io_util_readdouble(f->file, &fposy, f->swapped);
                        io_util_readdouble(f->file, &fposz, f->swapped);
                }
                /* STEP 2:  Detect extreme positions */
                if (isless(fposx, f->minpos[0]))                        f->minpos[0] = fposx;
                if (isless(fposy, f->minpos[1]))                        f->minpos[1] = fposy;
                if (isless(fposz, f->minpos[2]))                        f->minpos[2] = fposz;
                if (isgreater(fposx, f->maxpos[0]))             f->maxpos[0] = fposx;
                if (isgreater(fposy, f->maxpos[1]))             f->maxpos[1] = fposy;
                if (isgreater(fposz, f->maxpos[2]))             f->maxpos[2] = fposz;
    
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.posx.val) = (float)fposx;
                        *((float *)strg.posy.val) = (float)fposy;
                        *((float *)strg.posz.val) = (float)fposz;
                } else {
                        *((double *)strg.posx.val) = fposx;
                        *((double *)strg.posy.val) = fposy;
                        *((double *)strg.posz.val) = fposz;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride);
                strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride);
                strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride);
    
        } /* End of particle position loop */

  
	/* Go to the end of the particle position block */ /////////// TODO: close file now? or after returning?
	/////////////        fseek(f->file, partsize*(f->no_part - (*pread + *pskip)), SEEK_CUR);
	//////////        SKIP2;
	////////////        CHECK_BLOCK;

        /* And return the number of read particles for error checking */
        return *pread;
}

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
        double fmomx, fmomy, fmomz;
        float dummy;
        uint64_t i;

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = 3*bytes_file;

	
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR);  /* fseek should be ok because not seeking far */
        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle velocities "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored.",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.));

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR);  
#endif
 
	
        /* Loop over the particle velocities */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data form the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fmomx = (double)dummy;
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fmomy = (double)dummy;
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fmomz = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &fmomx, f->swapped);
                        io_util_readdouble(f->file, &fmomy, f->swapped);
                        io_util_readdouble(f->file, &fmomz, f->swapped);
                }
    
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.momx.val) = (float)fmomx;
                        *((float *)strg.momy.val) = (float)fmomy;
                        *((float *)strg.momz.val) = (float)fmomz;
                } else {
                        *((double *)strg.momx.val) = fmomx;
                        *((double *)strg.momy.val) = fmomy;
                        *((double *)strg.momz.val) = fmomz;
                }

                /* STEP 4:  Increment the pointers to the next particle */
                strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride);
                strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride);
                strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride);		
		
        } /* End of particle velocity loop */

        /* And return the number of read particles for error checking */
        return *pread;
}

static uint64_t
local_get_block_mass(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
	double fweight, oldfweight;
        float dummy;
        uint64_t i;

/* We are going to need that a few times for book-keeping */
#	define BOOK_KEEPING {\
        if (   isgreater(fweight, oldfweight) || isless(fweight, oldfweight)) { \
            f->no_species++; \
            oldfweight = fweight; \
            if (isless(fweight,f->minweight)) f->minweight = fweight; \
            if (isgreater(fweight, f->maxweight)) f->maxweight = fweight; } \
	}
	/***** 	                        if (    (i >= f->header->np[0]) \
                             && (i < f->header->np[0]+f->header->np[1]) \
                             && isless(fweight, f->mmass) ) \
			     f->mmass = fweight; \  ***/         ///// mgadget bits  TODO: debugging delete these gadget lines //////////////

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = bytes_file; 
       
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR); 
        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle masses "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored.",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.));

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR); 
#endif

        /* Initialize some things */
        f->sumweight = 0.0;
        oldfweight = 0.0;
        f->no_species = 0;

        /* Loop over the particle masses */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data form the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fweight = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &fweight, f->swapped);
                }
    
                /* STEP 2:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.weight.val) = (float)fweight;
                } else {
                        *((double *)strg.weight.val) = fweight;
                }
		/* STEP 3:  Do the book-keeping */
                BOOK_KEEPING;
		f->sumweight += fweight;

                /* STEP 4:  Increment the pointers to the next particle */
                strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride);

        } /* End of particle mass loop */

	/* hack -- mmass seems to be defined in ahf as the minimum dark matter particle mass.  */
	/*         but we store the minimum mass for every file, then only use it later if dark file */
	
	f->mmass = f->minweight;  
	
#       undef BOOK_KEEPING
	
        /* And return the number of read particles for error checking */
        return *pread;
}

/////// TODO: set u to -1 for dark.  -4 for stars.
static uint64_t
local_get_block_u(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
        double fu;
        float dummy;
        uint64_t i;

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = bytes_file; 
       
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR); 
        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle energies "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored.",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.));

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR); 
#endif
 
        /* Loop over the particle masses */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data form the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        fu = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &fu, f->swapped);
                }
    
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.u.val) = (float)fu;
                } else {
                        *((double *)strg.u.val) = fu;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
    
        } /* End of particle velocity loop */

        /* And return the number of read particles for error checking */
        return *pread;
}




static uint64_t
local_get_block_density(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
        double frho;
        float dummy;
        uint64_t i;

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = bytes_file; 
       
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR); 
        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle densities "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored.",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.));

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR); 
#endif
 
        /* Loop over the particle densities */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data form the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        frho = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &frho, f->swapped);
                }
    
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.rho.val) = (float)frho;
                } else {
                        *((double *)strg.rho.val) = frho;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride);
    
        } /* End of particle density loop */

        /* And return the number of read particles for error checking */
        return *pread;
}



static uint64_t
local_get_block_softening(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_file;
        double feps;
        float dummy;
        uint64_t i;

        /* Figure out how many bytes are used for float storage */
	if (f->header->code == 9) { 
	  bytes_file = 4;
	}
	if (f->header->code == 10) { 
	  bytes_file = 8;
	}
        /* Set the particle size */
        partsize = bytes_file; 
       
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR);  
        CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
        io_logging_msg(log, INT32_C(1),
                       "A total of %" PRIu64 " particle softenings "
                       "with %" PRIi32 " bytes per float (%f MB total) "
                       "are stored.",
                       f->no_part, bytes_file,
                       (float)(f->no_part * partsize/1024./1024.));

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR); 
#endif
 
        /* Loop over the particle densities */
        for (i=0; i<*pread; i++) {
                /* STEP 1:  Read the particle data form the file */
                if (bytes_file == sizeof(float)) {
                        io_util_readfloat(f->file, &dummy, f->swapped);
                        feps = (double)dummy;
                } else if (bytes_file == sizeof(double)){
                        io_util_readdouble(f->file, &feps, f->swapped);
                }
    
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.eps.val) = (float)feps;
                } else {
                        *((double *)strg.eps.val) = feps;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.eps.val = (void *)(((char *)strg.eps.val) + strg.eps.stride);
    
        } /* End of particle density loop */

        /* And return the number of read particles for error checking */
        return *pread;
}

static uint64_t
local_get_block_id(io_logging_t log,
                    io_nchilada_file_t f, 
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
#ifdef FSEEKO
        int64_t blocksize, blocksize2, partsize;
#else
        long blocksize, blocksize2, partsize;
#endif
        uint32_t bytes_int_file;
        uint64_t fid;
        uint32_t dummy_int;
        uint64_t i;
	int64_t signed_fid;
	int32_t dummy_signed_int;
	
        /* Figure out how many bytes are used for float storage */
	// // TODO: assuming ids are signed ints given by "code" but ahf ids are unsigned.
	// //     TODO      This will break IDs if there are any coded negative IDs.
	// //     TODO      if Npart > 2**31 (32 bit) or Npart > 2**63 (64 bit)
	// // TODO: can nchilada IDs be unsigned int type?
	if (f->header->code == 5) { // int32 (signed)
	  bytes_int_file = 4;
	}
	else if (f->header->code == 7) { // int64
	  bytes_int_file = 8;
	}
	else {
	  io_logging_fatal(log,
			   "Can't handle reading ids "
			   "NChilada code %" PRIi32 " not understood.",
			   f->header->code);
	  return pread;
	}
        /* Set the particle size */
        partsize = bytes_int_file; 
       
	/* skip nchilada field min and max.  would be better to read and then verify later */
	fseek(f->file, partsize*2, SEEK_CUR);  

        /* Go to the first particle we want to read */
#ifdef FSEEKO
	fseeko(f->file, partsize*(*pskip), SEEK_CUR);
#else
        fseek(f->file, partsize*(*pskip), SEEK_CUR); 
#endif
 
        /* See if we have to read the IDs */
        if (strg.id.val == NULL) {
                io_logging_warn(log, INT32_C(1),
                            "Discarding IDs as no storage for the IDs "
                                "has been specified.");
#ifdef FSEEKO
                fseeko(f->file, partsize*(*pread), SEEK_CUR);
#else
                fseek(f->file, partsize*(*pread), SEEK_CUR);
#endif
		
        } else {
                /* Loop over the particle IDs */
                for (i=0; i<*pread; i++) {
                        /* STEP 1:  Read the ID from the file */
                        if (bytes_int_file == sizeof(int32_t)) {
                                io_util_readint32(f->file, &dummy_signed_int, f->swapped);
                                fid = (uint64_t)dummy_signed_int;
                        } else  {
                                io_util_readint64(f->file, &signed_fid, f->swapped);
				fid = (uint64_t)signed_fid;
                        }
                        /* STEP 2:  Store the ID in the array */
                        if (strg.bytes_int == 4) {
			        *((uint32_t *)strg.id.val) = (uint32_t)fid; 
                        } else {
			        *((uint64_t *)strg.id.val) = fid;
                        }
                        /* STEP 3:  Increment the pointers to the next particle */
                        strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride);
                } /* End of particle ID loop */
        } /* End of catch NULL-Id */

        /* And return the number of read particles for error checking */
        return *pread;
}

static uint64_t
local_get_block_u_nongas_label(io_logging_t log,
		      uint64_t *pread,
		      double tag,
		      io_file_strg_struct_t strg)
{
        uint64_t i;

        io_logging_msg(log, INT32_C(1),
                       "Set a total of %" PRIu64 " nongas particle energies tags ",
                       *pread);

        /* Loop over the particles */
        for (i=0; i<*pread; i++) {
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.u.val) = (float)tag;
                } else {
                        *((double *)strg.u.val) = tag;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
    
        } /* End of particle loop */

        /* And return the number of read particles for error checking */
        return *pread;
}


//  initialize the densities (dark and star densities not part of NChilada) // is this needed?
static uint64_t
local_get_block_density_nongas_label(io_logging_t log,
		      uint64_t *pread,
		      double tag,
		      io_file_strg_struct_t strg)
{
        uint64_t i;

        io_logging_msg(log, INT32_C(1),
                       "Set a total of %" PRIu64 " nongas particle densities tags",
                       *pread);

        /* Loop over the particles */
        for (i=0; i<*pread; i++) {
                /* STEP 3:  Store the particle in the array */
                if (strg.bytes_float == sizeof(float)) {
                        *((float *)strg.rho.val) = (float)tag;
                } else {
                        *((double *)strg.rho.val) = tag;
                }
                /* STEP 4:  Increment the pointers to the next particle */
                strg.rho.val = (void *)(((char *)strg.rho.val) + strg.rho.stride);
    
        } /* End of particle loop */

        /* And return the number of read particles for error checking */
        return *pread;
}

/* Getting rid of the macros */
#undef SKIP
#undef SKIP2
#undef CHECK_BLOCK
#undef DESCRIBE_BLOCK
#undef CHECK_FLOATBYTES

