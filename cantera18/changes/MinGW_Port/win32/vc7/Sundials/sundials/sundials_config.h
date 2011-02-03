/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006/11/13 19:46:18 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 *------------------------------------------------------------------
 * SUNDIALS configuration header file
 *------------------------------------------------------------------
 */

/* Define SUNDIALS version number */
#define SUNDIALS_PACKAGE_VERSION "2.2.0"

/* FCMIX: Define Fortran name-mangling macro */
#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _

/* FCMIX: Define case of function names */


/* FCMIX: Define number of underscores to append to function names */


/* Define precision of SUNDIALS data type 'realtype' */
#define SUNDIALS_DOUBLE_PRECISION 1

/* Use generic math functions */
#define SUNDIALS_USE_GENERIC_MATH 1

/* FNVECTOR: Allow user to specify different MPI communicator */
#define SUNDIALS_MPI_COMM_F2C 0
