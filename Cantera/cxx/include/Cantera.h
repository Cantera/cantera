/**
 * @file Cantera.h
 *  Basic include file to be used in all Cantera application environments.
 */

/* 
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

/*
 *  Note, this include should be the first include that code containing the
 *  Cantera namespace sees when in the Cantera application environment.
 */

#ifndef CANTERA_H_INCL
#define CANTERA_H_INCL

// If we are using this file, then we are in the Cantera Apps environment.
// Define a variable to signify this fact. 
#ifndef CANTERA_APP
#define CANTERA_APP
#endif

// define the presence of the Cantera_CXX namespace
namespace Cantera_CXX{ }

// Include global typedefs and values for physical constants using SI units
#include "kernel/ct_defs.h"

// some useful functions
#include "kernel/global.h"

// the CanteraError exception class
#include "kernel/ctexceptions.h"

// The Cantera logger class
#include "kernel/logger.h"

// Include the timer
#include "kernel/clockWC.h"

// Include routines for reading and writing XML files
#include "kernel/xml.h"

// Include string utility routines
#include "kernel/stringUtils.h"

// Include the array object
#include "kernel/Array.h"

#endif


