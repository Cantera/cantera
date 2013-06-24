/**
 * @file Cantera.h
 *  Basic include file to be used in all Cantera application environments.
 *  @deprecated To be removed in Cantera 2.2. Applications should include
 *     relevant headers directly.
 */

/* 
 * $Revision: 923 $
 * $Date: 2012-01-03 10:05:28 -0700 (Tue, 03 Jan 2012) $
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
#include "base/ct_defs.h"

// some useful functions
#include "base/global.h"

// the CanteraError exception class
#include "base/ctexceptions.h"

// The Cantera logger class
#include "base/logger.h"

// Include the timer
#include "base/clockWC.h"

// Include routines for reading and writing XML files
#include "base/xml.h"

// Include string utility routines
#include "base/stringUtils.h"

// Include the array object
#include "base/Array.h"

#endif


