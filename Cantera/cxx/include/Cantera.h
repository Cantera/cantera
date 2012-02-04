/**
 * @file Cantera.h
 *  Basic include file to be used in all Cantera application
 *  environments.
 */

/* 
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CANTERA_H_INCL
#define CANTERA_H_INCL

// definitions
#ifndef CANTERA_APP
#define CANTERA_APP
#endif

namespace Cantera_CXX{ }

using namespace Cantera_CXX;

#include "ct_defs.h"

// some useful functions
#include "global.h"

// the CanteraError exception class
#include "ctexceptions.h"

//
//#include "kernel/importCTML.h"

// The Cantera logger class
#include "logger.h"

// Include the timer
#include "clockWC.h"

// Include routines for reading and writing XML files
#include "xml.h"

// Include string utility routines
#include "stringUtils.h"

#endif






