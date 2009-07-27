/**
 * @file Cantera.h
 *  Basic include file to be used in all Cantera application
 *  environments.
 */

/* 
 * $Revision: 1.11 $
 * $Date: 2009/01/05 23:34:40 $
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

#include "kernel/ct_defs.h"

// some useful functions
#include "kernel/global.h"

// the CanteraError exception class
#include "kernel/ctexceptions.h"

//
//#include "kernel/importCTML.h"

// The Cantera logger class
#include "kernel/logger.h"

// Include the timer
#include "kernel/clockWC.h"

// Include routines for reading and writing XML files
#include "kernel/xml.h"

// Include string utility routines
#include "kernel/stringUtils.h"

#endif






