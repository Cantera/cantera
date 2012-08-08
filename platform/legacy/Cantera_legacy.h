#ifndef CANTERA_H_INCL
#define CANTERA_H_INCL

// Legacy 'Cantera.h' header for applications that have not
// updated to the new header files.

#warning "Use of Cantera.h is deprecated. Please include required headers directly instead."

// definitions
#ifndef CANTERA_APP
#define CANTERA_APP
#endif

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

#endif
