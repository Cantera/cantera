#ifndef CANTERA_H_INCL
#define CANTERA_H_INCL

// Current 'Cantera.h' header 


// definitions
#ifndef CANTERA_APP
#define CANTERA_APP
#endif

namespace Cantera_CXX{ }

using namespace Cantera_CXX;

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

#endif
