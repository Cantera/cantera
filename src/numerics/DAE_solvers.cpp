/**
 *  @file DAE_solvers.cpp
 *       Factory routine for picking the DAE solver package
 */
/*
 * $Revision: 725 $
 * $Date: 2011-05-16 18:45:08 -0600 (Mon, 16 May 2011) $
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#include "cantera/base/ct_defs.h"
#include "cantera/numerics/DAE_Solver.h"
#include "cantera/numerics/IDA_Solver.h"

// DAE_DEVEL is turned off at the current time
#define DAE_DEVEL
#ifdef DAE_DEVEL


namespace Cantera {

  DAE_Solver* newDAE_Solver(std::string itype, ResidJacEval& f) {
        if (itype == "IDA") {
#ifdef HAS_SUNDIALS
            return new IDA_Solver(f);
#else
            throw CanteraError("newDAE_Solver","IDA solver requires sundials"
			       " package, but Cantera was not built with sundials.");
#endif
        }
        else {
            throw CanteraError("newDAE_Solver",
			       "unknown DAE solver: "+itype);
        }
    }
}

#
#endif
