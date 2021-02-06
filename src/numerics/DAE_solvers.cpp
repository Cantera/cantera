//! @file DAE_solvers.cpp Factory routine for picking the DAE solver package

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/DAE_Solver.h"
#include "cantera/numerics/IDA_Solver.h"

// DAE_DEVEL is turned off at the current time
#define DAE_DEVEL
#ifdef DAE_DEVEL

namespace Cantera
{
DAE_Solver* newDAE_Solver(const std::string& itype, ResidJacEval& f)
{
    if (itype == "IDA") {
        return new IDA_Solver(f);
    } else {
        throw CanteraError("newDAE_Solver",
                           "unknown DAE solver: "+itype);
    }
}
}

#endif
