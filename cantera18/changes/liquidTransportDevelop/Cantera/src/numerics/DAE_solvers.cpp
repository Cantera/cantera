
#include "ct_defs.h"
#include "DAE_Solver.h"

// DAE_DEVEL is turned off at the current time
#define DAE_DEVEL
#ifdef DAE_DEVEL

#ifdef HAS_SUNDIALS
#include "IDA_Solver.cpp"
#endif

namespace Cantera {

  DAE_Solver* newDAE_Solver(string itype, ResidJacEval& f) {
        if (itype == "IDA") {
#ifdef HAS_SUNDIALS
            return new IDA_Solver(f);
#else
            raise CanteraError("newDAE_Solver","IDA solver requires sundials"
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
