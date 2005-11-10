#include "ct_defs.h"
#include "Integrator.h"

#ifdef HAS_SUNDIALS
#include "CVodesIntegrator.cpp"
#else
#include "CVode.cpp"
#endif

namespace Cantera {

    Integrator* newIntegrator(string itype) {
        if (itype == "CVODE") {
#ifdef HAS_SUNDIALS
            return new CVodesIntegrator();
#else
            return new CVodeInt();
#endif
        }
        else {
            throw CanteraError("newIntegrator",
                "unknown ODE integrator: "+itype);
        }
    }
}
