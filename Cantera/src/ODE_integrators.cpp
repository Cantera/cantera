#include "ct_defs.h"
#include "Integrator.h"

#ifdef NO_SUNDIALS
#undef HAS_SUNDIALS
#endif

#ifdef HAS_SUNDIALS
#include "CVodesIntegrator.cpp"
#else
#include "CVode.cpp"
#endif

namespace Cantera {

    Integrator* newIntegrator(std::string itype) {
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
