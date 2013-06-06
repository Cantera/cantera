//! @file ODE_integrators.cpp
#include "cantera/base/ct_defs.h"
#include "cantera/numerics/Integrator.h"


#ifdef HAS_SUNDIALS
#include "cantera/numerics/CVodesIntegrator.h"
#else
#include "CVodeInt.h"
#endif

namespace Cantera
{

Integrator* newIntegrator(const std::string& itype)
{
    if (itype == "CVODE") {
#ifdef HAS_SUNDIALS
        return new CVodesIntegrator();
#else
        return new CVodeInt();
#endif
    } else {
        throw CanteraError("newIntegrator",
                           "unknown ODE integrator: "+itype);
    }
    return 0;
}

}
