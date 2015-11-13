//! @file ODE_integrators.cpp
#include "cantera/base/ct_defs.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/CVodesIntegrator.h"

namespace Cantera
{

Integrator* newIntegrator(const std::string& itype)
{
    if (itype == "CVODE") {
        return new CVodesIntegrator();
    } else {
        throw CanteraError("newIntegrator",
                           "unknown ODE integrator: "+itype);
    }
    return 0;
}

}
