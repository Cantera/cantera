//! @file ODE_integrators.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
