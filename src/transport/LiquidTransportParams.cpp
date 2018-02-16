/**
 *  @file LiquidTransportParams.cpp
 *  Source code for liquid mixture transport property evaluations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/LiquidTransportParams.h"
using namespace std;

namespace Cantera
{

LiquidTransportParams::LiquidTransportParams() :
    viscosity(0),
    ionConductivity(0),
    thermalCond(0),
    speciesDiffusivity(0),
    electCond(0),
    hydroRadius(0),
    model_viscosity(LTI_MODEL_NOTSET),
    model_speciesDiffusivity(LTI_MODEL_NOTSET),
    model_hydroradius(LTI_MODEL_NOTSET),
    compositionDepTypeDefault_(LTI_MODEL_NOTSET)
{
    warn_deprecated("Class LiquidTransportParams", "To be removed after Cantera 2.4");
}

LiquidTransportParams::~LiquidTransportParams()
{
    delete viscosity;
    delete ionConductivity;
    delete thermalCond;
    delete speciesDiffusivity;
    delete electCond;
    delete hydroRadius;
}

} //namespace Cantera
