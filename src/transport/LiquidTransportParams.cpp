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

LiquidTransportParams::LiquidTransportParams(const LiquidTransportParams& right) :
    viscosity(0),
    thermalCond(0),
    speciesDiffusivity(0),
    electCond(0),
    hydroRadius(0),
    model_viscosity(LTI_MODEL_NOTSET),
    model_speciesDiffusivity(LTI_MODEL_NOTSET),
    model_hydroradius(LTI_MODEL_NOTSET),
    compositionDepTypeDefault_(LTI_MODEL_NOTSET)
{
    operator=(right);
}

LiquidTransportParams& LiquidTransportParams::operator=(const LiquidTransportParams& right)
{
    if (&right != this) {
        return *this;
    }

    LTData = right.LTData;

    delete viscosity;
    if (right.viscosity) {
        viscosity = new LiquidTranInteraction(*(right.viscosity));
    }
    delete ionConductivity;
    if (right.ionConductivity) {
        ionConductivity = new LiquidTranInteraction(*(right.ionConductivity));
    }
    deepStdVectorPointerCopy<LiquidTranInteraction>(right.mobilityRatio, mobilityRatio);
    deepStdVectorPointerCopy<LiquidTranInteraction>(right.selfDiffusion, selfDiffusion);

    delete thermalCond;
    if (right.thermalCond) {
        thermalCond = new LiquidTranInteraction(*(right.thermalCond));
    }
    delete speciesDiffusivity;
    if (right.speciesDiffusivity) {
        speciesDiffusivity = new LiquidTranInteraction(*(right.speciesDiffusivity));
    }

    delete electCond;
    if (right.electCond) {
        electCond = new LiquidTranInteraction(*(right.electCond));
    }
    delete hydroRadius;
    if (right.hydroRadius) {
        hydroRadius = new LiquidTranInteraction(*(right.hydroRadius));
    }
    model_viscosity = right.model_viscosity;
    model_ionConductivity = right.model_ionConductivity;
    deepStdVectorPointerCopy<LiquidTranMixingModel>(right.model_mobilityRatio, model_mobilityRatio);
    deepStdVectorPointerCopy<LiquidTranMixingModel>(right.model_selfDiffusion, model_selfDiffusion);
    thermalCond_Aij = right.thermalCond_Aij;
    model_speciesDiffusivity = right.model_speciesDiffusivity;
    diff_Dij = right.diff_Dij;
    model_hydroradius = right.model_hydroradius;
    radius_Aij = right.radius_Aij;
    compositionDepTypeDefault_ = right.compositionDepTypeDefault_;

    throw CanteraError("LiquidTransportParams(const LiquidTransportParams &right)", "not tested");
    return *this;
}

} //namespace Cantera
