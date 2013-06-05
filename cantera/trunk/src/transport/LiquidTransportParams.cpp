/**
 *  @file LiquidTransportParams.cpp
 *  Source code for liquid mixture transport property evaluations.
 */

#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/MargulesVPSSTP.h"
using namespace std;

namespace Cantera
{
//! Exception thrown if an error is encountered while reading the transport database.
class LTPmodelError : public CanteraError
{
public:
    explicit LTPmodelError(const std::string& msg) :
        CanteraError("LTPspecies", "error parsing transport data: " + msg + "\n") {
    }
};

LiquidTransportParams::LiquidTransportParams() :
    TransportParams(),
    LTData(0),
    viscosity(0),
    ionConductivity(0),
    mobilityRatio(0),
    selfDiffusion(0),
    thermalCond(0),
    speciesDiffusivity(0),
    electCond(0),
    hydroRadius(0),
    model_viscosity(LTI_MODEL_NOTSET),
    model_speciesDiffusivity(LTI_MODEL_NOTSET),
    model_hydroradius(LTI_MODEL_NOTSET)
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
    TransportParams(),
    LTData(0),
    viscosity(0),
    thermalCond(0),
    speciesDiffusivity(0),
    electCond(0),
    hydroRadius(0),
    model_viscosity(LTI_MODEL_NOTSET),
    model_speciesDiffusivity(LTI_MODEL_NOTSET),
    model_hydroradius(LTI_MODEL_NOTSET)
{
    operator=(right);
}

LiquidTransportParams&  LiquidTransportParams::operator=(const LiquidTransportParams& right)
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

    throw CanteraError("LiquidTransportParams(const LiquidTransportParams &right)", "not tested");

    return *this;
}

} //namespace Cantera
