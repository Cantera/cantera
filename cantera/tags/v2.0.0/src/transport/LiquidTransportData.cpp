/**
 *  @file LiquidTransportData.cpp
 *  Source code for liquid transport property evaluations.
 */

#include "cantera/transport/LiquidTransportData.h"
using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete (x); x = 0; }
#endif
namespace Cantera
{

//====================================================================================================================
LiquidTransportData::LiquidTransportData() :
    speciesName("-"),
    hydroRadius(0),
    viscosity(0),
    ionConductivity(0),
    mobilityRatio(0),
    selfDiffusion(0),
    thermalCond(0),
    electCond(0),
    speciesDiffusivity(0)
{

}
//====================================================================================================================
// Copy constructor
LiquidTransportData::LiquidTransportData(const LiquidTransportData& right) :
    speciesName("-"),
    hydroRadius(0),
    viscosity(0),
    ionConductivity(0),
    mobilityRatio(0),
    selfDiffusion(0),
    thermalCond(0),
    electCond(0),
    speciesDiffusivity(0)
{
    *this = right; //use assignment operator to do other work
}
//====================================================================================================================
// Assignment operator
LiquidTransportData& LiquidTransportData::operator=(const LiquidTransportData& right)
{
    if (&right != this) {
        // These are all shallow pointer copies - yes, yes, yes horrible crime.
        speciesName        = right.speciesName;
        if (right.hydroRadius) {
            hydroRadius = (right.hydroRadius)->duplMyselfAsLTPspecies();
        }
        if (right.viscosity) {
            viscosity = (right.viscosity)->duplMyselfAsLTPspecies();
        }
        if (right.ionConductivity) {
            ionConductivity = (right.ionConductivity)->duplMyselfAsLTPspecies();
        }

        mobilityRatio = right.mobilityRatio;
        for (size_t k = 0; k < mobilityRatio.size(); k++) {
            if (right.mobilityRatio[k]) {
                mobilityRatio[k] = (right.mobilityRatio[k])->duplMyselfAsLTPspecies();
            }
        }

        selfDiffusion = right.selfDiffusion;
        for (size_t k = 0; k < selfDiffusion.size(); k++) {
            if (right.selfDiffusion[k]) {
                selfDiffusion[k] = (right.selfDiffusion[k])->duplMyselfAsLTPspecies();
            }
        }

        if (right.thermalCond) {
            thermalCond  = (right.thermalCond)->duplMyselfAsLTPspecies();
        }
        if (right.electCond) {
            electCond = (right.electCond)->duplMyselfAsLTPspecies();
        }
        if (right.speciesDiffusivity) {
            speciesDiffusivity = (right.speciesDiffusivity)->duplMyselfAsLTPspecies();
        }
    }
    return *this;
}
//====================================================================================================================
LiquidTransportData::~LiquidTransportData()
{

    SAFE_DELETE(hydroRadius);
    SAFE_DELETE(viscosity);
    SAFE_DELETE(ionConductivity);

    for (size_t k = 0; k < mobilityRatio.size(); k++) {
        if (mobilityRatio[k]) {
            SAFE_DELETE(mobilityRatio[k]);
        }
    }
    for (size_t k = 0; k < selfDiffusion.size(); k++) {
        if (selfDiffusion[k]) {
            SAFE_DELETE(selfDiffusion[k]);
        }
    }

    SAFE_DELETE(thermalCond);
    SAFE_DELETE(electCond);
    SAFE_DELETE(speciesDiffusivity);

}
//====================================================================================================================
}
