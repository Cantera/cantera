/**
 *  @file LiquidTransportData.cpp
 *  Source code for liquid transport property evaluations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/LiquidTransportData.h"
using namespace std;

namespace Cantera
{
LiquidTransportData::LiquidTransportData() :
    speciesName("-"),
    hydroRadius(0),
    viscosity(0),
    ionConductivity(0),
    thermalCond(0),
    electCond(0),
    speciesDiffusivity(0)
{
    warn_deprecated("class LiquidTransportData", "To be removed after Cantera 2.4");
}

LiquidTransportData::LiquidTransportData(const LiquidTransportData& right) :
    speciesName("-"),
    hydroRadius(0),
    viscosity(0),
    ionConductivity(0),
    thermalCond(0),
    electCond(0),
    speciesDiffusivity(0)
{
    *this = right; //use assignment operator to do other work
}

LiquidTransportData& LiquidTransportData::operator=(const LiquidTransportData& right)
{
    if (&right != this) {
        // These are all shallow pointer copies - yes, yes, yes horrible crime.
        speciesName = right.speciesName;
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
            thermalCond = (right.thermalCond)->duplMyselfAsLTPspecies();
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

LiquidTransportData::~LiquidTransportData()
{
    delete hydroRadius;
    delete viscosity;
    delete ionConductivity;

    for (size_t k = 0; k < mobilityRatio.size(); k++) {
        delete mobilityRatio[k];
    }
    for (size_t k = 0; k < selfDiffusion.size(); k++) {
        delete selfDiffusion[k];
    }

    delete thermalCond;
    delete electCond;
    delete speciesDiffusivity;
}

}
