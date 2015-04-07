/**
 *  @file SolidTransportData.cpp
 *  Source code for solid transport property evaluations.
 */

#include "cantera/transport/SolidTransportData.h"

using namespace std;

namespace Cantera
{
SolidTransportData::SolidTransportData() :
    speciesName("-"),
    ionConductivity(0),
    thermalConductivity(0),
    electConductivity(0),
    defectDiffusivity(0),
    defectActivity(0)
{

}

SolidTransportData::SolidTransportData(const SolidTransportData& right) :
    speciesName("-"),
    ionConductivity(0),
    thermalConductivity(0),
    electConductivity(0),
    defectDiffusivity(0),
    defectActivity(0)
{
    *this = right; //use assignment operator to do other work
}

SolidTransportData& SolidTransportData::operator=(const SolidTransportData& right)
{
    if (&right != this) {
        // These are all shallow pointer copies - yes, yes, yes horrible crime.
        speciesName        = right.speciesName;
        if (right.ionConductivity) {
            ionConductivity = (right.ionConductivity)->duplMyselfAsLTPspecies();
        }

        if (right.thermalConductivity) {
            thermalConductivity  = (right.thermalConductivity)->duplMyselfAsLTPspecies();
        }
        if (right.electConductivity) {
            electConductivity = (right.electConductivity)->duplMyselfAsLTPspecies();
        }
        if (right.defectDiffusivity) {
            defectDiffusivity = (right.defectDiffusivity)->duplMyselfAsLTPspecies();
        }
        if (right.defectActivity) {
            defectActivity = (right.defectActivity)->duplMyselfAsLTPspecies();
        }
    }
    return *this;
}

SolidTransportData::~SolidTransportData()
{
    delete ionConductivity;
    delete thermalConductivity;
    delete electConductivity;
    delete defectDiffusivity;
    delete defectActivity;
}

}
