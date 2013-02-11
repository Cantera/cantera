/**
 *  @file SolidTransportData.cpp
 *  Source code for solid transport property evaluations.
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2010-07-13 13:22:30 -0600 (Tue, 13 Jul 2010) $
 *  $Revision: 507 $
 */

#include "cantera/transport/SolidTransportData.h"

using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete (x); x = 0; }
#endif
namespace Cantera
{

//====================================================================================================================
SolidTransportData::SolidTransportData() :
    speciesName("-"),
    ionConductivity(0),
    thermalConductivity(0),
    electConductivity(0),
    defectDiffusivity(0),
    defectActivity(0)
{

}
//====================================================================================================================
// Copy constructor
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
//====================================================================================================================
// Assignment operator
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
//====================================================================================================================
SolidTransportData::~SolidTransportData()
{

    SAFE_DELETE(ionConductivity);
    SAFE_DELETE(thermalConductivity);
    SAFE_DELETE(electConductivity);
    SAFE_DELETE(defectDiffusivity);
    SAFE_DELETE(defectActivity);

}
//====================================================================================================================
}
