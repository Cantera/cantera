/**
 *  @file KineticsFactory.cpp
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/KineticsFactory.h"

#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/kinetics/AqueousKinetics.h"
#include "cantera/base/xml.h"

using namespace std;

namespace Cantera
{

KineticsFactory* KineticsFactory::s_factory = 0;
mutex_t KineticsFactory::kinetics_mutex;

Kinetics* KineticsFactory::newKinetics(XML_Node& phaseData,
                                       vector<ThermoPhase*> th)
{
    /*
     * Look for a child of the XML element phase called
     * "kinetics". It has an attribute name "model".
     * Store the value of that attribute in the variable kintype
     */
    string kintype = phaseData.child("kinetics")["model"];

    // Create a kinetics object of the desired type
    Kinetics* k = newKinetics(kintype);
    // Now that we have the kinetics manager, we can
    // import the reaction mechanism into it.
    importKinetics(phaseData, th, k);

    // Return the pointer to the kinetics manager
    return k;
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{
    string lcmodel = lowercase(model);
    if (lcmodel == "none") {
        return new Kinetics();
    } else if (lcmodel == "gaskinetics") {
        return new GasKinetics();
    } else if (lcmodel == "interface") {
        return new InterfaceKinetics();
    } else if (lcmodel == "edge") {
        return new EdgeKinetics();
    } else if (lcmodel == "aqueouskinetics") {
        return new AqueousKinetics();
    } else {
        throw UnknownKineticsModel("KineticsFactory::newKinetics", model);
    }
}

}
