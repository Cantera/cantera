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
#include "cantera/kinetics/ElectrodeKinetics.h"
#include "cantera/base/xml.h"

using namespace std;

namespace Cantera
{

KineticsFactory* KineticsFactory::s_factory = 0;
mutex_t KineticsFactory::kinetics_mutex;

static int ntypes = 6;
static string _types[] = {"none", "GasKinetics", "Interface", "Edge", "AqueousKinetics", "ElectrodeKinetics"};
static int _itypes[]   = {0, cGasKinetics, cInterfaceKinetics, cEdgeKinetics, cAqueousKinetics, cElectrodeKinetics};

Kinetics* KineticsFactory::newKinetics(XML_Node& phaseData,
                                       vector<ThermoPhase*> th)
{
    /*
     * Look for a child of the xml element phase called
     * "kinetics". It has an attribute name "model".
     * Store the value of that attribute in the variable kintype
     */
    string kintype = phaseData.child("kinetics")["model"];
    /*
     * look up the string kintype in the list of known
     * kinetics managers (list is kept at the top of this file).
     * Translate it to an integer value, ikin.
     */
    int ikin=-1;
    int n;
    for (n = 0; n < ntypes; n++) {
        if (kintype == _types[n]) {
            ikin = _itypes[n];
        }
    }
    /*
     * Assign the kinetics manager based on the value of ikin.
     * Kinetics managers are classes derived from the base
     * Kinetics class. Unknown kinetics managers will throw a
     * CanteraError here.
     */
    Kinetics* k=0;
    switch (ikin) {

    case 0:
        k = new Kinetics();
        break;

    case cGasKinetics:
        k = new GasKinetics();
        break;

    case cInterfaceKinetics:
        k = new InterfaceKinetics();
        break;

    case cEdgeKinetics:
        k = new EdgeKinetics();
        break;

    case cAqueousKinetics:
        k = new AqueousKinetics();
        break;

     case cElectrodeKinetics:
        k = new ElectrodeKinetics();
        break;

    default:
        throw UnknownKineticsModel("KineticsFactory::newKinetics",
                                   kintype);
    }

    // Now that we have the kinetics manager, we can
    // import the reaction mechanism into it.
    importKinetics(phaseData, th, k);

    // Return the pointer to the kinetics manager
    return k;
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{

    int ikin = -1;
    int n;
    for (n = 0; n < ntypes; n++) {
        if (model == _types[n]) {
            ikin = _itypes[n];
        }
    }
    Kinetics* k=0;
    switch (ikin) {

    case cGasKinetics:
        k = new GasKinetics();
        break;

    case cInterfaceKinetics:
        k = new InterfaceKinetics();
        break;

    case cElectrodeKinetics:
        k = new ElectrodeKinetics();

    case cAqueousKinetics:
	k = new AqueousKinetics();

    default:
        throw UnknownKineticsModel("KineticsFactory::newKinetics",
                                   model);
    }
    return k;
}

}
