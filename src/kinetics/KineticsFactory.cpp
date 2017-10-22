/**
 *  @file KineticsFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
std::mutex KineticsFactory::kinetics_mutex;

Kinetics* KineticsFactory::newKinetics(XML_Node& phaseData,
                                       vector<ThermoPhase*> th)
{
    // Look for a child of the XML element phase called "kinetics". It has an
    // attribute name "model". Store the value of that attribute in the variable
    // kintype
    string kintype = phaseData.child("kinetics")["model"];

    // Create a kinetics object of the desired type
    Kinetics* k = newKinetics(kintype);
    // Now that we have the kinetics manager, we can import the reaction
    // mechanism into it.
    importKinetics(phaseData, th, k);

    // Return the pointer to the kinetics manager
    return k;
}

KineticsFactory::KineticsFactory() {
    reg("none", []() { return new Kinetics(); });
    reg("gaskinetics", []() { return new GasKinetics(); });
    reg("interface", []() { return new InterfaceKinetics(); });
    reg("edge", []() { return new EdgeKinetics(); });
    reg("aqueouskinetics", []() { return new AqueousKinetics(); });
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{
    return create(toLowerCopy(model));
}

}
