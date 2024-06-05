//! @file Interface.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Interface.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"

namespace Cantera
{

void Interface::setThermo(shared_ptr<ThermoPhase> thermo) {
    Solution::setThermo(thermo);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(thermo);
    if (!surf) {
        throw CanteraError("Interface::setThermo",
            "Thermo object of type '{}' does not descend from SurfPhase.",
            thermo->type());
    }
    m_surf = surf;
}

void Interface::setKinetics(shared_ptr<Kinetics> kinetics) {
    Solution::setKinetics(kinetics);
    auto surfkin = std::dynamic_pointer_cast<InterfaceKinetics>(kinetics);
    if (!surfkin) {
        throw CanteraError("Interface::setKinetics",
            "Kinetics object of type '{}' does not descend from InterfaceKinetics.",
            kinetics->kineticsType());
    }
    m_surfkin = surfkin;
}

shared_ptr<Interface> newInterface(const string& infile,
    const string& name, const vector<string>& adjacent)
{
    auto sol = newSolution(infile, name, "default", adjacent);
    auto iface = std::dynamic_pointer_cast<Interface>(sol);
    if (!iface) {
        auto rootNode = AnyMap::fromYamlFile(infile);
        AnyMap& phaseNode = rootNode["phases"].getMapWhere("name", name);
        throw InputFileError("newInterface", phaseNode,
            "Phase definition does not define a surface phase");
    }
    return iface;
}

shared_ptr<Interface> newInterface(const string& infile,
    const string& name, const vector<shared_ptr<Solution>>& adjacent)
{
    auto rootNode = AnyMap::fromYamlFile(infile);
    AnyMap& phaseNode = rootNode["phases"].getMapWhere("name", name);
    return newInterface(phaseNode, rootNode, adjacent);
}

shared_ptr<Interface> newInterface(AnyMap& phaseNode, const AnyMap& rootNode,
    const vector<shared_ptr<Solution>>& adjacent)
{
    auto sol = newSolution(phaseNode, rootNode, "default", adjacent);
    auto iface = std::dynamic_pointer_cast<Interface>(sol);
    if (!iface) {
        throw InputFileError("newInterface", phaseNode,
            "Phase definition does not define a surface phase");
    }
    return iface;
}

} // end namespace Cantera
