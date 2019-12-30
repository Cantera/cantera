/**
 *  @file PlasmaElectronFactory.cpp
 *     Definitions for the factory class that can create known PlasmaElectron objects
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/PlasmaElectronFactory.h"
#include "cantera/plasma/PlasmaElectron.h"
#include "cantera/plasma/WeakIonGasElectron.h"


using namespace std;

namespace Cantera
{

PlasmaElectronFactory* PlasmaElectronFactory::s_factory = 0;
std::mutex PlasmaElectronFactory::electron_mutex;

PlasmaElectronFactory::PlasmaElectronFactory()
{
    reg("WeaklyIonizedGas", []() { return new WeakIonGasElectron(); });
}

PlasmaElectron* PlasmaElectronFactory::newPlasmaElectron(const std::string& model)
{
    return create(model);
}

unique_ptr<PlasmaElectron> newPlasmaElectron(const AnyMap& rootNode, thermo_t* phase)
{
    unique_ptr<PlasmaElectron> electron(newPlasmaElectron(rootNode["electron"].asString()));
    addElectronCrossSections(*electron, rootNode["cross_section"]);
    electron->init(phase);
    return electron;
}

void addElectronCrossSections(PlasmaElectron& electron, const AnyValue& cross_section)
{
    for (const auto& item : cross_section.asVector<AnyMap>()) {
        electron.addElectronCrossSection(newElectronCrossSection(item));
    }
}

}
