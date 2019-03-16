/**
 *  @file ElectronFactory.cpp
 *     Definitions for the factory class that can create known Electron objects
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/ElectronFactory.h"
#include "cantera/electron/Electron.h"
#include "cantera/electron/WeakIonGasElectron.h"


using namespace std;

namespace Cantera
{

ElectronFactory* ElectronFactory::s_factory = 0;
std::mutex ElectronFactory::electron_mutex;

ElectronFactory::ElectronFactory()
{
    reg("WeaklyIonizedGas", []() { return new WeakIonGasElectron(); });
}

Electron* ElectronFactory::newElectron(const std::string& model)
{
    return create(model);
}

unique_ptr<Electron> newElectron(const AnyMap& rootNode, thermo_t* phase)
{
    unique_ptr<Electron> electron(newElectron(rootNode["electron"].asString()));
    electron->init(phase);
    addElectronCrossSections(*electron, rootNode["cross_section"]);
    return electron;
}

void addElectronCrossSections(Electron& electron, const AnyValue& cross_section)
{
    for (const auto& item : cross_section.asVector<AnyMap>()) {
        electron.addElectronCrossSection(newElectronCrossSection(item));
    }
}

}
