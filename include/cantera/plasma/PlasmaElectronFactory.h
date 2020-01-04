/**
 *  @file PlasmaElectronFactory.h
 *     Headers for the factory class that can create known PlasmaElectron objects
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef PLASMAELECTRON_FACTORY_H
#define PLASMAELECTRON_FACTORY_H

#include "cantera/plasma/PlasmaElectron.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{
//! Factory class for plasmaelectron data managers.
/*!
 * This class keeps a list of the known PlasmaElectron classes, and is
 * used to create new instances of these classes.
 */
class PlasmaElectronFactory : public Factory<PlasmaElectron>
{
public:
    //! Static function that creates a static instance of the factory.
    static PlasmaElectronFactory* factory() {
        std::unique_lock<std::mutex> lock(electron_mutex);
        if (!s_factory) {
            s_factory = new PlasmaElectronFactory;
        }
        return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(electron_mutex);
        delete s_factory;
        s_factory = 0;
    }

    virtual PlasmaElectron* newPlasmaElectron(const std::string& model);

private:
    //! static member of a single instance
    static PlasmaElectronFactory* s_factory;

    //! Private constructors prevents usage
    PlasmaElectronFactory();

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex electron_mutex;
};

//! Create a new electron manager instance.
/*!
 * @param model   String to look up the model against
 * @returns a pointer to a new PlasmaElectron instance matching the model string.
 *   Returns NULL if something went wrong. Throws an exception
 *   UnknownThermoPhaseModel if the string wasn't matched.
 */
inline PlasmaElectron* newPlasmaElectron(const std::string& model)
{
    return PlasmaElectronFactory::factory()->create(model);
}

unique_ptr<PlasmaElectron> newPlasmaElectron(const AnyMap& phaseNode=AnyMap(),
                                             const AnyMap& rootNode=AnyMap(),
                                             thermo_t* phase = 0);

void addElectronCrossSections(PlasmaElectron& electron, const AnyValue& cross_section, const AnyValue& names);

}

#endif
