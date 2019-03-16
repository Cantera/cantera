/**
 *  @file ElectronFactory.h
 *     Headers for the factory class that can create known Electron objects
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef ELECTRON_FACTORY_H
#define ELECTRON_FACTORY_H

#include "cantera/electron/Electron.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{
//! Factory class for electron data managers.
/*!
 * This class keeps a list of the known Electron classes, and is
 * used to create new instances of these classes.
 */
class ElectronFactory : public Factory<Electron>
{
public:
    //! Static function that creates a static instance of the factory.
    static ElectronFactory* factory() {
        std::unique_lock<std::mutex> lock(electron_mutex);
        if (!s_factory) {
            s_factory = new ElectronFactory;
        }
        return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(electron_mutex);
        delete s_factory;
        s_factory = 0;
    }

    virtual Electron* newElectron(const std::string& model);

private:
    //! static member of a single instance
    static ElectronFactory* s_factory;

    //! Private constructors prevents usage
    ElectronFactory();

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex electron_mutex;
};

//! Create a new electron manager instance.
/*!
 * @param model   String to look up the model against
 * @returns a pointer to a new Electron instance matching the model string.
 *   Returns NULL if something went wrong. Throws an exception
 *   UnknownThermoPhaseModel if the string wasn't matched.
 */
inline Electron* newElectron(const std::string& model)
{
    return ElectronFactory::factory()->create(model);
}

unique_ptr<Electron> newElectron(const AnyMap& rootNode=AnyMap(), thermo_t* phase = 0);

void addElectronCrossSections(Electron& electron, const AnyValue& cross_section);

}

#endif
