/**
 *  @file KineticsFactory.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef KINETICS_FACTORY_H
#define KINETICS_FACTORY_H

#include "Kinetics.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

//! @deprecated Unused. To be removed after Cantera 2.6.
class UnknownKineticsModel : public CanteraError
{
public:
    UnknownKineticsModel(const std::string& proc, const std::string& kineticsModel) :
        CanteraError(proc, "Specified Kinetics model "
                     + kineticsModel +
                     " does not match any known type.") {
        warn_deprecated("class UnknownKineticsModel",
            "Unused. To be removed after Cantera 2.6.");
    }
};

/**
 * Factory for kinetics managers.
 */
class KineticsFactory : public Factory<Kinetics>
{
public:
    static KineticsFactory* factory() {
        std::unique_lock<std::mutex> lock(kinetics_mutex);
        if (!s_factory) {
            s_factory = new KineticsFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(kinetics_mutex);
        delete s_factory;
        s_factory = 0;
    }

    /**
     * Return a new, empty kinetics manager.
     */
    virtual Kinetics* newKinetics(const std::string& model);

private:
    static KineticsFactory* s_factory;
    KineticsFactory();
    static std::mutex kinetics_mutex;
};

/**
 *  Create a new kinetics manager.
 */
inline Kinetics* newKineticsMgr(const std::string& model)
{
    return KineticsFactory::factory()->newKinetics(model);
}

/**
 *  Create a new Kinetics instance.
 */
inline shared_ptr<Kinetics> newKinetics(const std::string& model)
{
    shared_ptr<Kinetics> kin(KineticsFactory::factory()->newKinetics(model));
    return kin;
}

/*!
 * Create a new kinetics manager, initialize it, and add reactions
 *
 * @param phases     Vector of phases containing species which participate in
 *     reactions, with the phase where the reactions occur (lowest-dimensional
 *     phase) listed first.
 * @param phaseNode  Phase entry for the phase where the reactions occur. This
 *     phase definition is used to determine the source of the reactions added
 *     to the Kinetics object.
 * @param rootNode   The root node of the file containing the phase definition,
 *     which will be treated as the default source for reactions
 */
unique_ptr<Kinetics> newKinetics(const std::vector<ThermoPhase*>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode=AnyMap());

/*!
 * Create a new kinetics manager, initialize it, and add reactions
 *
 * @param phases      Vector of phases containing species which participate in
 *     reactions, with the phase where the reactions occur (lowest-dimensional
 *     phase) listed first.
 * @param filename    File containing the phase definition for the phase where
 *     the reactions occur. Searches the Cantera data for this file.
 * @param phase_name  The name of the reacting phase in the input file (that is, the
 *     name of the first phase in the `phases` vector)
 */
unique_ptr<Kinetics> newKinetics(const std::vector<ThermoPhase*>& phases,
                                 const std::string& filename,
                                 const std::string& phase_name);

/*!
 * Add reactions to a Kinetics object
 *
 * @param kin        The Kinetics object to be initialized
 * @param phaseNode  Phase entry for the phase where the reactions occur. This
 *     phase definition is used to determine the source of the reactions added
 *     to the Kinetics object.
 * @param rootNode   The root node of the file containing the phase definition,
 *     which will be treated as the default source for reactions
 */
void addReactions(Kinetics& kin, const AnyMap& phaseNode,
                  const AnyMap& rootNode=AnyMap());

}

#endif
