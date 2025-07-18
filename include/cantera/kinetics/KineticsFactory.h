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

/**
 * Factory for kinetics managers.
 */
class KineticsFactory : public Factory<Kinetics>
{
public:
    static KineticsFactory* factory();

    void deleteFactory() override;

    /**
     * Return a new, empty kinetics manager.
     */
    Kinetics* newKinetics(const string& model);

private:
    static KineticsFactory* s_factory;
    KineticsFactory();
    static std::mutex kinetics_mutex;
};

//! @addtogroup kineticsmgr
//! @{

/**
 *  Create a new Kinetics instance.
 */
shared_ptr<Kinetics> newKinetics(const string& model);

//! Create a new kinetics manager, initialize it, and add reactions.
/*!
 * @param phases     Vector of phases containing species which participate in
 *     reactions, with the phase where the reactions occur (lowest-dimensional
 *     phase) listed first.
 * @param phaseNode  Phase entry for the phase where the reactions occur. This
 *     phase definition is used to determine the source of the reactions added
 *     to the Kinetics object.
 * @param rootNode   The root node of the file containing the phase definition,
 *     which will be treated as the default source for reactions
 * @param soln       The Solution object that this Kinetics object is being added to.
 */
shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode=AnyMap(),
                                 shared_ptr<Solution> soln={});

//! Create a new kinetics manager, initialize it, and add reactions.
/*!
 * @param phases      Vector of phases containing species which participate in
 *     reactions, with the phase where the reactions occur (lowest-dimensional
 *     phase) listed first.
 * @param filename    File containing the phase definition for the phase where
 *     the reactions occur. Searches the %Cantera data for this file.
 */
shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const string& filename);

/**
 * Add reactions to a Kinetics object.
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

//! @}

}

#endif
