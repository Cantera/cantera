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

    virtual void deleteFactory();

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
 *  @deprecated  To be removed after Cantera 3.0; superseded by newKinetics.
 */
Kinetics* newKineticsMgr(const string& model);

/**
 *  Create a new Kinetics instance.
 */
shared_ptr<Kinetics> newKinetics(const string& model);

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
 * @param soln       The Solution object that this Kinetics object is being added to.
 */
shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode=AnyMap(),
                                 shared_ptr<Solution> soln={});

//! @see newKinetics(const vector<shared_ptr<ThermoPhase>>&, const AnyMap&, const AnyMap&, shared_ptr<Solution>)
//! @deprecated  To be removed after Cantera 3.0;
//!     superseded by newKinetics() returning shared_ptr
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
 * @deprecated The 'phase_name' argument is deprecated and will be removed after
 *     Cantera 3.0.
 * @since  Starting with Cantera 3.0, if the reacting phase is not the first item in the
 *     `phases` vector, a deprecation warning will be issued. In Cantera 3.1, this
 *     warning will become an error.
 */
shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const string& filename,
                                 const string& phase_name="");

//! @copydoc newKinetics(const vector<shared_ptr<ThermoPhase>>&, const string&, const string&)
//! @deprecated  To be removed after Cantera 3.0;
//!     superseded by newKinetics() returning shared_ptr
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
