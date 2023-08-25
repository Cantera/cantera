/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known ThermoPhase objects
 *     (see @ref thermoprops and class @link Cantera::ThermoFactory ThermoFactory@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef THERMO_FACTORY_H
#define THERMO_FACTORY_H

#include "ThermoPhase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

/**
 * @addtogroup thermoprops
 */
//! @{


//! Factory class for thermodynamic property managers.
/*!
 * This class keeps a list of the known ThermoPhase classes, and is
 * used to create new instances of these classes.
 */
class ThermoFactory : public Factory<ThermoPhase>
{
public:
    //! Static function that creates a static instance of the factory.
    static ThermoFactory* factory();

    //! delete the static instance of this factory
    void deleteFactory() override;

private:
    //! static member of a single instance
    static ThermoFactory* s_factory;

    //! Private constructors prevents usage
    ThermoFactory();

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex thermo_mutex;
};

//! Create a new ThermoPhase instance.
 /*!
  * @param model   String to look up the model against
  * @returns a shared pointer to a new ThermoPhase object of the type specified. Throws a
  *     CanteraError if the named model is not registered with ThermoFactory.
  * @since New in %Cantera 3.0. Replaces newThermo
  */
 shared_ptr<ThermoPhase> newThermoModel(const string& model);

//! Create a new ThermoPhase object and initialize it
/*!
 * @param phaseNode  The node containing the phase definition (that is, thermo
 *     model, list of species, and initial state)
 * @param rootNode   The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 * @since New in %Cantera 3.0.
 */
shared_ptr<ThermoPhase> newThermo(const AnyMap& phaseNode,
                                  const AnyMap& rootNode=AnyMap());

//! Create and Initialize a ThermoPhase object from an input file.
/*!
 * This function uses AnyMap::fromYamlFile() to read the input file, newThermo()
 * to create an empty ThermoPhase of the appropriate type, and setupPhase() to
 * initialize the phase.
 *
 * @param infile name of the input file
 * @param id     name (id) of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @returns an initialized ThermoPhase object.
 * @since Changed in %Cantera 3.0. Prior to %Cantera 3.0, the function used a single
 *      argument specifying the thermo model, which is now deprecated.
 */
shared_ptr<ThermoPhase> newThermo(const string& infile, const string& id="");

//! Initialize a ThermoPhase object
/*!
 *  @param phase      The ThermoPhase object to be initialized
 *  @param phaseNode  The node containing the phase definition (that is, thermo
 *     model, list of species, and initial state)
 * @param rootNode    The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 */
void setupPhase(ThermoPhase& phase, const AnyMap& phaseNode,
                const AnyMap& rootNode=AnyMap());

//! @}

}

#endif
