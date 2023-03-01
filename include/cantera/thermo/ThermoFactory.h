/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef THERMO_FACTORY_H
#define THERMO_FACTORY_H

#include "ThermoPhase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

/*!
 * @addtogroup thermoprops
 *
 * Standard ThermoPhase objects may be instantiated by calling the main %Cantera
 * factory class for ThermoPhase objects; This class is called ThermoFactory.
 */
//@{


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
    virtual void deleteFactory();

    //! Create a new thermodynamic property manager.
    /*!
     * @param model  The name of the thermo model
     * @returns a pointer to a new ThermoPhase object of the type specified. Throws a
     *     CanteraError if the named model isn't registered with ThermoFactory.
     * @deprecated  To be removed after Cantera 3.0; superseded by newThermo()
     */
    virtual ThermoPhase* newThermoPhase(const std::string& model);

private:
    //! static member of a single instance
    static ThermoFactory* s_factory;

    //! Private constructors prevents usage
    ThermoFactory();

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex thermo_mutex;
};

//! @copydoc ThermoFactory::newThermo(const string&)
//! @deprecated  To be removed after Cantera 3.0; superseded by newThermo()
ThermoPhase* newThermoPhase(const string& model);

//! Create a new ThermoPhase instance.
/*!
 * @param model   String to look up the model against
 * @returns a shared pointer to a new ThermoPhase object of the type specified. Throws a
 *     CanteraError if the named model is not registered with ThermoFactory.
 */
shared_ptr<ThermoPhase> newThermo(const string& model);

//! Create a new ThermoPhase object and initialize it
/*!
 * @param phaseNode  The node containing the phase definition (that is, thermo
 *     model, list of species, and initial state)
 * @param rootNode   The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 */
shared_ptr<ThermoPhase> newThermoPhase(const AnyMap& phaseNode,
                                       const AnyMap& rootNode=AnyMap());

//! @copydoc ThermoFactory::newThermoPhase(const AnyMap&, const AnyMap&)
//! @deprecated  To be removed after Cantera 3.0; superseded by newThermoPhase()
unique_ptr<ThermoPhase> newPhase(const AnyMap& phaseNode,
                                 const AnyMap& rootNode=AnyMap());

//! Create and Initialize a ThermoPhase object from an input file.
/*!
 * This function uses AnyMap::fromYamlFile() to read the input file, newThermoPhase()
 * to create an empty ThermoPhase of the appropriate type, and setupPhase() to
 * initialize the phase.
 *
 * @param infile name of the input file
 * @param id     name (id) of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @returns an initialized ThermoPhase object.
 */
shared_ptr<ThermoPhase> newThermoPhase(const string& infile, const string& id);

//! @copydoc ThermoFactory::newThermoPhase(const string&, const string&)
//! @deprecated  To be removed after Cantera 3.0; superseded by newThermoPhase()
ThermoPhase* newPhase(const std::string& infile, std::string id="");

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

//@}

}

#endif
