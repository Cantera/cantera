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

//! Specific error to be thrown if the type of Thermo manager is unrecognized.
/*!
 * This particular error class may be caught, if the application may have other
 * models that the main Cantera application doesn't know about.
 *
 * @deprecated Unused. To be removed after Cantera 2.6.
 */
class UnknownThermoPhaseModel : public CanteraError
{
public:
    //! Constructor
    /*!
     * @param proc Function name where the error occurred.
     * @param thermoModel Sting name of ThermoPhase which didn't match
     */
    UnknownThermoPhaseModel(const std::string& proc,
                            const std::string& thermoModel) :
        CanteraError(proc, "Specified ThermoPhase model "
                     + thermoModel +
                     " does not match any known type.") {
        warn_deprecated("class UnknownThermoPhaseModel",
            "Unused. To be removed after Cantera 2.6.");
    }
};


//! Factory class for thermodynamic property managers.
/*!
 * This class keeps a list of the known ThermoPhase classes, and is
 * used to create new instances of these classes.
 */
class ThermoFactory : public Factory<ThermoPhase>
{
public:
    //! Static function that creates a static instance of the factory.
    static ThermoFactory* factory() {
        std::unique_lock<std::mutex> lock(thermo_mutex);
        if (!s_factory) {
            s_factory = new ThermoFactory;
        }
        return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(thermo_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Create a new thermodynamic property manager.
    /*!
     * @param model  The name of the thermo model
     * @returns a pointer to a new ThermoPhase object of the type specified. Throws a
     *     CanteraError if the named model isn't registered with ThermoFactory.
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

//! @copydoc ThermoFactory::newThermoPhase
inline ThermoPhase* newThermoPhase(const std::string& model)
{
    return ThermoFactory::factory()->create(model);
}

//! Create a new ThermoPhase instance.
/*!
 * @param model   String to look up the model against
 * @returns a shared pointer to a new ThermoPhase instance matching the model string.
 */
inline shared_ptr<ThermoPhase> newThermo(const std::string& model)
{
    ThermoPhase* tptr = ThermoFactory::factory()->create(model);
    return shared_ptr<ThermoPhase> (tptr);
}

//! Create a new ThermoPhase object and initialize it
/*!
 * @param phaseNode  The node containing the phase definition (that is, thermo
 *     model, list of species, and initial state)
 * @param rootNode   The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 */
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
