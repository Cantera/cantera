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
#include "cantera/base/xml.h"
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
                     " does not match any known type.") {}
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
     * @param model  String to look up the model against
     * @returns a pointer to a new ThermoPhase instance matching the model
     *   string. Returns NULL if something went wrong. Throws an exception
     *   UnknownThermoPhaseModel if the string wasn't matched.
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

//! Create a new thermo manager instance.
/*!
 * @param model   String to look up the model against
 * @returns a pointer to a new ThermoPhase instance matching the model string.
 *   Returns NULL if something went wrong. Throws an exception
 *   UnknownThermoPhaseModel if the string wasn't matched.
 */
inline ThermoPhase* newThermoPhase(const std::string& model)
{
    return ThermoFactory::factory()->create(model);
}

//! Create a new ThermoPhase object and initializes it according to the XML tree
/*!
 * This routine first looks up the identity of the model for the solution
 * thermodynamics in the model attribute of the thermo child of the XML phase
 * node. Then, it does a string lookup using Cantera's internal ThermoPhase
 * Factory routines on the model to figure out what ThermoPhase derived class
 * should be assigned. It creates a new instance of that class, and then calls
 * importPhase() to populate that class with the correct parameters from the
 * XML tree.
 *
 * @param phase XML_Node reference pointing to the phase XML element.
 * @return  A pointer to the completed and initialized ThermoPhase object.
 *
 * @ingroup inputfiles
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
ThermoPhase* newPhase(XML_Node& phase);

//! Create a new ThermoPhase object and initialize it
/*!
 * @param phaseNode  The node containing the phase definition (i.e. thermo
 *     model, list of species, and initial state)
 * @param rootNode   The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 */
unique_ptr<ThermoPhase> newPhase(AnyMap& phaseNode,
                                 const AnyMap& rootNode=AnyMap());

//! Create and Initialize a ThermoPhase object from an input file.
/*!
 * For YAML input files, this function uses AnyMap::fromYamlFile() to read the
 * input file, newThermoPhase() to create an empty ThermoPhase of the
 * appropriate type, and setupPhase() to initialize the phase.
 *
 * For CTI and XML input files, this function uses get_XML_File() to read the
 * input file and newPhase(XML_Node) to create and initialize the phase.
 *
 * @param infile name of the input file
 * @param id     name (id) of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @returns an initialized ThermoPhase object.
 */
ThermoPhase* newPhase(const std::string& infile, std::string id="");

//! Import a phase information into an empty ThermoPhase object
/*!
 * Here we read an XML description of the thermodynamic information for a phase.
 * At the end of this routine, the phase should be ready to be used within
 * applications. This routine contains some key routines that are used as pass
 * back routines so that the phase (and the contents of the XML file) may
 * contain variable parameterizations for the specification of the species
 * standard states, the equation of state, and the specification of other
 * nonidealities. Below, a description is presented of the main algorithm for
 * bringing up a ThermoPhase object, with care to present points where
 * customizations occur.
 *
 * Before invoking this routine, either the ThermoPhase Factory routines are
 * called or direct constructor routines are called that instantiate an
 * inherited ThermoPhase object. This object is input to this routine, and
 * therefore contains inherited routines that drive the customization of the
 * initialization process.
 *
 * At the start of the routine, we import descriptions of the elements that make
 * up the species in a phase.
 *
 * We call setParametersFromXML(eos) to read parameters about the thermo phase
 * before the species are read in.
 *
 * We call addElementsFromXML() to add elements into the description of the
 * phase.
 *
 * We create a new species thermo manager. Function 'newSpeciesThermoMgr' looks
 * at the species in the database to see what thermodynamic property
 * parameterizations are used, and selects a class that can handle the
 * parameterizations found.
 *
 * We import information about the species, including their reference state
 * thermodynamic polynomials. We then freeze the state of the species in the
 * element.
 *
 * Finally, we call initThermoXML(), a member function of the ThermoPhase
 * object, to "finish" the description. Now that the species are known,
 * additional information may be read in about the thermodynamics of the phase,
 * (e.g.,  virial coefficients, which are binary or ternary interaction
 * parameters between species).
 *
 * @param phase This object must be the phase node of a complete XML tree
 *              description of the phase, including all of the species data. In
 *              other words while "phase" must point to an XML phase object, it
 *              must have sibling nodes "speciesData" that describe the species
 *              in the phase.
 * @param th    Pointer to the ThermoPhase object which will handle the
 *              thermodynamics for this phase. We initialize part of the
 *              ThermoPhase object here, especially for those objects which are
 *              part of the Cantera Kernel.
 * @ingroup thermoprops
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
void importPhase(XML_Node& phase, ThermoPhase* th);

//! Initialize a ThermoPhase object
/*!
 *  @param phase      The ThermoPhase object to be initialized
 *  @param phaseNode  The node containing the phase definition (i.e. thermo
 *     model, list of species, and initial state)
 * @param rootNode    The root node of the tree containing the phase definition,
 *     which will be used as the default location from which to read species
 *     definitions.
 */
void setupPhase(ThermoPhase& phase, AnyMap& phaseNode,
                const AnyMap& rootNode=AnyMap());

//! Add the elements given in an XML_Node tree to the specified phase
//!
//! @deprecated The XML input format is deprecated and will be removed in
//!     Cantera 3.0.
void installElements(Phase& th, const XML_Node& phaseNode);

//!  Search an XML tree for species data.
/*!
 * This utility routine will search the XML tree for the species named by the
 * string, kname. It will return the XML_Node pointer to the species data for
 * that species. Failures of any kind return the null pointer.
 *
 * @param kname String containing the name of the species.
 * @param phaseSpeciesData   Pointer to the XML speciesData element
 *              containing the species data for that phase.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
const XML_Node* speciesXML_Node(const std::string& kname,
                                const XML_Node* phaseSpeciesData);

//@}

}

#endif
