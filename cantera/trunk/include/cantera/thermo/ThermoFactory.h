/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 */
// Copyright 2001  California Institute of Technology


#ifndef THERMO_FACTORY_H
#define THERMO_FACTORY_H

#include "ThermoPhase.h"
#include "cantera/base/xml.h"
#include "cantera/base/ct_thread.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

class SpeciesThermoFactory;
class VPSSMgr;

/*!
 *  @addtogroup thermoprops
 *
 *  Standard %ThermoPhase objects may be instantiated by calling
 *  the main %Cantera factory class for %ThermoPhase objects; This class is called ThermoFactory.
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
class ThermoFactory : public FactoryBase
{

public:

    //! Static function that creates a static instance of the factory.
    static ThermoFactory* factory() {
        ScopedLock lock(thermo_mutex);
        if (!s_factory) {
            s_factory = new ThermoFactory;
        }
        return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
        ScopedLock lock(thermo_mutex);
        if (s_factory) {
            delete s_factory;
            s_factory = 0;
        }
    }

    //! Create a new thermodynamic property manager.
    /*!
     * @param model  String to look up the model against
     *
     * @return
     *   Returns a pointer to a new ThermoPhase instance matching the
     *   model string. Returns NULL if something went wrong.
     *   Throws an exception UnknownThermoPhaseModel if the string
     *   wasn't matched.
     */
    virtual ThermoPhase* newThermoPhase(const std::string& model);

private:
    //! static member of a single instance
    static ThermoFactory* s_factory;

    //! Private constructors prevents usage
    ThermoFactory() {};

    //! Decl for locking mutex for thermo factory singleton
    static mutex_t thermo_mutex;
};

//!  Create a new thermo manager instance.
/*!
 * @param model   String to look up the model against
 * @param f       ThermoFactor instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new ThermoPhase instance matching the
 *   model string. Returns NULL if something went wrong.
 *   Throws an exception UnknownThermoPhaseModel if the string
 *   wasn't matched.
 */
inline ThermoPhase* newThermoPhase(const std::string& model,
                                   ThermoFactory* f=0)
{
    if (f == 0) {
        f = ThermoFactory::factory();
    }
    return f->newThermoPhase(model);
}

//! Translate the eosType id into a string
/*!
 *  Returns a string representation of the eosType id for a phase.
 *  @param ieos  eosType id of the phase. This is unique for the phase
 *  @param length maximum length of the return string. Defaults to 100
 *
 *  @return returns a string representation.
 */
std::string eosTypeString(int ieos, int length = 100);

//! Create a new ThermoPhase object and initializes it according to the XML
//! tree.
/*!
 * This routine first looks up the identity of the model for the solution
 * thermodynamics in the model attribute of the thermo child of the xml phase
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
 */
ThermoPhase* newPhase(XML_Node& phase);

//! Create and Initialize a ThermoPhase object from an XML input file.
/*!
 *  This routine is a wrapper around the newPhase(XML_Node) routine
 *  which does the work. The wrapper locates the input phase XML_Node
 *  in a file, and then instantiates the object, returning the pointer
 *  to the ThermoPhase object.
 *
 * @param infile name of the input file
 * @param id     name of the phase id in the file.
 *               If this is blank, the first phase in the file is used.
 *
 * @return
 *   Returns an initialized ThermoPhase object.
 */
ThermoPhase* newPhase(const std::string& infile, std::string id);

//! Import a phase information into an empty thermophase object
/*!
 *   Here we read an XML description of the thermodynamic information
 *   for a phase. At the end of this routine, the phase should
 *   be ready to be used within applications. This routine contains
 *   some key routines that are used as pass back routines so that
 *   the phase (and the contents of the XML file) may contain
 *   variable parameterizations for the specification of the
 *   species standard states, the equation of state, and the
 *   specification of other nonidealities. Below, a description
 *   is presented of the main algorithm for bringing up a %ThermoPhase
 *   object, with care to present points where customizations
 *   occur.
 *
 *   Before invoking this routine, either the ThermoPhase Factory routines
 *   are called or direct constructor routines are called that
 *   instantiate an inherited ThermoPhase object. This object is input
 *   to this routine, and therefore contains inherited routines that
 *   drive the customization of the initialization process.
 *
 *   At the start of the routine, we import descriptions of the elements
 *   that make up the species in a phase.
 *
 *   We call setParametersFromXML(eos) to read parameters about
 *   the thermo phase before the species are read in.
 *
 *   We call addElementsFromXML() to add elements into the
 *   description of the phase.
 *
 *   We create a new species thermo manager.  Function
 *   'newSpeciesThermoMgr' looks at the species in the database
 *   to see what thermodynamic property parameterizations are
 *   used, and selects a class that can handle the
 *   parameterizations found.
 *
 *   We import information about the species, including their
 *   reference state thermodynamic polynomials. We then freeze
 *   the state of the species in the element.
 *
 *   Finally, we call initThermoXML(),
 *   a member function of the ThermoPhase object, to "finish"
 *   the description. Now that the species are known,
 *   additional information may be read in about the thermodynamics
 *   of the phase, (e.g.,  virial coefficients, which are
 *   binary or ternary interaction parameters between species).
 *
 * @param phase This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * @param th   Pointer to the ThermoPhase object which will
 *             handle the thermodynamics for this phase.
 *             We initialize part of the Thermophase object
 *             here, especially for those objects which are
 *             part of the Cantera Kernel.
 *
 * @param spfactory species Thermo factory pointer, if
 *                  available. If not available, one will be
 *                  created.
 *
 * @ingroup thermoprops
 */
bool importPhase(XML_Node& phase, ThermoPhase* th, SpeciesThermoFactory* spfactory = 0);

//! Install a species into a ThermoPhase object, which defines
//! the phase thermodynamics and speciation.
/*!
 *  This routine first gathers the information from the Species XML
 *  tree and calls addUniqueSpecies() to add it to the
 *  ThermoPhase object, p.
 *  This information consists of:
 *         ecomp[] = element composition of species.
 *         chgr    = electric charge of species
 *         name    = string name of species
 *         sz      = size of the species
 *                 (option double used a lot in thermo)
 *
 *  Then, the routine processes the "thermo" XML element and
 *  calls underlying utility routines to read the XML elements
 *  containing the thermodynamic information for the reference
 *  state of the species. Failures or lack of information trigger
 *  an "UnknownSpeciesThermoModel" exception being thrown.
 *
 * @param k     Species Index in the phase
 * @param s     XML_Node containing the species data for this species.
 * @param p     Reference to the ThermoPhase object.
 * @param spthermo_ptr Reference to the SpeciesThermo object, where
 *              the standard state thermo properties for this
 *              species will be installed.
 * @param rule  Parameter that handles what to do with species
 *              who have elements that aren't declared.
 *              Check that all elements in the species
 *              exist in 'p'. If rule != 0, quietly skip
 *              this species if some elements are undeclared;
 *              otherwise, throw an exception
 * @param phaseNode_ptr Pointer to the XML_Node for this phase
 *              (defaults to 0)
 * @param vpss_ptr pointer to the Manager that calculates standard
 *                state thermo properties
 * @param factory Pointer to the SpeciesThermoFactory .
 *              (defaults to 0)
 *
 * @return
 *  Returns true if everything is ok, false otherwise.
 */
bool installSpecies(size_t k, const XML_Node& s, thermo_t& p,
                    SpeciesThermo* spthermo_ptr, int rule,
                    XML_Node* phaseNode_ptr = 0,
                    VPSSMgr* vpss_ptr = 0,
                    SpeciesThermoFactory* factory = 0);

//!  Search an XML tree for species data.
/*!
 *   This utility routine will search the XML tree for the species
 *   named by the string, kname. It will return the XML_Node
 *   pointer to the species data for that species.
 *   Failures of any kind return the null pointer.
 *
 * @param kname String containing the name of the species.
 * @param phaseSpeciesData   Pointer to the XML speciesData element
 *              containing the species data for that phase.
 *
 *
 */
const XML_Node* speciesXML_Node(const std::string& kname,
                                const XML_Node* phaseSpeciesData);

//@}

}

#endif


