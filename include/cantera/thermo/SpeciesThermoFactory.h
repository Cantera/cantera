/**
 *  @file SpeciesThermoFactory.h
 *    Header for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species
 *    (see \ref spthermo and class
 *     \link Cantera::SpeciesThermoFactory SpeciesThermoFactory\endlink);
 */
// Copyright 2001  California Institute of Technology

#ifndef SPECIESTHERMO_FACTORY_H
#define SPECIESTHERMO_FACTORY_H

#include "SpeciesThermo.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/base/ct_thread.h"

namespace Cantera
{

class XML_Node;
class VPStandardStateTP;
class VPSSMgr;
class ThermoPhase;

/**
 * Throw a named error for an unknown or missing species thermo model.
 *
 * @ingroup thermoprops
 */
class UnknownSpeciesThermoModel: public CanteraError
{
public:
    //! constructor
    /*!
     * @param proc Function name error occurred.
     * @param spName Species Name that caused the error
     * @param speciesThermoModel Unrecognized species thermo model name
     */
    UnknownSpeciesThermoModel(const std::string& proc, const std::string& spName,
                              const std::string& speciesThermoModel) :
        CanteraError(proc, "species " + spName +
                     ": Specified speciesThermoPhase model "
                     + speciesThermoModel +
                     " does not match any known type.") {}
};

//! Factory to build instances of classes that manage the
//! standard-state thermodynamic properties of a set of species.
/*!
 *  This class is responsible for making the decision concerning
 *  which derivative of SpeciesThermo object to use.
 *  The SpeciesThermo object is used to calculate
 *  thermodynamic functions for the reference state.
 *  It queries the database of species to understand what
 *  the requirements are for the submodels for all of the
 *  species in the phase. Then, it picks the SpeciesThermo
 *  object to use, and passes it back to the calling routine.
 *  It doesn't load any of the data into the derived
 *  SpeciesThermo object.
 *
 *  Making the choice of SpeciesThermo types is the only
 *  thing this class does.
 *
 * This class is implemented as a singleton -- one in which
 * only one instance is needed.  The recommended way to access
 * the factory is to call this static method, which
 * instantiates the class if it is the first call, but
 * otherwise simply returns the pointer to the existing
 * instance.
 *
 * @ingroup thermoprops
 */
class SpeciesThermoFactory : public FactoryBase
{

public:

    //! Static method to return an instance of this class
    /*!
     * This class is implemented as a singleton -- one in which only one
     * instance is needed.  The recommended way to access the factory is to
     * call this static method, which instantiates the class if it is the
     * first call, but otherwise simply returns the pointer to the existing
     * instance.
     */
    static SpeciesThermoFactory* factory();

    //! Delete static instance of this class
    /**
     * If it is necessary to explicitly delete the factory before
     * the process terminates (for example, when checking for
     * memory leaks) then this method can be called to delete it.
     */
    void deleteFactory();

    //! Create a new species property manager for the reference state.
    /*!
     *  @param type the integer type to be created.
     *  @return  Returns the pointer to the newly allocated species property
     *           manager for the reference state
     */
    SpeciesThermo* newSpeciesThermo(int type) const;

    //! Create a new species thermo property manager given a string
    /*!
     *  @param stype  String name for the species thermo type
     *  @return       Returns the pointer to the newly malloced species
     *                property manager for the reference state
     */
    SpeciesThermo* newSpeciesThermoManager(std::string& stype) const;

    //! Create a new species property manager for the reference
    //! state for a group of species
    /*!
     * This routine will look through species nodes. It will discover what
     * each species needs for its species property managers. Then,
     * it will malloc and return the proper species property manager to use.
     *
     * @param spDataNodeList  This vector contains a list of species XML
     *                        nodes that will be in the phase
     * @return  Returns the pointer to the newly malloced species property
     *          manager for the reference state
     */
    SpeciesThermo* newSpeciesThermo(std::vector<XML_Node*> & spDataNodeList) const;

    //! Install a species thermodynamic property parameterization
    //! for the reference state for one species into a species thermo manager.
    /*!
     * @param k             Species number
     * @param speciesNode   Reference to the XML node specifying the species
     *                      standard state information
     * @param th_ptr        Pointer to the %ThermoPhase object for the species
     * @param spthermo      Species reference state thermo manager
     * @param phaseNode_ptr Optional pointer to the XML phase information for
     *                      the phase in which the species resides
     */
    void installThermoForSpecies(size_t k, const XML_Node& speciesNode,
                                 ThermoPhase* th_ptr, SpeciesThermo& spthermo,
                                 const XML_Node* phaseNode_ptr = 0) const;

    //! Install a species thermodynamic property parameterization
    //! for the standard state for one species into a species thermo manager, VPSSMgr
    /*!
     * This is a wrapper around the createInstallVPSS() function in the
     * VPStandardStateTP object.
     *
     * This serves to install the species into vpss_ptr, create a PDSS file. We also
     * read the xml database to extract the constants for these steps.
     *
     * @param k             species number
     * @param speciesNode   Reference to the XML node specifying the species
     *                      standard state information
     * @param vp_ptr        variable pressure ThermoPhase object
     * @param vpss_ptr      Pointer to the Manager for calculating variable
     *                      pressure substances.
     * @param spthermo_ptr  Species reference state thermo manager
     * @param phaseNode_ptr Optional Pointer to the XML phase information for
     *                      the phase in which the species resides
     */
    void installVPThermoForSpecies(size_t k, const XML_Node& speciesNode,
                                   VPStandardStateTP* vp_ptr,
                                   VPSSMgr* vpss_ptr,
                                   SpeciesThermo* spthermo_ptr,
                                   const XML_Node* phaseNode_ptr) const;

private:

    //! Pointer to the sole instance of this class, which is static
    static SpeciesThermoFactory* s_factory;

    //! Decl of the static mutex variable that locks the %SpeciesThermo factory singleton
    static mutex_t species_thermo_mutex;

    //! Constructor. This is made private, so that only the static
    //! method factory() can instantiate the class.
    SpeciesThermoFactory() {}
};


////////////////////// Convenience functions ////////////////////
//
//  These functions allow using a different factory class that
//  derives from SpeciesThermoFactory.
//
//////////////////////////////////////////////////////////////////


//! Create a new species thermo manager instance, by specifying the type and
//! (optionally) a pointer to the factory to use to create it.
/*!
 * This utility program  will look through species nodes. It will discover what
 * each species needs for its species property managers. Then,
 * it will malloc and return the proper species property manager to use.
 *
 * These functions allow using a different factory class that
 * derives from SpeciesThermoFactory.
 *
 * @param type         Species thermo type.
 * @param f            Pointer to a SpeciesThermoFactory. optional parameter.
 *                     Defaults to NULL.
 */
SpeciesThermo* newSpeciesThermoMgr(int type, SpeciesThermoFactory* f=0);

//! Create a new species thermo manager instance, by specifying the type and
//! (optionally) a pointer to the factory to use to create it.
/*!
 * This utility program is a basic factory operation for spawning a
 * new species reference-state thermo manager
 *
 * These functions allows for using a different factory class that
 * derives from SpeciesThermoFactory. However, no applications of this
 * have been done yet.
 *
 * @param stype       String specifying the species thermo type
 * @param f           Pointer to a SpeciesThermoFactory. optional parameter.
 *                    Defaults to NULL.
 */
SpeciesThermo* newSpeciesThermoMgr(std::string& stype,
                                   SpeciesThermoFactory* f=0);

//! Function to return SpeciesThermo manager
/*!
 * This utility program  will look through species nodes. It will discover what
 * each species needs for its species property managers. Then,
 * it will malloc and return the proper species reference state manager to use.
 *
 * These functions allow using a different factory class that
 * derives from SpeciesThermoFactory.
 *
 * @param spDataNodeList This vector contains a list of species XML nodes that
 *                       will be in the phase
 * @param f              Pointer to a SpeciesThermoFactory. optional
 *                       parameter. Defaults to NULL.
 */
SpeciesThermo* newSpeciesThermoMgr(std::vector<XML_Node*> spDataNodeList,
                                   SpeciesThermoFactory* f=0);

}

#endif


