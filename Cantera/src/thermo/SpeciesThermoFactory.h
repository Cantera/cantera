/**
 *  @file SpeciesThermoFactory.h
 *    Header for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species 
 *    (see \ref spthermo and class \link Cantera::SpeciesThermoFactory SpeciesThermoFactory\endlink);
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef SPECIESTHERMO_FACTORY_H
#define SPECIESTHERMO_FACTORY_H

#include "SpeciesThermo.h"
#include "ctexceptions.h"
#include "FactoryBase.h"


namespace Cantera {

  class XML_Node;
  class VPStandardStateTP;
  class VPSSMgr;

  /**
   * Throw a named error for an unknown or missing species thermo model. 
   *
   * @ingroup thermoprops
   */
  class UnknownSpeciesThermoModel: public CanteraError {
  public:
    //! constructor
    /*!
     * @param proc Function name error occurred.
     * @param spName Species Name that caused the error
     * @param speciesThermoModel Unrecognized species thermo model name
     */
    UnknownSpeciesThermoModel(std::string proc, std::string spName,
			      std::string speciesThermoModel) :
      CanteraError(proc, "species " + spName + 
		   ": Specified speciesThermoPhase model "   
		   + speciesThermoModel + 
		   " does not match any known type.") {}
    //! destructor
    virtual ~UnknownSpeciesThermoModel() {}
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
   *  object to use, and passies it back to the calling routine.
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
  class SpeciesThermoFactory : public FactoryBase {

  public:

    //! Static method to return an instance of this class
    /*!
     * This class is implemented as a singleton -- one in which
     * only one instance is needed.  The recommended way to access
     * the factory is to call this static method, which
     * instantiates the class if it is the first call, but
     * otherwise simply returns the pointer to the existing
     * instance.
     */
    static SpeciesThermoFactory* factory() {
#if defined(THREAD_SAFE_CANTERA)
        boost::mutex::scoped_lock lock(species_thermo_mutex);
#endif
      if (!s_factory) s_factory = new SpeciesThermoFactory;
      return s_factory;
    }

    //! Delete static instance of this class
    /**
     * If it is necessary to explicitly delete the factory before
     * the process terminates (for example, when checking for
     * memory leaks) then this method can be called to delete it.
     */
    void deleteFactory() {
#if defined(THREAD_SAFE_CANTERA)
      boost::mutex::scoped_lock lock(species_thermo_mutex);
#endif
      if (s_factory) {
	delete s_factory;
	s_factory = 0;
      }
    }
	
    //! Destructor
    /**
     * Doesn't do anything. We do not delete statically
     * created single instance of this class here, because it would
     * create an infinite loop if destructor is called for that
     * single instance.
     */
    virtual ~SpeciesThermoFactory() {
    }

    //! Create a new species property manager.
    /*!
     * @param type the integer type to be created.
     */ 
    SpeciesThermo* newSpeciesThermo(int type);

    //! Create a new species property manager.
    /*!
     * This routine will look through species nodes. It will discover what
     * each species needs for its species property managers. Then,
     * it will malloc and return the proper species property manager to use.
     *
     * @param spData_node  Pointer to a speciesData XML Node.
     *                     Each speciesData node contains a list of XML species elements
     *                      e.g., \<speciesData id="Species_Data"\>
     */ 
    SpeciesThermo* newSpeciesThermo(XML_Node* spData_node);

    //! Create a new species property manager for a group of species
    /*!
     * This routine will look through species nodes. It will discover what
     * each species needs for its species property managers. Then,
     * it will malloc and return the proper species property manager to use.
     *
     * @param spData_nodes Vector of XML_Nodes, each of which is a speciesData XML Node.
     *                     Each speciesData node contains a list of XML species elements
     *                      e.g., \<speciesData id="Species_Data"\>
     */ 
    SpeciesThermo* newSpeciesThermo(std::vector<XML_Node*> spData_nodes);

    //! Create a new species property manager.
    /*!
     * This routine will look through species nodes. It will discover what
     * each species needs for its species property managers. Then,
     * it will malloc and return the proper species property manager to use.
     *
     *
     * @param spData_nodes Vector of XML_Nodes, each of which is a speciesData XML Node.
     *                     Each %speciesData node contains a list of XML species elements
     *                      e.g., \<speciesData id="Species_Data"\>
     *
     *  @todo is this used? 
     */ 
    SpeciesThermo* newSpeciesThermoOpt(std::vector<XML_Node*> spData_nodes);

    //! Install a species thermodynamic property parameterization
    //! for the reference state for one species into a species thermo manager.
    /*!
     * @param k species number
     * @param speciesNode  Reference to the XML node specifying the species standard
     *           state information
     * @param spthermo Species reference state thermo manager
     * @param phaseNode_ptr Optional Pointer to the XML phase
     *                      information for the phase in which the species
     *                      resides
     */
    void installThermoForSpecies(int k, const XML_Node& speciesNode, 
				 SpeciesThermo& spthermo,
				 const XML_Node *phaseNode_ptr = 0);

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
     * @param speciesNode   Reference to the XML node specifying the species standard
     *                      state information
     * @param vp_ptr        variable pressure ThermoPhase object 
     * @param vpss_ptr      Pointer to the Manager for calculating variable pressure
     *                      substances.
     * @param spthermo_ptr  Species reference state thermo manager
     * @param phaseNode_ptr Optional Pointer to the XML phase
     *                      information for the phase in which the species
     *                      resides
     */
    void installVPThermoForSpecies(int k, const XML_Node& speciesNode, 
				   VPStandardStateTP *vp_ptr,
				   VPSSMgr *vpss_ptr,
				   SpeciesThermo *spthermo_ptr,
				   const XML_Node *phaseNode_ptr);

  private:

    //! pointer to the sole instance of this class
    static SpeciesThermoFactory* s_factory;

#if defined(THREAD_SAFE_CANTERA)
    //! Decl of the static mutex variable that locks the %SpeciesThermo factory singelton
    static boost::mutex species_thermo_mutex;
#endif

    //! Constructor. This is made private, so that only the static
    //! method factory() can instantiate the class.
    SpeciesThermoFactory(){}
  };


  ////////////////////// Convenience functions ////////////////////
  //
  //  These functions allow using a different factory class that
  //  derives from SpeciesThermoFactory.
  //
  //////////////////////////////////////////////////////////////////


  //! Create a new species thermo manager instance, by specifying
  //!the type and (optionally) a pointer to the factory to use to create it.
  /*!
   * This utility program  will look through species nodes. It will discover what
   * each species needs for its species property managers. Then,
   * it will malloc and return the proper species property manager to use.
   *
   *  These functions allow using a different factory class that
   *  derives from SpeciesThermoFactory.
   *
   * @param type         Species thermo type.
   * @param f            Pointer to a SpeciesThermoFactory. optional parameter. 
   *                    Defautls to NULL.
   */
  inline SpeciesThermo* newSpeciesThermoMgr(int type, 
					    SpeciesThermoFactory* f=0) {
    if (f == 0) {
      f = SpeciesThermoFactory::factory();
    }
    SpeciesThermo* sptherm = f->newSpeciesThermo(type);
    return sptherm;
  }

  //! Function to return SpeciesThermo manager
  /*!
   * This utility program  will look through species nodes. It will discover what
   * each species needs for its species property managers. Then,
   * it will malloc and return the proper species property manager to use.
   *
   *  These functions allow using a different factory class that
   *  derives from SpeciesThermoFactory.
   *
   * @param spData_node Vector of XML_Nodes, each of which is a speciesData XML Node.
   *                     Each %speciesData node contains a list of XML species elements
   *                      e.g., \<speciesData id="Species_Data"\>
   * @param f            Pointer to a SpeciesThermoFactory. optional parameter. 
   *                    Defautls to NULL.
   */
  inline SpeciesThermo* newSpeciesThermoMgr(XML_Node* spData_node, 
					    SpeciesThermoFactory* f=0) {
    if (f == 0) {
      f = SpeciesThermoFactory::factory();
    }
    SpeciesThermo* sptherm = f->newSpeciesThermo(spData_node);
    return sptherm;
  }

  //! Function to return SpeciesThermo manager
  /*!
   * This utility program  will look through species nodes. It will discover what
   * each species needs for its species property managers. Then,
   * it will malloc and return the proper species property manager to use.
   *
   *  These functions allow using a different factory class that
   *  derives from SpeciesThermoFactory.
   *
   * @param spData_nodes Vector of XML_Nodes, each of which is a speciesData XML Node.
   *                     Each %speciesData node contains a list of XML species elements
   *                      e.g., \<speciesData id="Species_Data"\>
   * @param f            Pointer to a SpeciesThermoFactory. optional parameter. 
   *                    Defautls to NULL.
   * @param opt         Boolean defaults to false.
   */
  inline SpeciesThermo* newSpeciesThermoMgr(std::vector<XML_Node*> spData_nodes, 
					    SpeciesThermoFactory* f=0, bool opt=false) {
    if (f == 0) {
      f = SpeciesThermoFactory::factory();
    }
    SpeciesThermo* sptherm;
    if (opt) {
      sptherm = f->newSpeciesThermoOpt(spData_nodes);
    } else { 
      sptherm = f->newSpeciesThermo(spData_nodes);
    }
    return sptherm;
  }

}

#endif


