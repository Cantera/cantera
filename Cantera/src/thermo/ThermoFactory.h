/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef THERMO_FACTORY_H
#define THERMO_FACTORY_H

#include "ThermoPhase.h"
#include "xml.h"

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

//#include "SpeciesThermoFactory.h"
#include "FactoryBase.h"

namespace Cantera {

    class SpeciesThermoFactory;

  /*!
   *  @addtogroup thermoprops
   *
   *  Standard %ThermoPhase objects may be instantiated by calling
   *  the main %Cantera factory class for %ThermoPhase objects; This class is called ThermoFactory.
   */
  //@{

  //! Specific error to be thrown if the type of Thermo mananger is unrecognized.
  /*!
   * This particular error class may be caught, if the application may have other
   * models that the main Cantera appliation doesn't know about.
   */
  class UnknownThermoPhaseModel : public CanteraError {
  public:
    //! Constructor
    /*!
     * @param proc Function name where the error occurred.
     * @param thermoModel Sting name of ThermoPhase which didn't match
     */
    UnknownThermoPhaseModel(std::string proc, std::string thermoModel) :
      CanteraError(proc, "Specified ThermoPhase model "   
		   + thermoModel + 
		   " does not match any known type.") {}
    //! destructor
    virtual ~UnknownThermoPhaseModel() {}
  };

  
  //! Factory class for thermodynamic property managers.
  /*!
   * This class keeps a list of the known ThermoPhase classes, and is
   * used to create new instances of these classes.
   */
    class ThermoFactory : public FactoryBase {

  public:

    //! Static function that creates a static instance of the factory.
    static ThermoFactory* factory() {
#if defined(THREAD_SAFE_CANTERA)
        boost::mutex::scoped_lock lock(thermo_mutex);
#endif
      if (!s_factory) s_factory = new ThermoFactory;
      return s_factory;
    }

      //! delete the static instance of this factory
      virtual void deleteFactory() {
#if defined(THREAD_SAFE_CANTERA)
          boost::mutex::scoped_lock lock(thermo_mutex);
#endif
          if (s_factory) {
            delete s_factory;
            s_factory = 0;
        }
    }
    
    //! Destructor doesn't do anything.
    /*!
     * We do not delete statically created single instance of this
     * class here, because it would create an infinite loop if
     * destructor is called for that single instance.
     */
    virtual ~ThermoFactory() { }

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
    virtual ThermoPhase* newThermoPhase(std::string model);

  private:
    //! static member of a single instance
    static ThermoFactory* s_factory;

    //! Private constructor prevents usage
      ThermoFactory(){}

#if defined(THREAD_SAFE_CANTERA)
        static boost::mutex thermo_mutex;
#endif

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
    inline ThermoPhase* newThermoPhase(std::string model,  
				     ThermoFactory* f=0) {
        if (f == 0) {
            f = ThermoFactory::factory();
        }
        return f->newThermoPhase(model);
    }


  /*!
   *  This routine first looks up the
   * identity of the model for the solution thermodynamics in the
   * model attribute of the thermo child of the xml phase
   * node. Then, it does a string lookup using Cantera's internal ThermoPhase Factory routines
   * on the model to figure out
   * what ThermoPhase derived class should be assigned. It creates a new
   * instance of that class, and then calls importPhase() to
   * populate that class with the correct parameters from the XML
   * tree.
   *
   * @param phase XML_Node reference pointing to the phase XML element.
   *
   * @return
   *    Returns a pointer to the completed and initialized ThermoPhase object. 
   *
   * @ingroup inputfiles
   */
    ThermoPhase* newPhase(XML_Node& phase);
    ThermoPhase* newPhase(std::string infile, std::string id);

  //! Import a phase information into an empty thermophase object
  /*!
   *   Here we read an XML description of the thermodynamic information
   *   for a phase. At the end of this routine, the phase should
   *   be ready to be used within applications. This routine contains
   *   some key routines that are used as pass back routines so that
   *   the phase (and the contents of the XML file) may contain
   *   variable paramerizations for the specification of the
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
   *   drive the custimation of the initialization process.
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
  bool importPhase(XML_Node& phase, ThermoPhase* th, 
		   SpeciesThermoFactory* spfactory = 0);

    bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
        SpeciesThermo& spthermo, int rule, 
        SpeciesThermoFactory* factory = 0);

    const XML_Node *speciesXML_Node(std::string kname,
        const XML_Node *phaseSpeciesData);
  //@}

}

#endif


