/**
 *  @file ThermoFactory.h
 *
 *         This file contains the definition for the factory
 *         class that can create know %ThermoPhase objects.
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


namespace Cantera {

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
   * This class keeps a list of the known ThermoPhase classes, and is used
   * to create new instances of these classes.
   */
  class ThermoFactory {

  public:

    //! Static function that creates a static instance of the factor.
    static ThermoFactory* factory() {
      if (!s_factory) s_factory = new ThermoFactory;
      return s_factory;
    }

    //! delete the static instance of this factory
    static void deleteFactory() {
      if (s_factory) {
	delete s_factory;
	s_factory = 0;
      }
    }
    
    //! Destructor doesn't do anything.
    /*!
     * We do not delete statically
     * created single instance of this class here, because it would
     * create an infinite loop if destructor is called for that
     * single instance.
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

  //@}

}

#endif


