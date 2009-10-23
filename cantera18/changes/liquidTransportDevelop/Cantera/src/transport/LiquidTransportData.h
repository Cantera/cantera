/**
 *  @file TransportFactory.h
 *  Header file defining class TransportFactory
 *     (see \link Cantera::TransportFactory TransportFactory\endlink)
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2008/12/24 18:19:01 $
 *  $Revision: 1.14 $
 *
 *
 *
 */

#ifndef CT_LIQUIDTRANSPORTDATA_H
#define CT_LIQUIDTRANSPORTDATA_H


// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <new>



// Cantera includes
#include "ct_defs.h"
#include "TransportBase.h"
#include "FactoryBase.h"


namespace Cantera {

  enum LiquidTR_Model {
    //! Temperature dependence type for pure (liquid) species properties
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    LTR_MODEL_NOTSET=-1,
    LTR_MODEL_CONSTANT, 
    LTR_MODEL_ARRHENIUS,
    LTR_MODEL_POLY
  };

  //! Class LiquidTransportData holds transport parameters for a 
  //! specific liquid-phase species.
  class LiquidTransportData {

  public:

    LiquidTransportData() : 
      speciesName("-"), 
      model_hydroradius(LTR_MODEL_NOTSET),
      model_viscosity(LTR_MODEL_NOTSET),
      model_thermalCond(LTR_MODEL_NOTSET),
      model_speciesDiffusivity(LTR_MODEL_NOTSET)
    {
    }

    std::string speciesName;
   
    //! Model type for the hydroradius
    LiquidTR_Model model_hydroradius;

    //! Ceofficients for the hydroradius model
    vector_fp  hydroRadiusCoeffs;

    //! Model type for the viscosity
    LiquidTR_Model model_viscosity;

    //! Ceofficients for the viscosity model
    vector_fp   viscCoeffs; 

    //! Model type for the thermal conductivity
    LiquidTR_Model model_thermalCond;
   
    //! Ceofficients for the thermal conductivity model
    vector_fp   thermalCondCoeffs;

    //! Model type for the speciesDiffusivity
    LiquidTR_Model model_speciesDiffusivity;
   
    //! Ceofficients for the species diffusivity model
    vector_fp   speciesDiffusivityCoeffs;
  };

}
#endif
