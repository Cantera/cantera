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
 *  Copyright 2001 California Institute of Technology
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
    LTR_MODEL_NOTSET=-1,
    LTR_MODEL_CONSTANT, 
    LTR_MODEL_ARRHENIUS,
    LTR_MODEL_COEFF
  };

  class LiquidTransportData {

  public:

    LiquidTransportData() : 
      speciesName("-"), 
      model_hydroradius(LTR_MODEL_NOTSET),
      hydroradius(-1.0),
      model_viscosity(LTR_MODEL_NOTSET),
      model_thermalCond(LTR_MODEL_NOTSET),
      model_speciesDiffusivity(LTR_MODEL_NOTSET)
    {
    }

    std::string speciesName;
   
    //! Model type for the hydroradius
    LiquidTR_Model model_hydroradius;

    //! Actual value of the hydroradius
    doublereal  hydroradius;

    //! Model type for the hydroradius
    LiquidTR_Model model_viscosity;
    vector_fp   viscCoeffs;

    //! Model type for the hydroradius
    LiquidTR_Model model_thermalCond;
   
    vector_fp   thermalCondCoeffs;

    //! Model type for the hydroradius
    LiquidTR_Model model_speciesDiffusivity;
   
    vector_fp   speciesDiffusivityCoeffs;
  };

}
#endif
