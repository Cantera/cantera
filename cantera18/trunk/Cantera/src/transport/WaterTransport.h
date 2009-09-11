/**
 *
 *  @file LiquidTransport.h
 *   Header file defining class LiquidTransport
 */
/*
 * $Revision: 1.9 $
 * $Date: 2009/03/27 18:24:39 $
 */

#ifndef CT_WATERTRAN_H
#define CT_WATERTRAN_H



// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

using namespace std;

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"
#include "LiquidTransportParams.h"

namespace Cantera {

  const int LVISC_CONSTANT     = 0;
  const int LVISC_WILKES       = 1;
  const int LVISC_MIXTUREAVG   = 2;

  const int LDIFF_MIXDIFF_UNCORRECTED     = 0;
  const int LDIFF_MIXDIFF_FLUXCORRECTED  = 1;
  const int LDIFF_MULTICOMP_STEFANMAXWELL  = 2;



  class TransportParams;


  class WaterTransport : public Transport {
  public:

    //! default constructor
    WaterTransport(thermo_t* thermo = 0, int ndim = 1);

    //!Copy Constructor for the %LiquidThermo object.
    /*!
     * @param right  ThermoPhase to be copied
     */
    WaterTransport(const WaterTransport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %ThermoPhase object to be copied into the
     *                 current one.
     */
    WaterTransport&  operator=(const  WaterTransport& right);
    
    //! Duplication routine for objects which inherit from
    //! %Transport
    /*!
     *  This virtual routine can be used to duplicate %Transport objects
     *  inherited from %Transport even if the application only has
     *  a pointer to %Transport to work with.
     *
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     */
    virtual Transport *duplMyselfAsTransport() const;

    //! virtual destructor
    virtual ~WaterTransport();

    //! Return the model id for this transport parameterization
    virtual int model() {
      return cWaterTransport; 
    }

    //! overloaded base class methods

    //! Returns the viscosity of the solution
    /*!
     * The viscosity is computed using the Wilke mixture rule.
     * \f[
     * \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
     * \f]
     * Here \f$ \mu_k \f$ is the viscosity of pure species \e k,
     * and 
     * \f[
     * \Phi_{k,j} = \frac{\left[1 
     * + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
     * {\sqrt{8}\sqrt{1 + M_k/M_j}}
     * \f] 
     * @see updateViscosity_T();
     *
     * Controlling update boolean m_viscmix_ok
     */ 
    virtual doublereal viscosity();

   

  };
}
#endif






