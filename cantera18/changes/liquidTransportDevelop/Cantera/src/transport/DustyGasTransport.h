/**
 *  @file DustyGasTransport.h
 *    Headers for the DustyGasTransport object, which models transport properties
 *    in porous media using the dusty gas approximation
 *    (see \ref tranprops and \link Cantera::DustyGasTransport DustyGasTransport \endlink) .
 *
 */
/*
 * $Revision$
 * $Date$
 */

// Copyright 2003  California Institute of Technology


#ifndef CT_DUSTYGASTRAN_H
#define CT_DUSTYGASTRAN_H

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"


namespace Cantera {

  //! Class DustyGasTransport implements the Dusty Gas model for transport in porous media.
  /*!
   *    As implemented here, only species transport is handled. The viscosity, thermal conductivity, and thermal 
   *    diffusion coefficients are not implemented.
   */
  class DustyGasTransport : public Transport {

  public:

    //! default constructor
    /*!
     *  @param thermo   Pointer to the %ThermoPhase object for this phase. Defaults to zero.
     */
    DustyGasTransport(thermo_t* thermo=0);

    //!   Copy Constructor for the %DustyGasTransport object.
    /*!
     *    @param right  %LiquidTransport to be copied
     */
    DustyGasTransport(const DustyGasTransport &right);

    //! Assignment operator
    /*!
     *
     *    Warning -> Shallow pointer copies are made of m_thermo and m_gastran.. gastran may not point to the correct 
     *               object after this copy. The routine initialize() must be called after this
     *               routine to complete the copy.
     *
     * @param right    Reference to %DustyGasTransport object to be copied 
     *                 into the current one.
     */
    DustyGasTransport& operator=(const  DustyGasTransport& right);
        
    //! Destructor.
    virtual ~DustyGasTransport();

    //! Duplication routine for objects which inherit from %Transport
    /*!
     *  This virtual routine can be used to duplicate %Transport objects
     *  inherited from %Transport even if the application only has
     *  a pointer to %Transport to work with.
     *
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     */
    virtual Transport *duplMyselfAsTransport() const;
        
    //---------------------------------------------------------
    // overloaded base class methods

    virtual int model() const { return cDustyGasTransport; }

    virtual void setParameters(const int type, const int k, const doublereal* const p);
        
    virtual void getMultiDiffCoeffs(const int ld, doublereal* const d);
        
    //! Get the molar fluxes [kmol/m^2/s], given the thermodynamic
    //! state at two nearby points. 
    /*!
     * @param state1 Array of temperature, density, and mass
     *               fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     *               fractions for state 2.  
     * @param delta  Distance from state 1 to state 2 (m).
     */ 
    virtual void getMolarFluxes(const doublereal * const state1,
				const doublereal* const state2, const doublereal delta, 
				doublereal* const fluxes);
        
    //-----------------------------------------------------------
    // new methods added in this class
        
    /// Set the porosity (dimensionless)
    void setPorosity(doublereal porosity) {
      m_porosity = porosity;
      m_knudsen_ok = false;
      m_bulk_ok = false;
    }

    /// Set the tortuosity (dimensionless)
    void setTortuosity(doublereal tort) {
      m_tortuosity = tort;
      m_knudsen_ok = false;
      m_bulk_ok = false;
    }

    /// Set the mean pore radius (m)
    void setMeanPoreRadius(doublereal rbar) {
      m_pore_radius = rbar;
      m_knudsen_ok = false;
    }

    /// Set the mean particle diameter
    void setMeanParticleDiameter(doublereal dbar) {
      m_diam = dbar;
    }
        
    //! Set the permeability of the media
    /*!
     * If not set, the value for close-packed spheres will be used by default. 
     *
     *  The value for close-packed spheres is given below, where p is the porosity,
     *  t is the tortuosity, and d is the diameter of the sphere
     *
     *  \f[
     *      \kappa = \frac{p^3 d^2}{72 t (1 - p)^2}
     *  \f]
     *
     * @param B  set the permeability of the media (units = m^2)
     */
    void setPermeability(doublereal B) {
      m_perm = B;
    }
        
    //! Return a reference to the transport manager used to compute the gas
    //! binary diffusion coefficients and the visdcosity.
    Transport& gasTransport() { return *m_gastran; }
        
        
    friend class TransportFactory;
        
        
  protected:

    //!  Initialization routine called by TransportFactory
    /*!
     *  The DustyGas model is a subordinate model to the gas phase transport model. Here we 
     *  set the gas phase models.
     *
     *  This is a protected routine, so that initialiation of the Model must occur within Cantera's setup
     *
     *   @param  phase           Pointer to the underlying ThermoPhase model for the gas phase
     *   @param  gastr           Pointer to the underlying Transport model for transport in the gas phse.
     */
    void initialize(ThermoPhase* phase, Transport* gastr);
        
        
  private:

    //! Update temperature-dependent quantities within the object
    /*!
     *  The object keeps a value m_temp, which is the temperature at which quantities were last evaluated
     *  at. If the temperature is changed, update Booleans are set false, triggering recomputation.
     */
    void updateTransport_T();
    void updateTransport_C();

    void updateBinaryDiffCoeffs();
    void updateMultiDiffCoeffs();
    void updateKnudsenDiffCoeffs();
    void eval_H_matrix();


    // gas attributes
    int m_nsp;
    doublereal m_tmin, m_tmax;
    vector_fp  m_mw;

  
    //! binary diffusion coefficients
    DenseMatrix                  m_d;

    //! mole fractions
    vector_fp                    m_x;

    //! Knudsen diffusion coefficients
    /*!
     *  The Knudsen diffusion coefficients are given by the following form
     *  
     *     \f[
     *        \mathcal{D}^{knud}_k =  \frac{2}{3} \frac{r_{pore} \phi}{\tau} \left( \frac{8 R T}{\pi W_k}  \right)^{1/2}
     *     \f]
     *
     */
    vector_fp                    m_dk;

    //! temperature
    doublereal                   m_temp;

    //! Multicomponent diffusion coefficients
    /*!
     *  The multicomponent diffusion matrix \f$  H_{k,l} \f$ is given by the following form
     *
     *     \f[
     *        H_{k,l} = - \frac{X_k}{D_{k,l}}
     *     \f]
     *     \f[
     *        H_{k,k} = \frac{1}{\mathcal(D)^{knud}_{k}} + \sum_{j \ne k}^N{ \frac{X_j}{D_{k,j}} }
     *     \f]
     */
    DenseMatrix                  m_multidiff;

    //!  work space of size m_nsp;
    vector_fp  m_spwork;
    
    //!  work space of size m_nsp;
    vector_fp  m_spwork2;

 
    //! Pressure Gradient
    doublereal m_gradP;

    //! Update-to-date variable for Knudsen diffusion coefficients
    bool m_knudsen_ok;

   //! Update-to-date variable for Binary diffusion coefficients
    bool m_bulk_ok;

    bool m_gradConc_set;

    bool m_gradP_set;

    //! Porosity
    doublereal m_porosity;   

    //! Tortuosity
    doublereal m_tortuosity;
    doublereal m_pore_radius;   /// pore radius (m)

    //! Particle diameter
    /*!
     *   The medium is assumed to consist of particles of size m_diam
     *   units =  m
     */
    doublereal m_diam;

    //! Permeability of the media
    /*!
     *  The permeability is the proportionality constant for Darcy's
     *  law which relates discharge rate and viscosity to the applied
     *  pressure gradient.
     *
     *   Below is Darcy's law, where \f$ \kappa \f$ is the permeability
     *
     *     \f[
     *         v = \frac{\kappa}{\mu} \frac{\delta P}{\delta x}
     *     \f]
     *  
     *  units are m2
     */
    doublereal m_perm; 

    Transport* m_gastran;       /// pointer to gas transport manager

  };
}
#endif






