/**
 *  @file TransportBase.h  
 *   Provides class Transport.
 */

/*
 * $Revision: 1.23 $
 * $Date: 2009/03/27 18:24:39 $
 */

// Copyright 2001-2003  California Institute of Technology


/**
 * @defgroup tranprops Transport Properties
 *
 * @ingroup phases
 *
 * These classes provide transport properties.
 */

#ifndef CT_TRANSPORTBASE_H
#define CT_TRANSPORTBASE_H

#include "ThermoPhase.h"


namespace Cantera {

  class TransportParams;

  const int CK_Mode = 10;

  // types of transport models that can be constructed
  const int None                 = 199;
  const int cMulticomponent      = 200;
  const int CK_Multicomponent    = 202;
  const int cMixtureAveraged     = 210;
  const int CK_MixtureAveraged   = 211;
  const int cSolidTransport      = 300;
  const int cDustyGasTransport   = 400;
  const int cUserTransport       = 500;
  const int cFtnTransport        = 600;
  const int cLiquidTransport     = 700;
  const int cAqueousTransport    = 750;
  const int cRadiativeTransport  = 800;

  // forward reference
  class XML_Writer;


  /**
   * Base class for transport property managers.  All classes that
   * compute transport properties derive from this class.  Class
   * Transport is meant to be used as a base class only. It is
   * possible to instantiate it, but its methods throw exceptions if
   * called.
   */
  class Transport {

  public:


    /**
     * Constructor. New transport managers should be created using
     * TransportFactory, not by calling the constructor directly.
     * @see TransportFactory
     */
    Transport(thermo_t* thermo=0, int ndim = 1);

    //! Destructor.
    virtual ~Transport();   

    //!  Copy Constructor for the %Transport  object.
    /*!
     * @param right  Transport to be copied
     */
    Transport(const Transport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to Transport object to be copied into the
     *                 current one.
     */
    Transport&  operator=(const Transport& right);

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
    // Note ->need working copy constructors and operator=() functions for all first
    virtual Transport *duplMyselfAsTransport() const;


    /**
     * Transport model. The transport model is the set of
     * equations used to compute the transport properties. This
     * virtual method returns an integer flag that identifies the
     * transport model implemented. The base class returns 0.
     */
    virtual int model() {return 0;}

    /**
     * Phase object. Every transport manager is designed to compute
     * properties for a specific phase of a mixture, which might be a
     * liquid solution, a gas mixture, a surface, etc. This method
     * returns a reference to the object representing the phase
     * itself.
     */        
    thermo_t& thermo() { return *m_thermo; }


    /**
     * Returns true if the transport manager is ready for use.
     */
    bool ready();

    /**
     * Returns an integer index number. This is for internal use
     * of Cantera, and may be removed in the future.
     */
    int index() const ;

    /**
     * Set an integer index number. This is for internal use of
     * Cantera, and may be removed in the future.
     */
    void setIndex(int i);

    //! Set the number of dimensions to be expected in flux expressions
    /*!
     * Internal memory will be set with this value
     */
    void setNDim(const int ndim);

    //! return the number of dimensions
    int nDim() const { return m_nDim; }

    /**
     * @name Transport Properties
     */
    //@{
         

    /**
     * The viscosity in Pa-s. 
     */
    virtual doublereal viscosity() 
    { return err("viscosity"); }


    /**
     * The bulk viscosity in Pa-s. The bulk viscosity is only
     * non-zero in rare cases. Most transport managers either
     * overload this method to return zero, or do not implement
     * it, in which case an exception is thrown if called.
     */
    virtual doublereal bulkViscosity()  
    { return err("bulkViscosity"); }

        
    /**
     * The thermal conductivity in W/m/K. 
     */
    virtual doublereal thermalConductivity()
    { return err("thermalConductivity"); }

    /**
     * The electrical conductivity (Siemens/m).
     */
    virtual doublereal electricalConductivity()
    { return err("electricalConductivity"); }

    /**
     * Electrical mobilities (m^2/V/s). Returns the mobilities of
     * the species in array \c mobil. The array must be
     * dimensioned at least as large as the number of species.
     */
    virtual void getMobilities(doublereal* const mobil)
    { err("getMobilities"); }


    //@}


    //! Get the species diffusive mass fluxes wrt to 
    //! the mass averaged velocity, 
    //! given the gradients in mole fraction and temperature
    /*!
     *  Units for the returned fluxes are kg m-2 s-1.
     * 
     *  @param ndim Number of dimensions in the flux expressions
     *  @param grad_T Gradient of the temperature
     *                 (length = ndim)
     * @param ldx  Leading dimension of the grad_X array 
     *              (usually equal to m_nsp but not always)
     * @param grad_X Gradients of the mole fraction
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array 
     *              (usually equal to m_nsp but not always)
     * @param fluxes  Output of the diffusive mass fluxes
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesFluxes(int ndim, 
				  const doublereal* grad_T, 
				  int ldx, 
				  const doublereal* grad_X,
				  int ldf, 
				  doublereal* fluxes) { 
      err("getSpeciesFluxes"); 
    }

    /** 
     * Get the molar fluxes [kmol/m^2/s], given the thermodynamic
     * state at two nearby points. 
     * @param state1 Array of temperature, density, and mass
     * fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     * fractions for state 2.  
     * @param delta Distance from state 1 to state 2 (m).
     */ 
    virtual void getMolarFluxes(const doublereal* state1,
				const doublereal* state2, doublereal delta, 
				doublereal* fluxes) { err("getMolarFluxes"); }

    /** 
     * Get the mass fluxes [kg/m^2/s], given the thermodynamic
     * state at two nearby points. 
     * @param state1 Array of temperature, density, and mass
     * fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     * fractions for state 2.  
     * @param delta Distance from state 1 to state 2 (m).
     */ 
    virtual void getMassFluxes(const doublereal* state1,
			       const doublereal* state2, doublereal delta, 
			       doublereal* fluxes) { err("getMassFluxes"); }

    /**
     * Thermal diffusion coefficients [kg/m/sec].
     * The thermal diffusion coefficient \f$ D^T_k \f$ is defined
     * so that the diffusive mass flux of species k induced by the
     * local temperature gradient is \f[ M_k J_k = -D^T_k \nabla
     * \ln T. \f]. The thermal diffusion coefficient can be either
     * positive or negative.
     * 
     * @param dt on return, dt will contain the species thermal
     * diffusion coefficients.  Dimension dt at least as large as
     * the number of species.
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt) 
    { err("getThermalDiffCoeffs"); }


    /**
     * Binary diffusion coefficients [m^2/s].
     */
    virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d) 
    { err("getBinaryDiffCoeffs"); }


    /**
     * Multicomponent diffusion coefficients. Units: [m^2/s].  If
     * the transport manager implements a multicomponent diffusion
     * model, then this method returns the array of multicomponent
     * diffusion coefficients. Otherwise it throws an exception.
     */
    virtual void getMultiDiffCoeffs(const int ld, doublereal* const d) 
    { err("getMultiDiffCoeffs"); }


    /**
     * Mixture-averaged diffusion coefficients [m^2/s].  If the
     * transport manager implements a mixture-averaged diffusion
     * model, then this method returns the array of
     * mixture-averaged diffusion coefficients. Otherwise it
     * throws an exception.
     */
    virtual void getMixDiffCoeffs(doublereal* const d) 
    { err("getMixDiffCoeffs"); }


    /**
     * Set transport model parameters. This method may be
     * overloaded in subclasses to set model-specific parameters.
     */
    virtual void setParameters(const int type, const int k,
			       const doublereal* const p); 
   

    friend class TransportFactory;


  protected:

    /**
     * @name Transport manager construction
     * These methods are used internally during construction.  
     * @{
     */

    /**
     * Called by TransportFactory to set parameters.
     */
    virtual bool init(TransportParams& tr)
    { err("init"); return false; }


    /**
     * Set the phase object. 
     */
    void setThermo(thermo_t& thermo);


    /** 
     * Enable for use. Once finalize() has been called, the
     * transport manager should be ready to compute any supported
     * transport property, and no further modifications to the
     * model parameters should be made.
     */
    void finalize();

    //@}


    thermo_t*  m_thermo;  ///< pointer to the object representing the phase 
    bool      m_ready;    ///< true if finalize has been called
    size_t    m_nmin;     ///< number of species
    int       m_index;

    //! Number of dimensions used in flux expresions
    int       m_nDim;


  private:

    /**
     * Throw an exception if a method of this class is
     * invoked. This probably indicates that a transport manager
     * is being used that does not implement all virtual methods,
     * and one of those methods was called by the application
     * program. For example, a transport manager that computes the
     * thermal conductivity of a solid may not define the
     * viscosity() method, since the viscosity is in this case
     * meaningless. If the application invokes the viscosity()
     * method, the base class method will be called, resulting in
     * an exception being thrown.
     */
    doublereal err(std::string msg) const;

  };

  typedef Transport transport_t;

}

#endif
