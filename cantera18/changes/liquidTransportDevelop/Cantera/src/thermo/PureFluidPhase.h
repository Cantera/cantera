/**
 *  @file PureFluidPhase.h
 *
 *   Header for a ThermoPhase class for a pure fluid phase consisting of 
 *   gas, liquid, mixed-gas-liquid and supercrit fluid (see \ref thermoprops 
 *   and class \link Cantera::PureFluidPhase PureFluidPhase\endlink).
 *
 * This class is only available if the WITH_PURE_FLUIDS optional compile
 * capability has been turned on in Cantera's makefile system. 
 * It inherits from ThermoPhase, but is built on top of the tpx package.
 */

/*  $Author: dggoodwin $
 *  $Date: 2007/12/24 15:32:30 $
 *  $Revision: 1.7 $
 *
 *  Copyright 2003 California Institute of Technology
 */

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"
/*
 * For doxygen to recognize the ifdef below, it seems necessary to 
 * include a specific include reference.
 */
/**
 * This object is only available if the WITH_PURE_FLUIDS optional compile
 * capability has been turned on in Cantera's makefile system.  
 */
#ifdef WITH_PURE_FLUIDS

#include "mix_defs.h"

namespace tpx {
    class Substance;
}

namespace Cantera {

  //!   This phase object consists of a single component that can be a
  //!   gas, a liquid, a mixed gas-liquid fluid, or a fluid beyond its
  //!   critical point
  /*!
   *  The object inherits from ThermoPhase. However, its build on top
   *  of the tpx package.
   *
   *
   * <H2> Specification of Species Standard State Properties </H2>
   *
   *
   * <H2> Application within %Kinetics Managers </H2>
   *
   *
   * <H2> XML Example </H2>
   *
   *
   * <H2> Instantiation of the Class </H2>
   *
   * @ingroup thermoprops
   */
  class PureFluidPhase  : public ThermoPhase {

  public:

    //! Empty Base Constructor
    PureFluidPhase();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    PureFluidPhase(const PureFluidPhase &right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    PureFluidPhase& operator=(const PureFluidPhase& right);

    //! Destructor
    virtual ~PureFluidPhase();

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase *duplMyselfAsThermoPhase() const;

    //! Equation of state type
    virtual int eosType() const { return cPureFluid; }

    /// Molar enthalpy. Units: J/kmol. 
    virtual doublereal enthalpy_mole() const;

    /// Molar internal energy. Units: J/kmol. 
    virtual doublereal intEnergy_mole() const;

    /// Molar entropy. Units: J/kmol/K. 
    virtual doublereal entropy_mole() const;

    /// Molar Gibbs function. Units: J/kmol. 
    virtual doublereal gibbs_mole() const;

      /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K. 
    virtual doublereal cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    /*!
     * This method calculates the current pressure consistent with the
     * independent variables, T, rho.
     */
    virtual doublereal pressure() const;

    //! sets the thermodynamic pressure (Pa).
    /*!
     * This method calculates the density that is consistent with the
     * desired pressure, given the temperature.
     *
     * @param p  Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the 
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical 
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const {
      mu[0] = gibbs_mole();
    }

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;
    
    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!      
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;

    //! Returns a reference to the substance object
    tpx::Substance& TPX_Substance();

    /// critical temperature 
    virtual doublereal critTemperature() const;
 
    /// critical pressure
    virtual doublereal critPressure() const;
        
    /// critical density
    virtual doublereal critDensity() const;
        
    /// saturation temperature
    /*!
     * @param p  Pressure (Pa)
     */
    virtual doublereal satTemperature(doublereal p) const;

    //! Set the internally storred specific enthalpy (J/kg) and pressure (Pa) of the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_HP(doublereal h, doublereal p, 
			     doublereal tol = 1.e-8);

    //! Set the specific internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific internal energy and specific volume have the value of the input parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_UV(doublereal u, doublereal v, 
			     doublereal tol = 1.e-8);

    //! Set the specific entropy (J/kg/K) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and specific volume have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_SV(doublereal s, doublereal v, 
			     doublereal tol = 1.e-8);

    //! Set the specific entropy (J/kg/K) and pressure (Pa).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and the pressure have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param p    specific pressure (Pa).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_SP(doublereal s, doublereal p, 
			     doublereal tol = 1.e-8);


      
    //! @name Saturation properties.
    /*!
     * These methods are only implemented by subclasses that 
     * implement full liquid-vapor equations of state. They may be
     * moved out of ThermoPhase at a later date.
     */
    //@{

    //! Return the saturation pressure given the temperatur
    /*!
     * @param t Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal t) const;

    //! Return the fraction of vapor at the current conditions
    virtual doublereal vaporFraction() const;

    //! Set the state to a saturated system at a particular temperature
    /*!
     * @param t  Temperature (kelvin)
     * @param x  Fraction of vapor
     */
    virtual void setState_Tsat(doublereal t, doublereal x);

    //! Set the state to a saturated system at a particular pressure
    /*!
     * @param p  Pressure (Pa)
     * @param x  Fraction of vapor
     */
    virtual void setState_Psat(doublereal p, doublereal x);
    //@}

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //! Set equation of state parameter values from XML entries.
    /*!
     *
     * This method is called by function importPhase() in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialzed with elements and/or species.
     *   
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);


    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     */
    virtual std::string report(bool show_thermo = true) const;

  protected:

    //! Main call to the tpx level to set the state of the system
    /*!
     * @param n  Integer indicating which 2 thermo components are held constant
     * @param x  Value of the first component
     * @param y  Value of the second component
     */
    void Set(int n, double x, double y) const;

    //! Sets the state using a TPX::TV call
    void setTPXState() const;

    //! Carry out a internal check on tpx, it may have thrown an error.
    /*!
     * @param v Defaults to zero
     */
    void check(doublereal v = 0.0) const;

    //! Report errors in the TPX level
    void reportTPXError() const;

  private:

    //! Pointer to the underlying tpx object Substance that does the work
    mutable tpx::Substance* m_sub;

    //! Int indicating the type of the fluid
    /*!
     * The tpx package uses an int to indicate what fluid is being sought.
     */
    int m_subflag;

    //! Molecular weight of the substance (kg kmol-1)
    doublereal m_mw;

    //! flag to turn on some printing.
    bool m_verbose;
  };

}

#endif
#endif

