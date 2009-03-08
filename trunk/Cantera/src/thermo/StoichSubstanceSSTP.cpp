/**
 *
 *  @file StoichSubstanceSSTP.cpp
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

/*
 * $Id: StoichSubstanceSSTP.cpp,v 1.1 2005/10/24 21:52:23 hkmoffa Exp $ 
 */

#include "ct_defs.h"
#include "mix_defs.h"
#include "StoichSubstanceSSTP.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /*
     * ----  Constructors -------
     */

    /**
     * Default Constructor for the StoichSubstanceSSTP class
     */
    StoichSubstanceSSTP::StoichSubstanceSSTP():
	SingleSpeciesTP()
    {
    }

    /**
     * Destructor for the routine (virtual)
     *        
     */
    StoichSubstanceSSTP::~StoichSubstanceSSTP() 
    {
    }

    /*
     * ---- Utilities -----
     */

    /**
     * Equation of state flag. Returns the value cStoichSubstance,
     * defined in mix_defs.h.
     */
    int StoichSubstanceSSTP::eosType() const {
	return cStoichSubstance; 
    }

    /*
     * ---- Molar Thermodynamic properties of the solution ----
     */

    /**
     * ----- Mechanical Equation of State ------
     */

    /**
     * Pressure. Units: Pa.
     * For an incompressible substance, the density is independent
     * of pressure. This method simply returns the stored
     * pressure value.
     */ 
    doublereal StoichSubstanceSSTP::pressure() const {
	return m_press;
    }
    
    /**
     * Set the pressure at constant temperature. Units: Pa.
     * For an incompressible substance, the density is 
     * independent of pressure. Therefore, this method only 
     * stores the specified pressure value. It does not 
     * modify the density.
     */
    void StoichSubstanceSSTP::setPressure(doublereal p) {
	m_press = p;
    }

    /**
     * The isothermal compressibility. Units: 1/Pa.
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *
     *  It's equal to zero for this model, since the molar volume
     *  doesn't change with pressure or temperature.
     */
    doublereal StoichSubstanceSSTP::isothermalCompressibility() const {
	return 0.0;
    } 
    
    /**
     * The thermal expansion coefficient. Units: 1/K.
     * The thermal expansion coefficient is defined as
     *
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     *
     *  It's equal to zero for this model, since the molar volume
     *  doesn't change with pressure or temperature.
     */
    doublereal StoichSubstanceSSTP::thermalExpansionCoeff() const {
	return 0.0;
    }
    
    /*
     * ---- Chemical Potentials and Activities ----
     */

    /**
     * This method returns the array of generalized
     * concentrations.  For a stoichiomeetric substance, there is
     * only one species, and the generalized concentration is 1.0.
     */
    void StoichSubstanceSSTP::
    getActivityConcentrations(doublereal* c) const {
	c[0] = 1.0;
    }

    /**
     * The standard concentration. This is defined as the concentration 
     * by which the generalized concentration is normalized to produce 
     * the activity. 
     */ 
    doublereal StoichSubstanceSSTP::standardConcentration(int k) const {
	return 1.0;
    }

    /**
     * Returns the natural logarithm of the standard 
     * concentration of the kth species
     */
    doublereal StoichSubstanceSSTP::logStandardConc(int k) const {
	return 0.0;
    }

    /**
     * Returns the units of the standard and generalized
     * concentrations Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     */
    void StoichSubstanceSSTP::
    getUnitsStandardConc(double *uA, int k, int sizeUA) {
	for (int i = 0; i < 6; i++) {
	  uA[i] = 0;
	}
    }

    /*
     *  ---- Partial Molar Properties of the Solution ----
     */



    /*
     * ---- Properties of the Standard State of the Species in the Solution
     * ----
     */

    /**
     * Get the array of chemical potentials at unit activity 
     * \f$ \mu^0_k \f$.
     *
     * For a stoichiometric substance, there is no activity term in 
     * the chemical potential expression, and therefore the
     * standard chemical potential and the chemical potential
     * are both equal to the molar Gibbs function.
     */
    void StoichSubstanceSSTP::
    getStandardChemPotentials(doublereal* mu0) const {
	getGibbs_RT(mu0);
	mu0[0] *= GasConstant * temperature();
    }
    
    /**
     * Get the nondimensional Enthalpy functions for the species
     * at their standard states at the current 
     * <I>T</I> and <I>P</I> of the solution.
     * Molar enthalpy. Units: J/kmol.  For an incompressible,
     * stoichiometric substance, the internal energy is
     * independent of pressure, and therefore the molar enthalpy
     * is \f[ \hat h(T, P) = \hat u(T) + P \hat v \f], where the
     * molar specific volume is constant.
     */
    void StoichSubstanceSSTP::getEnthalpy_RT(doublereal* hrt) const {
	getEnthalpy_RT_ref(hrt);
	double RT = GasConstant * temperature();
	double presCorrect = (m_press - m_p0) /  molarDensity();
	hrt[0] += presCorrect / RT;
    }

    /**
     * Get the array of nondimensional Entropy functions for the
     * standard state species
     * at the current <I>T</I> and <I>P</I> of the solution.
     */
    void StoichSubstanceSSTP::getEntropy_R(doublereal* sr) const {
	getEntropy_R_ref(sr);
    }

    /**
     * Get the nondimensional Gibbs functions for the species
     * at their standard states of solution at the current T and P
     * of the solution
     */
    void StoichSubstanceSSTP::getGibbs_RT(doublereal* grt) const {
	getEnthalpy_RT(grt);
	grt[0] -= m_s0_R[0];
    }

    /**
     * Get the nondimensional Gibbs functions for the standard
     * state of the species at the current T and P.
     */
    void StoichSubstanceSSTP::getCp_R(doublereal* cpr) const {
	_updateThermo();
	cpr[0] = m_cp0_R[0];
    }
   
    /**
     * Molar internal energy (J/kmol).
     * For an incompressible,
     * stoichiometric substance, the molar internal energy is
     * independent of pressure. Since the thermodynamic properties
     * are specified by giving the standard-state enthalpy, the
     * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
     * enthalpy to compute the molar internal energy.
     */
    void StoichSubstanceSSTP::getIntEnergy_RT(doublereal* urt) const {
	_updateThermo();
	double RT = GasConstant * temperature();
	double PV = m_press / molarDensity();
	urt[0] = m_h0_RT[0] - PV / RT;
    }

    /*
     * ---- Thermodynamic Values for the Species Reference States ----
     */
    /**
     * Molar internal energy or the reference state at the current
     * temperature, T  (J/kmol).  
     * For an incompressible,
     * stoichiometric substance, the molar internal energy is
     * independent of pressure. Since the thermodynamic properties
     * are specified by giving the standard-state enthalpy, the
     * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
     * enthalpy to compute the molar internal energy.
     *
     * Note, this is equal to the standard state internal energy
     * evaluated at the reference pressure.
     */
    void StoichSubstanceSSTP::getIntEnergy_RT_ref(doublereal* urt) const {
	_updateThermo();
	double RT = GasConstant * temperature();
	double PV = m_p0 / molarDensity();
	urt[0] = m_h0_RT[0] - PV / RT;
    }
 
    /*
     * ---- Critical State Properties
     */
    /// Critical temperature (K).
    doublereal StoichSubstanceSSTP::critTemperature() const {
	return -1.0;
    }
    
    /// Critical pressure (Pa).
   doublereal StoichSubstanceSSTP::critPressure() const {
       return -1.0;
    }
    
    /// Critical density (kg/m3).
   doublereal StoichSubstanceSSTP::critDensity() const {
       return -1.0;
    }
             
    /*
     * ---- Saturation Properties
     */
    
    doublereal StoichSubstanceSSTP::satTemperature(doublereal p) const {
	return (-1.0);
    }       
    doublereal StoichSubstanceSSTP::satPressure(doublereal t) const {
	return 0.0;
    }
    doublereal StoichSubstanceSSTP::vaporFraction() const {
	return 0.0;
    }
    void StoichSubstanceSSTP::setState_Tsat(doublereal t, doublereal x) {
	setTemperature(t);
    }
    void StoichSubstanceSSTP::setState_Psat(doublereal p, doublereal x) {
	setPressure(p);
    }

    /*
     * ---- Initialization and Internal functions
     */

    /**
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    void StoichSubstanceSSTP::initThermo() {
	/*
	 * Make sure there is one and only one species in this phase.
	 */
        m_kk = nSpecies();
        if (m_kk != 1) {
            throw CanteraError("initThermo",
                "stoichiometric substances may only contain one species.");
        }
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
	/*
	 * Store the reference pressure in the variables for the class.
	 */
        m_p0 = refPressure();

	/*
	 * Resize temporary arrays.
	 */
        int leng = 1;
        m_h0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
	/*
	 * Call the base class thermo initializer
	 */
	SingleSpeciesTP::initThermo();
    }

    /**
     * setParameters:
     *
     *   Generic routine that is used to set the parameters used
     *   by this model.
     *        C[0] = density of phase [ kg/m3 ]
     */
    void StoichSubstanceSSTP::setParameters(int n, double * c) {
        double rho = c[0];
        setDensity(rho);
    }

    /**
     * getParameters:
     *
     *   Generic routine that is used to get the parameters used
     *   by this model.
     *        n = 1
     *        C[0] = density of phase [ kg/m3 ]
     */
    void StoichSubstanceSSTP::getParameters(int &n, double * const c) {
        double rho = density();
	n = 1;
        c[0] = rho;
    }

    /**
     * Reads an xml data block for the parameters needed by this
     * routine. eosdata is a reference to the xml thermo block, and looks
     * like this:
     * 
     *   <phase id="stoichsolid" >
     *     <thermo model="StoichSubstance">
     *         <density units="g/cm3">3.52</density>
     *     </thermo>
     *   </phase>
     */
    void StoichSubstanceSSTP::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","StoichSubstanceSSTP");
        doublereal rho = getFloat(eosdata, "density", "-");
        setDensity(rho);
    }

}


