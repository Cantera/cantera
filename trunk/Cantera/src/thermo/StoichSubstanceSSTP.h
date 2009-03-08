/**
 *
 *  @file StoichSubstanceSSTP.h
 * 
 * Header file for the StoichSubstanceSSTP class
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

/*  $Author: hkmoffa $
 *  $Date: 2005/11/14 18:49:56 $
 *  $Revision: 1.2 $
 *
 */

#ifndef CT_STOICHSUBSTANCESSTP_H
#define CT_STOICHSUBSTANCESSTP_H

#include "mix_defs.h"
#include "SingleSpeciesTP.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * @ingroup thermoprops
     *
     * Class StoichSubstance represents a stoichiometric (fixed composition) 
     * incompressible substance.
     *
     */
    class StoichSubstanceSSTP : public SingleSpeciesTP {

    public:
	/**
	 * Default Constructor for the StoichSubstanceSSTP class
	 */
        StoichSubstanceSSTP();

	/**
	 * Destructor for the routine (virtual)
	 *        
	 */
        virtual ~StoichSubstanceSSTP();

	/**
         *   
         * @name  Utilities  
         * @{
         */

        /**
         * Equation of state flag.
	 *
	 * Returns the value cStoichSubstance, defined in mix_defs.h.
         */
        virtual int eosType() const;

        /**
	 *  @}
         *  @name Molar Thermodynamic Properties of the Solution
         *  @{
         */

        /**
	 * @}
         * @name Mechanical Equation of State
         * @{
         */

        /**
         * Pressure. Units: Pa.
         * For an incompressible substance, the density is independent
         * of pressure. This method simply returns the stored
         * pressure value.
         */ 
        virtual doublereal pressure() const;

        /**
         * Set the pressure at constant temperature. Units: Pa.
         * For an incompressible substance, the density is 
         * independent of pressure. Therefore, this method only 
         * stores the specified pressure value. It does not 
         * modify the density.
         */
        virtual void setPressure(doublereal p);

        /**
         * The isothermal compressibility. Units: 1/Pa.
         * The isothermal compressibility is defined as
         * \f[
         * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
         * \f]
         */
        virtual doublereal isothermalCompressibility() const; 

        /**
         * The thermal expansion coefficient. Units: 1/K.
         * The thermal expansion coefficient is defined as
         *
         * \f[
         * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
         * \f]
         */
        virtual doublereal thermalExpansionCoeff() const ;


	/**
	 * @}
         * @name Activities, Standard States, and Activity Concentrations
	 *
	 *  This section is largely handled by parent classes, since there
	 *  is only one species. Therefore, the activity is equal to one.
	 * @{
         */

        /**
         * This method returns the array of generalized
         * concentrations.  For a stoichiomeetric substance, there is
         * only one species, and the generalized concentration is 1.0.
         */
        virtual void getActivityConcentrations(doublereal* c) const;

        /**
         * The standard concentration. This is defined as the concentration 
         * by which the generalized concentration is normalized to produce 
         * the activity. 
         */ 
	virtual doublereal standardConcentration(int k=0) const;

        /**
	 * Returns the natural logarithm of the standard 
	 * concentration of the kth species
	 */
	virtual doublereal logStandardConc(int k=0) const;

	/**
	 * Get the array of chemical potentials at unit activity 
         * \f$ \mu^0_k \f$.
	 *
         * For a stoichiometric substance, there is no activity term in 
         * the chemical potential expression, and therefore the
         * standard chemical potential and the chemical potential
         * are both equal to the molar Gibbs function.
         */
        virtual void getStandardChemPotentials(doublereal* mu0) const;

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
	 *  uA[0] = kmol units - default  = 0
	 *  uA[1] = m    units - default  = 0
	 *  uA[2] = kg   units - default  = 0;
	 *  uA[3] = Pa(pressure) units - default = 0;
	 *  uA[4] = Temperature units - default = 0;
	 *  uA[5] = time units - default = 0
	 */
	virtual void getUnitsStandardConc(double *uA, int k = 0,
					  int sizeUA = 6);

	//@}
        /// @name  Partial Molar Properties of the Solution
	///
	///        These properties are handled by the parent class,
	///        SingleSpeciesTP
        //@{


	//@}
        /// @name  Properties of the Standard State of the Species in the Solution 
        //@{

	/**
         * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current 
	 * <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const;

        /**
         * Get the array of nondimensional Entropy functions for the
         * standard state species
         * at the current <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEntropy_R(doublereal* sr) const;

        /**
         * Get the nondimensional Gibbs functions for the species
         * at their standard states of solution at the current T and P
         * of the solution
         */
        virtual void getGibbs_RT(doublereal* grt) const;

        /**
         * Get the nondimensional Gibbs functions for the standard
         * state of the species at the current T and P.
         */
        virtual void getCp_R(doublereal* cpr) const;


	/**
	 * Molar internal energy. J/kmol.  For an incompressible,
	 * stoichiometric substance, the molar internal energy is
	 * independent of pressure. Since the thermodynamic properties
	 * are specified by giving the standard-state enthalpy, the
	 * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
	 * enthalpy to compute the molar internal energy.
	 */
	virtual void getIntEnergy_RT(doublereal* urt) const;

	//@}
        /// @name Thermodynamic Values for the Species Reference States
	//@{

	/**
	 *  Returns the vector of nondimensional
         *  internal Energies of the reference state at the current temperature
         *  of the solution and the reference pressure for each species.
         */
        virtual void getIntEnergy_RT_ref(doublereal *urt) const;

	/*
	 * ---- Critical State Properties
	 */
	/// Critical temperature (K).
        virtual doublereal critTemperature() const;
        /// Critical pressure (Pa).
        virtual doublereal critPressure() const;
        /// Critical density (kg/m3).
        virtual doublereal critDensity() const;

	/*
	 * ---- Saturation Properties
	 */
	virtual doublereal satTemperature(doublereal p) const;
        virtual doublereal satPressure(doublereal t) const;
        virtual doublereal vaporFraction() const;
	virtual void setState_Tsat(doublereal t, doublereal x);
        virtual void setState_Psat(doublereal p, doublereal x);

	/*
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
        virtual void initThermo();

	/*
	 * setParameters:
	 *
	 *   Generic routine that is used to set the parameters used
	 *   by this model.
	 *        C[0] = density of phase [ kg/m3 ]
	 */
        virtual void setParameters(int n, double *c);
	/*
	 * getParameters:
	 *
	 *   Generic routine that is used to get the parameters used
	 *   by this model.
	 *        n = 1
	 *        C[0] = density of phase [ kg/m3 ]
	 */
        virtual void getParameters(int &n, double * const c);

	/*
	 * Reads an xml data block for the parameters needed by this
	 * routine. eosdata points to the thermo block, and looks
	 * like this:
	 * 
	 *   <phase id="stoichsolid" >
	 *     <thermo model="StoichSubstance">
	 *         <density units="g/cm3">3.52</density>
	 *     </thermo>
	 *   </phase>
	 */
	virtual void setParametersFromXML(const XML_Node& eosdata);

    protected:

    };
    
}
        
#endif
