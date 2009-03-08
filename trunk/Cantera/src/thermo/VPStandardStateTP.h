/**
 *  @file VPStandardStateTP.h
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties. These include most of the
 * methods for calculating liquid electrolyte thermodynamics.
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2006/06/13 16:02:42 $
 *  $Revision: 1.3 $
 */

#ifndef CT_VPSTANDARDSTATETP_H
#define CT_VPSTANDARDSTATETP_H

#include "ThermoPhase.h"

namespace Cantera {

    class XML_Node;

    /**
     * @ingroup thermoprops
     *
     *  This is a filter class for ThermoPhase that implements
     *  a variable pressure standard state for ThermoPhase objects.
     *
     *  In addition support for the molality unit scale is provided.
     *
     *   Currently, it really is just a shell. The ThermoPhase object
     *   itself is based around the general concepts of
     *   VPStandardStateTP. Therefore, there really isn't much going
     *   on here.  However, this may change. The ThermoPhase object
     *   itself could change. Additionally, this object may revolve
     *   around the molality unit scale in the near future. We will
     *   have to see how things fare.
     */

    class VPStandardStateTP : public ThermoPhase {

    public:
        
        /// Constructor. 
        VPStandardStateTP();

	/// Copy Constructor.
	VPStandardStateTP(const VPStandardStateTP &);

	/// Assignment operator
	VPStandardStateTP& operator=(const VPStandardStateTP &);

        /// Destructor. 
        virtual ~VPStandardStateTP();

	/*
	 * Duplication routine
	 */
	virtual ThermoPhase *duplMyselfAsThermoPhase();

        /**
         *   
         * @name  Utilities  
         * @{
         */

        /** 
         * Equation of state type flag. The base class returns
         * zero. Subclasses should define this to return a unique
         * non-zero value. Constants defined for this purpose are
         * listed in mix_defs.h.
         */
        virtual int eosType() const { return 0; }


        /**
         * @} 
         * @name  Molar Thermodynamic Properties of the Solution
         * @{
         */

	/*
	 * These are handled by inherited objects. At this level,
	 * this pass-through routine doesn't add anything to the
	 * ThermoPhase description.
	 */


        /**
         * @}
         * @name Mechanical Properties
         * @{
         */

	/*
	 * These are handled by inherited objects. At this level,
	 * this pass-through routine doesn't add anything to the
	 * ThermoPhase description.
	 */

	/**
         * @} 
         * @name Electric Potential
         * 
         * The phase may be at some non-zero electrical
	 * potential. These methods set or get the value of the
	 * electric potential.
         * @{
	 */

	/*
	 * These are handled by inherited objects. At this level,
	 * this pass-through routine doesn't add anything to the
	 * ThermoPhase description.
	 */

        /**
         * @}
         * @name Activities and Activity Concentrations
         *
         * The activity \f$a_k\f$ of a species in solution is
         * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
         * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T)\f$ is
         * the chemical potential at unit activity, which depends only
         * on temperature.
	 * @{
         */

      
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
	virtual void getUnitsStandardConc(double *uA, int k = 0,
					  int sizeUA = 6);

        //@}
        /// @name  Partial Molar Properties of the Solution 
        //@{

	/**
         * Get the array of non-dimensional species chemical potentials
	 * These are partial molar Gibbs free energies.
         * \f$ \mu_k / \hat R T \f$.
	 * Units: unitless
	 *
	 * We close the loop on this function, here, calling
	 * getChemPotentials() and then dividing by RT.
         */
        virtual void getChemPotentials_RT(doublereal* mu) const;

  
        //@}
        /// @name  Properties of the Standard State of the Species in the Solution
        //@{

	/*
	 * These are handled by inherited objects. At this level,
	 * this pass-through routine doesn't add anything to the
	 * ThermoPhase description.
	 *
	 * However, we assume these methods exist for inherited objects.
	 * Therefore, we will bring the error routines up to this object
	 */

	/**
         * Get the array of chemical potentials at unit activity.
	 * These
         * are the standard state chemical potentials \f$ \mu^0_k(T,P)
         * \f$.. The values are evaluated at the current
         * temperature and pressure.
         */
        virtual void getStandardChemPotentials(doublereal* mu) const {
            err("getStandardChemPotentials");
        }

        /**
	 * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current
         * <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const {
            err("getEnthalpy_RT");
        }

        /**
	 * Get the array of nondimensional Enthalpy functions for the
         * standard state species
         * at the current <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEntropy_R(doublereal* sr) const {
            err("getEntropy_R");
        }

        /**
         * Get the nondimensional Gibbs functions for the species
         * at their standard states of solution at the current T and P
         * of the solution.
         */
        virtual void getGibbs_RT(doublereal* grt) const {
            err("getGibbs_RT");
        }

	/**
	 * Get the nondimensional Gibbs functions for the standard
         * state of the species at the current T and P.
	 */
        virtual void getPureGibbs(doublereal* gpure) const {
            err("getPureGibbs");
        }

	/**
	 *  Returns the vector of nondimensional
         *  internal Energies of the standard state at the current temperature
         *  and pressure of the solution for each species.
         */
        virtual void getIntEnergy_RT(doublereal *urt) const {
            err("getIntEnergy_RT");
        }

        /**
         * Get the nondimensional Heat Capacities at constant
         * pressure for the standard state of the species 
         * at the current T and P. 
         */
        virtual void getCp_R(doublereal* cpr) const {
            err("getCp_R");
        }

	/**
         * Get the molar volumes of each species in their standard
         * states at the current
         * <I>T</I> and <I>P</I> of the solution.
         * units = m^3 / kmol
         */
        virtual void getStandardVolumes(doublereal *vol) const {
            err("getStandardVolumes");
        }

       //@}
        /// @name Thermodynamic Values for the Species Reference States --------------------
        //@{

	/**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEnthalpy_RT_ref(doublereal *hrt) const;
     
        /**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getGibbs_RT_ref(doublereal *grt) const;
                   
        /**
         *  Returns the vector of the
         *  gibbs function of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         *  units = J/kmol
         */
        virtual void getGibbs_ref(doublereal *g) const;
      
        /**
         *  Returns the vector of nondimensional
         *  entropies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEntropy_R_ref(doublereal *er) const;
                 
        /**
         *  Returns the vector of nondimensional
         *  constant pressure heat capacities of the reference state
         *  at the current temperature of the solution
         *  and reference pressure for the species.
         */
        virtual void getCp_R_ref(doublereal *cprt) const;

        ///////////////////////////////////////////////////////
        //
        //  The methods below are not virtual, and should not
        //  be overloaded.
        //
        //////////////////////////////////////////////////////

        /**
         * @name Specific Properties
         * @{
         */


        /**
         * @name Setting the State
         *
         * These methods set all or part of the thermodynamic
         * state.
         * @{
         */

        //@}

        /**
         * @name Chemical Equilibrium
         * Chemical equilibrium.
         * @{
         */

        //@}

	
        /**
         * Set equation of state parameter values from XML
         * entries. This method is called by function importPhase in
         * file importCTML.cpp when processing a phase definition in
         * an input file. It should be overloaded in subclasses to set
         * any parameters that are specific to that particular phase
         * model. 
         *   
         * @param eosdata An XML_Node object corresponding to
         * the "thermo" entry for this phase in the input file.
         */
        virtual void setParametersFromXML(const XML_Node& eosdata) {}
  

        //---------------------------------------------------------
        /// @name Critical state properties.
        /// These methods are only implemented by some subclasses.
        
        //@{
             
        //@}
        
        /// @name Saturation properties.
        /// These methods are only implemented by subclasses that 
        /// implement full liquid-vapor equations of state.
        ///
   

        //@}

        /// The following methods are used in the process of constructing
        /// the phase and setting its parameters from a specification in an 
        /// input file. They are not normally used in application programs.
        /// To see how they are used, see files importCTML.cpp and 
        /// ThermoFactory.cpp.

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
        virtual void initThermo();

        /**
	 *   Import and initialize a ThermoPhase object
	 *
	 * @param phaseNode This object must be the phase node of a
	 *             complete XML tree
	 *             description of the phase, including all of the
	 *             species data. In other words while "phase" must
	 *             point to an XML phase object, it must have
	 *             sibling nodes "speciesData" that describe
	 *             the species in the phase.
	 * @param id   ID of the phase. If nonnull, a check is done
	 *             to see if phaseNode is pointing to the phase
	 *             with the correct id. 
	 */
        void initThermoXML(XML_Node& phaseNode, string id);

    private:
      void initLengths();

    protected:
	/*
	 * The last temperature at which the reference thermodynamic
	 * properties were calculated at.
	 */
	mutable doublereal    m_tlast;
	/**
	 * Vector containing the species reference enthalpies at T = m_tlast
	 */
        mutable vector_fp      m_h0_RT;

	/**
	 * Vector containing the species reference constant pressure
	 * heat capacities at T = m_tlast
	 */
        mutable vector_fp      m_cp0_R;

	/**
	 * Vector containing the species reference Gibbs functions
	 * at T = m_tlast
	 */
        mutable vector_fp      m_g0_RT;

	/**
	 * Vector containing the species reference entropies
	 * at T = m_tlast
	 */
        mutable vector_fp      m_s0_R;
  
    private:

	/**
	 * VPStandardStateTP has its own err routine
	 *
	 */
        doublereal err(string msg) const;

	/**
	 * This function gets called for every call to functions in this
	 * class. It checks to see whether the temperature has changed and
	 * thus the reference thermodynamics functions for all of the species
	 * must be recalculated.
	 * If the temperature has changed, the species thermo manager is called
	 * to recalculate G, Cp, H, and S at the current temperature.
	 */
        void _updateRefStateThermo() const;
    };

}
        
#endif





