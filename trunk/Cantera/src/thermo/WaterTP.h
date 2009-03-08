/**
 *  @file WaterTP.h
 *
 * Declares a ThermoPhase class consisting of 
 * pure water.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterTP.h,v 1.1 2006/07/04 00:01:54 hkmoffa Exp $
 */

#ifndef CT_WATERTP_H
#define CT_WATERTP_H

#include "ThermoPhase.h"

class WaterPropsIAPWS;

namespace Cantera {


    /**
     *Class for single-component water. This is designed to cover just the
     * liquid part of water.
     *
     *
     * Notes:
     *   Base state for thermodynamic properties:
     * 
     *   The thermodynamic base state for water is set to the NIST basis here
     *   by specifying constants EW_Offset and SW_Offset. These offsets are
     *   specified so that the following properties hold:
     *
     *   Delta_Hfo_gas(298.15) = -241.826 kJ/gmol
     *   So_gas(298.15, 1bar)  = 188.835 J/gmolK
     *
     *           (http://webbook.nist.gov)
     *
     *   The "o" here refers to a hypothetical ideal gas state. The way
     *   we achieve this in practice is to evaluate at a very low pressure
     *   and then use the theoretical ideal gas results to scale up to
     *   higher pressures:
     *
     *   Ho(1bar) = H(P0)
     *
     *   So(1bar) = S(P0) + RT ln(1bar/P0)
     *
     *   The offsets used in the steam tables are different than NIST's. 
     *   They assume u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the
     *   triple point conditions.
     *
     *
     */
    class WaterTP : public ThermoPhase {

    public:

	/**
	 * Basic list of constructors and duplicators
	 */
        WaterTP(); 
        WaterTP(const WaterTP &b);
        WaterTP& operator=(const WaterTP&b);
        WaterTP(string inputFile, string id = "");
        WaterTP(XML_Node& phaseRef, string id = "");
        virtual ~WaterTP();
        ThermoPhase *duplMyselfAsThermoPhase();
        
       /**
         *   
         * @name  Utilities  
         * @{
         */
        virtual int eosType() const { return -1; }

       /**
         * @} 
         * @name  Molar Thermodynamic Properties of the Solution --------------
         * @{
         */
        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual doublereal entropy_mole() const;
        virtual doublereal gibbs_mole() const;
        virtual doublereal cp_mole() const;
        virtual doublereal cv_mole() const;

	//@}
        /// @name Mechanical Equation of State Properties ---------------------
        //@{

        virtual doublereal pressure() const;
        virtual void setPressure(doublereal p);

	/**
         * @} 
         * @name Potential Energy
         * @{
	 */

	/**
         * @}
         * @name Activities, Standard States,  and Activity Concentrations
	 * @{
	 */

	//@}
        /// @name  Partial Molar Properties of the Solution -----------------
        //@{

        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = gibbs_mole();
        }

	//@}
        /// @name  Properties of the Standard State of the Species
        //          in the Solution --
        //@{
    

        /// critical temperature 
        virtual doublereal critTemperature() const;
 
        /// critical pressure
        virtual doublereal critPressure() const;
        
        /// critical density
        virtual doublereal critDensity() const;
        
        /// saturation temperature
        //virtual doublereal satTemperature(doublereal p) const;
        
    

        /// saturation pressure
        virtual doublereal satPressure(doublereal t);
    

	virtual void setTemperature(double temp);
    
	virtual void constructPhase();
	virtual void constructPhaseFile(string inputFile, string id);
	virtual void constructPhaseXML(XML_Node& phaseNode, string id);


	virtual void initThermoXML(XML_Node& eosdata, string id);
        virtual void initThermo();
        virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
        
        void Set(int n, double x, double y) const;
        void setTPXState() const;
        void check(doublereal v = 0.0) const;
        void reportTPXError() const;

private:
	mutable WaterPropsIAPWS *m_sub;
        int m_subflag;
        doublereal m_mw;

	/**
	 * Offset constants used to obtain consistency with the NIST database.
	 * This is added to all internal energy and enthalpy results.
	 *  units = J kmol-1.
	 */
	double EW_Offset;

	/*
	 * Offset constant used to obtain consistency with NIST convention.
	 * This is added to all internal entropy results.
	 *  units = J kmol-1 K-1.
	 */
	double SW_Offset;

        bool m_verbose;

	/**
	 *  Since this phase represents a liquid phase, it's an error to 
	 *  return a gas-phase answer. However, if the below is true, then
	 *  a gas-phase answer is allowed. This is used to check the thermodynamic
	 *  consistency with ideal-gas thermo functions for example.
	 */
	bool m_allowGasPhase;
    };

}

#endif



