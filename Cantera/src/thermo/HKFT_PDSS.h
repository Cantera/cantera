/**
 *  @file HKFT_PDSS.h
 *
 * Declares class PDSS pressure dependent standard state
 * for a single species
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 * 
 */

#ifndef CT_HKFT_PDSS_H
#define CT_HKFT_PDSS_H
#include "ct_defs.h"

class XML_Node;
class ThermoPhase;

class WaterPropsIAPWS;
#include "PDSS.h"



namespace Cantera {


  class WaterProps;
  class WaterPDSS;

  /**
   * Class for pressure dependent standard states corresponding to 
   * ionic solutes in electrolyte water.
   *
   * NOTE: This is largely not done or not complete. 
   *
   */
  class HKFT_PDSS : public PDSS {

  public:

    /**
     * Basic list of constructors and duplicators
     */
    HKFT_PDSS(ThermoPhase *tp, int spindex);
    HKFT_PDSS(const HKFT_PDSS &b);
    HKFT_PDSS& operator=(const HKFT_PDSS&b);
    HKFT_PDSS(ThermoPhase *tp, int spindex,
	      std::string inputFile, std::string id = "");
    HKFT_PDSS(ThermoPhase *tp, int spindex, 
	      XML_Node& phaseRef, std::string id = "");
    virtual ~HKFT_PDSS();
        
    /**
     *   
     * @name  Utilities  
     * @{
     */
    virtual int pdssType() const { return -1; }

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

    /*
     * Get the difference in the standard state thermodynamic properties
     * between the reference pressure, po, and the current pressure.
     */
    virtual doublereal enthalpyDelp_mole() const;
    virtual doublereal intEnergyDelp_mole() const;
    virtual doublereal entropyDelp_mole() const;
    virtual doublereal gibbsDelp_mole() const;
    virtual doublereal cpDelp_mole() const;
    virtual doublereal cvDelp_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties ---------------------
    //@{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal p);

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
    
    virtual void setDensity(double dens);
    double density() const;
    virtual void setTemperature(double temp);
    double temperature() const;
    virtual void setState_TP(double temp, double pres);

    doublereal molecularWeight() const;
    void setMolecularWeight(double mw);
    
    void constructHKFT_PDSS(ThermoPhase *tp, int spindex);
    void constructHKFT_PDSSFile(ThermoPhase *tp, int spindex, 
				   std::string inputFile, std::string id);
    void constructHKFT_PDSSXML(ThermoPhase *tp, int spindex, 
				  XML_Node& phaseNode, std::string id);
    virtual void initThermoXML(XML_Node& eosdata, std::string id);
    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

    double deltaG();
    double electrostatic_radii_calc();

  private:
    double ag(const double temp, const int ifunc = 0) const;
    double bg(const double temp, const int ifunc = 0) const;
    double g(const double temp, const double pres, const int ifunc = 0);
    double f(const double temp, const double pres, const int ifunc = 0);
    double gstar(const double temp, const double pres, const int ifunc = 0);

  protected:

  private:
    //!  Water standard state calculator
    /*!
     *  derived from the equation of state for water.
     */
    WaterPDSS *m_waterSS;

    //! Current value of the pressure for this object
    doublereal m_pres;

    //! density of standard-state water
    /*!
     * internal temporary variable
     */
    double m_densWaterSS;

    /**
     *  Pointer to the water property calculator
     */
    WaterProps *m_waterProps;


    //! Born coefficient for the current ion or species

    doublereal m_born_coeff_j;

    //! Electrostatic radii
    doublereal r_e_j;


    //! Value of deltaG at Tr and Pr    (cal gmol-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     */
    doublereal m_deltaG_tr_pr;

    //! Value of S_j at Tr and Pr    (cal gmol-1 K-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     */
    doublereal m_Entrop_tr_pr;

    //! a1 coefficient (cal gmol-1 bar-1)
    doublereal m_a1;

    //! a2 coefficient (cal gmol-1)
    doublereal m_a2;

    //! c1 coefficient (cal gmol-1 K-1)
    doublereal m_c1;
    
    //! c2 coefficient (cal K gmol-1)
    doublereal m_c2;

    //! a3 coefficient (cal K gmol-1 bar-1)
    doublereal m_a3;

    //! a4 coefficient (cal K gmol-1)
    doublereal m_a4;

    //! omega_pr_tr coefficient(cal gmol-1)
    doublereal m_omega_pr_tr;

    //! y = dZdT = 1/(esp*esp) desp/dT
    double m_Y_pr_tr;
    //double m_Y_pr_tr = -5.799E-5;
    double m_Z_pr_tr;
    //double m_Z_pr_tr = -0.0127803;
    //! Reference pressure is 1 atm in units of bar= 1.0132
    doublereal m_presR_bar;


    //! Charge of the ion
    doublereal m_charge_j;

    WaterProps *m_wprops;

  };

}

#endif



