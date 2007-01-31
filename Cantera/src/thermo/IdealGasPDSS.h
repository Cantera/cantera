/**
 *  @file IDEALGASPDSS.h
 *
 * Declares class PDSS pressure dependent standard state
 * for a single species
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id$
 */

#ifndef CT_IDEALGASPDSS_H
#define CT_IDEALGASPDSS_H

#include "PDSS.h"

class XML_Node;
class ThermoPhase;

namespace Cantera {


  /**
   * Class for pressure dependent standard states.
   * This class is for a single Ideal Gas species.
   *
   */
  class IdealGasPDSS : public PDSS {

  public:

    /**
     * Basic list of constructors and duplicators
     */
    IdealGasPDSS(ThermoPhase *tp, int spindex);
    IdealGasPDSS(const IdealGasPDSS &b);
    IdealGasPDSS& operator=(const IdealGasPDSS&b);
    IdealGasPDSS(ThermoPhase *tp, int spindex, std::string inputFile, std::string id = "");
    IdealGasPDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRef, std::string id = "");
    virtual ~IdealGasPDSS();
        
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
    
    virtual void constructPDSS(ThermoPhase *tp, int spindex);
    virtual void constructPDSSFile(ThermoPhase *tp, int spindex, 
				   std::string inputFile, std::string id);
    virtual void constructPDSSXML(ThermoPhase *tp, int spindex, 
				  XML_Node& phaseNode, std::string id);
    virtual void initThermoXML(XML_Node& eosdata, std::string id);
    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

  protected:

    int m_kk, m_mm;
    doublereal m_tmin, m_tmax, m_p0;
    
  

  };

}

#endif



