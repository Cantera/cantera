/**
 *  @file PDSS.h
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
 *  $Id: PDSS.h,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */

#ifndef CT_PDSS_H
#define CT_PDSS_H
#include "ct_defs.h"

class XML_Node;
class ThermoPhase;

class WaterPropsIAPWS;

namespace Cantera {


  /**
   * Class for pressure dependent standard states.
   *
   *
   */
  class PDSS {

  public:

    /**
     * Basic list of constructors and duplicators
     */
    PDSS(ThermoPhase *tp, int spindex);
    PDSS(const PDSS &b);
    PDSS& operator=(const PDSS&b);
    PDSS(ThermoPhase *tp, int spindex, string inputFile, string id = "");
    PDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRef, string id = "");
    virtual ~PDSS();
        
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
				   string inputFile, string id);
    virtual void constructPDSSXML(ThermoPhase *tp, int spindex, 
				  XML_Node& phaseNode, string id);
    virtual void initThermoXML(XML_Node& eosdata, string id);
    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

  protected:

    /**
     * state of the system (temperature and density);
     * This may redundant and may go away. Should be able to
     * get this information from owning ThermoPhase object.
     */
    mutable doublereal m_temp;

    /**
     * state of the system (temperature and density);
     * This may redundant and may go away. Should be able to
     * get this information from owning ThermoPhase object.
     */
    doublereal m_dens;

    /**
     * Thermophase which this species belongs to. Note, in some
     * applications (i.e., mostly testing applications, this may be a null
     * value. Applications should test whether this is null before usage. 
     */
    ThermoPhase *m_tp;

    /**
     * Species index in the thermophase corresponding to this species.
     */
    int m_spindex;

    /**
     * Molecular Weight of the species
     */
    doublereal m_mw;

  };

}

#endif



