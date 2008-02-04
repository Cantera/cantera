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
 *  $Id$
 */

#ifndef CT_PDSS_H
#define CT_PDSS_H
#include "ct_defs.h"

class XML_Node;
class ThermoPhase;

 class WaterPropsIAPWS;


namespace Cantera {

  class SpeciesThermo;
 
  //! Virtual base class for a species with a pressure dependent 
  //! standard state 
  /*!
   * Virtual base class for calculation of  the
   * pressure dependent standard state for a single species
   * 
   * Class %PDSS is the base class
   * for a family of classes that compute properties of a set of
   * species in their standard states at a range of temperatures
   * and pressures. The independent variables for this object
   * are temperature and pressure.
   * The class may mave a reference to a SpeciesThermo object
   * which handles the calculation of the reference state temperature
   * behavior of a subset of species.
   *
   * This class is analagous to the SpeciesThermoInterpType
   * class, except that the standard state inherently incorporates
   * the pressure dependence.
   *
   * The class operates on a setState temperature and pressure basis.
   * It only recalculates the standard state when the setState functions
   * for temperature and pressure are called. 
   *
   */
  class PDSS {

  public:

    /**
     * Empty Constructor
     */
    PDSS();

    //! Constructor that initializes the object by examining the XML entries
    //! from the ThermoPhase object
    /*!
     *  This function calls the constructPDSS member function.
     * 
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS(ThermoPhase *tp, int spindex);

    //! Copy Constructor
    /*!
     * @param b object to be copied
     */
    PDSS(const PDSS &b);

    //! Assignment operator
    /*!
     * @param b Object to be copied
     */
    PDSS& operator=(const PDSS&b);

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSFile member function.
     * 
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param inputFile String name of the input file
     *  @param id        String name of the phase in the input file. The default
     *                   is the empty string, in which case the first phase in the
     *                   file is used.
     */
    PDSS(ThermoPhase *tp, int spindex, std::string inputFile, std::string id = "");

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSXML member function.
     * 
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param phaseRef  Reference to the XML tree containing the phase information.
     *  @param id        String name of the phase in the input file. The default
     *                   is the empty string, in which case the first phase in the
     *                   file is used.
     */
    PDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRef, std::string id = "");

    //! Destructor for the phase
    virtual ~PDSS();

    //! Duplication routine for objects which inherit from %PDSS
    /*!
     *  This virtual routine can be used to duplicate %PDSS  objects
     *  inherited from %PDSS even if the application only has
     *  a pointer to %PDSS to work with.
     *
     * @return returns a pointer to the base %PDSS object type
     */
    virtual PDSS *duplMyselfAsPDSS() const;
        
    /**
     *   
     * @name  Utilities
     * @{
     */

    //! Returns the type of the standard state parameterization
    /*!
     * @return Returns the integer # of the parameterization
     */
    virtual int pdssType() const { return -1; }

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Species Standard State
     * @{
     */
   
    //! Return the molar enthalpy in units of J kmol-1
    /*!
     * Returns the species standard state enthalpy in J kmol-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state enthalpy in  J kmol-1
     */
    virtual doublereal enthalpy_mole() const;

    //! Return the molar internal Energy in units of J kmol-1
    /*!
     * Returns the species standard state internal Energy in J kmol-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state internal Energy in  J kmol-1
     */
    virtual doublereal intEnergy_mole() const;

    //! Return the molar entropy in units of J kmol-1 K-1
    /*!
     * Returns the species standard state entropy in J kmol-1 K-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state entropy in J kmol-1 K-1
     */
    virtual doublereal entropy_mole() const;

    //! Return the molar gibbs free energy in units of J kmol-1
    /*!
     * Returns the species standard state gibbs free energy in J kmol-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state gibbs free energy in  J kmol-1
     */
    virtual doublereal gibbs_mole() const;

    //! Return the molar const pressure heat capacity in units of J kmol-1 K-1
    /*!
     * Returns the species standard state Cp in J kmol-1 K-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state Cp in J kmol-1 K-1
     */
    virtual doublereal cp_mole() const;

    //! Return the molar const volume heat capacity in units of J kmol-1 K-1
    /*!
     * Returns the species standard state Cv in J kmol-1 K-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state Cv in J kmol-1 K-1
     */
    virtual doublereal cv_mole() const;

    /*
     * Get the difference in the standard state thermodynamic properties
     * between the current pressure. and the reference pressure, p0
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
     * Molecular Weight of the species
     */
    doublereal m_mw;

    /**
     * Species index in the thermophase corresponding to this species.
     */
    int m_spindex;

    /**
     * Pointer to the species thermodynamic property manager.
     * This is a copy of the pointer in the ThermoPhase object.
     * Note, this object doesn't own the pointer.
     * If the SpeciesThermo ThermoPhase object doesn't know 
     * or doesn't control the calculation, this will be 
     * set to zero.
     */
    SpeciesThermo* m_spthermo;

    doublereal *m_cp0_R_ptr;
    doublereal *m_h0_RT_ptr;
    doublereal *m_s0_R_ptr;
    doublereal *m_g0_RT_ptr;


  };

}

#endif



