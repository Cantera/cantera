/**
 *  @file PDSS.h
 *    Declarations for the virtual base class PDSS (pressure dependent standard state)
 *    which handles calculations for a single species in a phase
 *    (see class \link Cantera::PDSS PDSS\endlink).
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
#include "mix_defs.h"

class WaterPropsIAPWS;


namespace Cantera {
 
  class XML_Node;
  class SpeciesThermo;
  class VPStandardStateTP;
  class VPSSMgr;
 
  //! Virtual base class for a species with a pressure dependent 
  //! standard state 
  /*!
   * Virtual base class for calculation of the
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
   * <H3> Thread Safety </H3>
   *
   *   These classes are designed such that they are not thread safe when
   *   called by themselves. The reason for this is that they sometimes use
   *   shared SpeciesThermo resources where they set the states. This condition
   *   may be remedied in the future if we get serious about employing 
   *   multithreaded capabilities by adding mutex locks to the 
   *   SpeciesThermo resources. 
   *
   *   However, in many other respects they can be thread safe. They use
   *   separate memory and hold intermediate data. 
   */
  class PDSS {

  public:

    /**
     * @name  Constructors
     * @{
     */

   
    //! Empty Constructor
    PDSS();

    //! Constructor that initializes the object by examining the XML entries
    //! from the ThermoPhase object
    /*!
     *  This function calls the constructPDSS member function.
     * 
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS(VPStandardStateTP *tp, int spindex);

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
     * @}   
     * @name  Utilities
     * @{
     */

    //! Returns the type of the standard state parameterization
    /*!
     * @return Returns the integer # of the parameterization
     */
    virtual PDSS_enumType reportPDSSType() const { return cPDSS_UNDEF; }

  private:

    //! Set an error within this object for an unhandled capability
    /*!
     * @param msg    Message string for this error
     */
    void err(std::string msg) const;
    
  public:

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Species Standard State
     *        in the Solution
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

    //! Return the standard state molar enthalpy divided by RT
    /*!
     * Returns the species standard state enthalpy divided by RT at the
     * current temperature and pressure.
     *
     * @return returns the species standard state enthalpy in unitless form
     */
    virtual doublereal enthalpy_RT() const;

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

    //! Return the standard state entropy divided by RT
    /*!
     * Returns the species standard state entropy divided by RT at the
     * current temperature and pressure.
     *
     * @return returns the species standard state entropy divided by RT
     */
    virtual doublereal entropy_R() const;

    //! Return the molar gibbs free energy in units of J kmol-1
    /*!
     * Returns the species standard state gibbs free energy in J kmol-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state gibbs free energy in  J kmol-1
     */
    virtual doublereal gibbs_mole() const;

    //! Return the molar gibbs free energy divided by RT
    /*!
     * Returns the species standard state gibbs free energy divided by RT at the
     * current temperature and pressure.
     *
     * @return returns the species standard state gibbs free energy divided by RT
     */
    virtual doublereal gibbs_RT() const;

    //! Return the molar const pressure heat capacity in units of J kmol-1 K-1
    /*!
     * Returns the species standard state Cp in J kmol-1 K-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state Cp in J kmol-1 K-1
     */
    virtual doublereal cp_mole() const;

    //! Return the molar const pressure heat capacity divided by RT
    /*!
     * Returns the species standard state Cp divided by RT at the
     * current temperature and pressure.
     *
     * @return returns the species standard state Cp divided by RT
     */
    virtual doublereal cp_R() const;

    //! Return the molar const volume heat capacity in units of J kmol-1 K-1
    /*!
     * Returns the species standard state Cv in J kmol-1 K-1 at the
     * current temperature and pressure.
     *
     * @return returns the species standard state Cv in J kmol-1 K-1
     */
    virtual doublereal cv_mole() const;

    //! Return the molar volume at standard state
    /*!
     * Returns the species standard state molar volume at the
     * current temperature and pressure
     *
     * @return returns the standard state molar volume divided by R
     *             units are m**3 kmol-1.
     */
    virtual doublereal molarVolume() const;

    //! Return the standard state density at standard state
    /*!
     * Returns the species standard state density at the
     * current temperature and pressure
     *
     * @return returns the standard state density
     *             units are kg m-3
     */
    virtual doublereal density() const;
 
    //! Get the difference in the standard state enthalpy
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal enthalpyDelp_mole() const;

    //! Get the difference in the standard state entropy between
    //! the current pressure and the reference pressure, p0
    virtual doublereal entropyDelp_mole() const;

    //! Get the difference in the standard state gibbs free energy
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal gibbsDelp_mole() const;

    //! Get the difference in standard state heat capacity
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal cpDelp_mole() const;

    /**
     * @} 
     * @name Properties of the Reference State of the Species
     *       in the Solution 
     * @{
     */

    //! Return the reference pressure for this phase.
    doublereal refPressure() const {
      return m_p0;
    }

    //! Return the molar gibbs free energy divided by RT at reference pressure
    /*!
     * Returns the species reference state gibbs free energy divided by RT at the
     * current temperature.
     *
     * @return returns the reference state gibbs free energy divided by RT
     */
    virtual doublereal gibbs_RT_ref() const;

    //! Return the molar enthalpy divided by RT at reference pressure
    /*!
     * Returns the species reference state enthalpy divided by RT at the
     * current temperature.
     *
     * @return returns the reference state enthalpy divided by RT
     */
    virtual doublereal enthalpy_RT_ref() const;

    //! Return the molar entropy divided by R at reference pressure
    /*!
     * Returns the species reference state entropy divided by R at the
     * current temperature.
     *
     * @return returns the reference state entropy divided by R
     */
    virtual doublereal entropy_R_ref() const;

    //! Return the molar heat capacity divided by R at reference pressure
    /*!
     * Returns the species reference state heat capacity divided by R at the
     * current temperature.
     *
     * @return returns the reference state heat capacity divided by R
     */
    virtual doublereal cp_R_ref() const;

    //! Return the molar volume at reference pressure
    /*!
     * Returns the species reference state molar volume at the
     * current temperature.
     *
     * @return returns the reference state molar volume divided by R
     *             units are m**3 kmol-1.
     */
    virtual doublereal molarVolume_ref() const;

    /**
     * @}
     *  @name Mechanical Equation of State Properties 
     * @{
     */

    //! Returns the pressure (Pa)
    virtual doublereal pressure() const;

    //! Sets the pressure in the object
    /*!
     * Currently, this sets the pressure in the PDSS object.
     * It is indeterminant what happens to the owning VPStandardStateTP
     * object and to the VPSSMgr object.
     *
     * @param   pres   Pressure to be set (Pascal)
     */
    virtual void setPressure(doublereal pres);

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;
 
    //@}
    /// @name  Partial Molar Properties of the Solution -----------------
    //@{

    //! Set the internal temperature
    /*!
     * @param temp Temperature (Kelvin)
     */
    virtual void setTemperature(doublereal temp);

    //! Return the current storred temperature
    doublereal temperature() const;
 
    //! Set the internal temperature and pressure
    /*!
     * @param  temp     Temperature (Kelvin)
     * @param  pres     pressure (Pascals)
     */
    virtual void setState_TP(doublereal temp, doublereal pres);

    //! Set the internal temperature and density
    /*! 
     * @param  temp     Temperature (Kelvin)
     * @param  rho      Density (kg m-3)
     */
    virtual void setState_TR(doublereal temp, doublereal rho);
  
    /**
     * @}
     *  @name  Miscellaneous properties of the standard state
     * @{
     */
  
    //! critical temperature 
    virtual doublereal critTemperature() const;
 
    //! critical pressure
    virtual doublereal critPressure() const;
        
    //! critical density
    virtual doublereal critDensity() const;
        
    //! saturation pressure
    /*!
     *  @param T Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal T);
    
    //! Return the molecular weight of the species
    //! in units of kg kmol-1
    doublereal molecularWeight() const;

    //! Set the molecular weight of the species
    /*!
     * @param mw Molecular Weight in kg kmol-1
     */
    void setMolecularWeight(doublereal mw);

    /**
     * @}
     *  @name  Initialization of the Object
     * @{
     */


    //! Initialization routine for all of the shallow pointers
    /*!
     *  This is a cascading call, where each level should call the
     *  the parent level.
     *
     *  The initThermo() routines get called before the initThermoXML() routines
     *  from the constructPDSSXML() routine.
     *
     *
     *  Calls initPtrs();
     */
    virtual void initThermo();

    //! Initialization routine for the PDSS object based on the phaseNode
    /*!
     *  This is a cascading call, where each level should call the
     *  the parent level.
     *
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     *
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    virtual void initThermoXML(const XML_Node& phaseNode, std::string& id);

  
    //! Initialize all of the internal shallow pointers that can be initialized
    /*!
     * This routine isn't virtual
     */
    void initPtrs();

    //! Initialize or Reinitialize all shallow pointers in the object
    /*!
     *  This command is called to reinitialize all shallow pointers in the
     *  object. It's needed for the duplicator capability
     *
     * @param vptp_ptr       Pointer to the Variable pressure %ThermoPhase object
     *                       This object must have already been malloced.
     *
     * @param vpssmgr_ptr    Pointer to the variable pressure standard state
     *                       calculator for this phase
     *
     * @param spthermo_ptr   Pointer to the optional SpeciesThermo object
     *                       that will handle the calculation of the reference
     *                       state thermodynamic coefficients.
     */
    void initAllPtrs(VPStandardStateTP *vptp_ptr, VPSSMgr *vpssmgr_ptr, 
		     SpeciesThermo* spthermo_ptr);

   //@}

  protected:

    //! Enumerated type describing the type of the PDSS object
    PDSS_enumType m_pdssType;

    //! Current temperature used by the PDSS object
    mutable doublereal m_temp;

    //! State of the system - pressure
    mutable doublereal m_pres;

    //! reference state pressure of the species.
    doublereal m_p0;

    //! Thermophase which this species belongs to. 
    /*!
     * Note, in some
     * applications (i.e., mostly testing applications, this may be a null
     * value. Applications should test whether this is null before usage. 
     */
    VPStandardStateTP *m_tp;

    //! Pointer to the VPSS manager for this object
    VPSSMgr *m_vpssmgr_ptr;

    /**
     * Molecular Weight of the species
     */
    doublereal m_mw;

    /**
     * Species index in the thermophase corresponding to this species.
     */
    int m_spindex;

    //! Pointer to the species thermodynamic property manager.
    /*!
     * This is a copy of the pointer in the ThermoPhase object.
     * Note, this object doesn't own the pointer.
     * If the SpeciesThermo ThermoPhase object doesn't know 
     * or doesn't control the calculation, this will be 
     * set to zero.
     */
    SpeciesThermo* m_spthermo;
  
    //!  Reference state enthalpy divided by RT.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and m_p0
     */
    doublereal *m_h0_RT_ptr;

    //!  Reference state heat capacity divided by R.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and m_p0
     */
    doublereal *m_cp0_R_ptr;

    //!  Reference state entropy divided by R.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and m_p0
     */
    doublereal *m_s0_R_ptr;

    //!  Reference state gibbs free energy divided by RT.
    /*!
     *  Calculated at the current value of T and m_p0
     */
    doublereal *m_g0_RT_ptr;

    //!  Reference state molar volume (m3 kg-1)
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and m_p0
     */
    doublereal *m_V0_ptr;

    //!  Standard state enthalpy divided by RT.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and P
     */
    doublereal *m_hss_RT_ptr;

    //!  Standard state heat capacity divided by R.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and P
     */
    doublereal *m_cpss_R_ptr;

    //!  Standard state entropy divided by R.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and P
     */
    doublereal *m_sss_R_ptr;

    //!  Standard state gibbs free energy divided by RT.
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and P
     */
    doublereal *m_gss_RT_ptr;

    //!  Standard State molar volume (m3 kg-1)
    /*!
     *  Storage for the thermo properties is provided by
     *  VPSSMgr.
     *  Calculated at the current value of T and P
     */
    doublereal *m_Vss_ptr;

  };

}

#endif



