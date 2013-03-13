/**
 *  @file PDSS_IdealGas.h
 *   Declarations for the class PDSS_IdealGas (pressure dependent standard state)
 *    which handles calculations for a single ideal gas species in a phase
 *    (see \ref pdssthermo and class \link Cantera::PDSS_IdealGas PDSS_IdealGas\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_PDSS_IDEALGAS_H
#define CT_PDSS_IDEALGAS_H

#include "PDSS.h"


namespace Cantera
{
class XML_Node;
class VPStandardStateTP;


//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * This class is for a single Ideal Gas species.
 *
 * @ingroup pdssthermo
 */
class PDSS_IdealGas : public PDSS
{

public:

    /**
      * @name  Constructors
      * @{
      */

    //! Constructor
    /*!
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS_IdealGas(VPStandardStateTP* tp, int spindex);

    //! Copy Constructur
    /*!
     * @param b Object to be copied
     */
    PDSS_IdealGas(const PDSS_IdealGas& b);

    //! Assignment operator
    /*!
     * @param b Object to be copeid
     */
    PDSS_IdealGas& operator=(const PDSS_IdealGas& b);

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
    PDSS_IdealGas(VPStandardStateTP* tp, int spindex,
                  std::string inputFile, std::string id = "");


    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSXML member function.
     *
     *  @param vptp_ptr    Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex     Species index of the species in the phase
     *  @param speciesNode Reference to the species XML tree.
     *  @param phaseRef    Reference to the XML tree containing the phase information.
     *  @param spInstalled Boolean indicating whether the species is installed yet
     *                     or not.
     */
    PDSS_IdealGas(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
                  const XML_Node& phaseRef, bool spInstalled);


    //! Destructor
    virtual ~PDSS_IdealGas();

    //! Duplicator
    virtual PDSS* duplMyselfAsPDSS() const;

    /**
     * @}
     * @name  Utilities
     * @{
     */

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

    /**
     * @}
     * @name Properties of the Reference State of the Species
     *       in the Solution
     * @{
     */

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

    /*
     * Get the difference in the standard state thermodynamic properties
     * between the reference pressure, po, and the current pressure.
     */

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

    //! Set the internal temperature
    /*!
     * @param temp Temperature (Kelvin)
     */
    virtual void setTemperature(doublereal temp);

    //! Return the current stored temperature
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
     * @param  rho      Density (Pascals)
     */
    virtual void setState_TR(doublereal temp, doublereal rho);

    /**
     * @}
     *  @name  Miscellaneous properties of the standard state
     * @{
     */

    /// critical temperature
    virtual doublereal critTemperature() const;

    /// critical pressure
    virtual doublereal critPressure() const;

    /// critical density
    virtual doublereal critDensity() const;

    /// saturation pressure
    /*!
     *  @param t  Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal t);

    /**
     * @}
     *  @name  Initialization of the Object
     * @{
     */

    //! Initialization of a PDSS object using an
    //! input XML file.
    /*!
     *
     * This routine is a precursor to constructPDSSXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param vptp_ptr    Pointer to the Variable pressure %ThermoPhase object
     *                    This object must have already been malloced.
     *
     * @param spindex     Species index within the phase
     *
     * @param inputFile   XML file containing the description of the
     *                    phase
     *
     * @param id          Optional parameter identifying the name of the
     *                    phase. If none is given, the first XML
     *                    phase element will be used.
     */
    void constructPDSSFile(VPStandardStateTP* vptp_ptr, size_t spindex,
                           std::string inputFile, std::string id);

    //!Initialization of a PDSS object using an xml tree
    /*!
     * This routine is a driver for the initialization of the
     * object.
     *
     *   basic logic:
     *       initThermo()                 (cascade)
     *       getStuff from species Part of XML file
     *       initThermoXML(phaseNode)      (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure %ThermoPhase object
     *                   This object must have already been malloced.
     *
     * @param spindex    Species index within the phase
     *
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     *
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& phaseNode, std::string id);

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

    //@}



protected:

    //! Maximum temperature the standard states are good for
    doublereal m_tmin;

    //! Minimum temperature the standard states are good for
    doublereal m_tmax;

};
}

#endif



