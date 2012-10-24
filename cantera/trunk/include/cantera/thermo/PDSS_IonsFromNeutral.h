/**
 *  @file PDSS_IonsFromNeutral.h
 *   Declarations for the class PDSS_IonsFromNeutral (
 *    which handles calculations for a single ion in a fluid, whose properties
 *    are calculated from another neutral molecule.
 *    (see \ref pdssthermo and class \link Cantera::PDSS_IonsFromNeutral PDSS_IonsFromNeutral\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_PDSS_IONSFROMNEUTRAL_H
#define CT_PDSS_IONSFROMNEUTRAL_H

#include "PDSS.h"


namespace Cantera
{
class XML_Node;
class VPStandardStateTP;
class ThermoPhase;


//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * This class is for a single Ideal Gas species.
 *
 * @ingroup pdssthermo
 */
class PDSS_IonsFromNeutral : public PDSS
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
    PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex);

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
    PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex,
                         const std::string& inputFile, const std::string& id = "");


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
    PDSS_IonsFromNeutral(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
                         const XML_Node& phaseRef, bool spInstalled);

    //! Copy Constructor
    /*!
     * @param b Object to be copied
     */
    PDSS_IonsFromNeutral(const PDSS_IonsFromNeutral& b);

    //! Assignment operator
    /*!
     * @param b Object to be copeid
     */
    PDSS_IonsFromNeutral& operator=(const PDSS_IonsFromNeutral& b);

    //! Destructor
    virtual ~PDSS_IonsFromNeutral();

    //! Duplicator
    virtual PDSS* duplMyselfAsPDSS() const;

    //! Initialize or Reinitialize all shallow pointers in the object
    /*!
     *  This command is called to reinitialize all shallow pointers in the
     *  object. It's needed for the duplicator capability.
     *  We need to have an inherited function here to set neutralMoleculePhase_ properly.
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
    virtual void initAllPtrs(VPStandardStateTP* vptp_ptr, VPSSMgr* vpssmgr_ptr,
                             SpeciesThermo* spthermo_ptr);

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

    //! Return the molar gibbs free energy divided by <I>RT</I>
    /*!
     * Returns the species standard state gibbs free energy divided by <I>RT</I> at the
     * current temperature and pressure.
     *
     *  \f[
     *    \frac{\mu^o_k}{RT} = \sum_{m}{ \alpha_{m , k} \frac{\mu^o_{m}}{RT}} + ( 1 - \delta_{k,sp}) 2.0 \ln{2.0}
     *  \f]
     *
     *  <I>m</I> is the neutral molecule species index. \f$ \alpha_{m , k} \f$ is the stoiciometric
     *  coefficient for the neutral molecule,  <I>m</I>, that creates the thermodynamics for the ionic species  <I>k</I>.
     *  A factor  \f$ 2.0 \ln{2.0} \f$ is added to all ions except for the species ionic species, which in this
     *  case is the single anion species, with species index <I>sp</I>.
     *
     * @return Returns the species standard state gibbs free energy divided by <I>RT</I>
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
                           const std::string& inputFile, const std::string& id);

    //!  Initialization of a PDSS object using an xml tree
    /*!
     * This routine is a driver for the initialization of the object.
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
     * @param speciesNode  Reference to the phase Information for the species
     *                     that this standard state refers to
     *
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     *
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& speciesNode,
                          const XML_Node& phaseNode, const std::string& id);

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


    //! Pointer to the Neutral Molecule ThermoPhase object
    /*!
     *  This is a shallow pointer.
     */
    ThermoPhase* neutralMoleculePhase_;

public:

    //! Number of neutral molecule species that make up the stoichiometric vector for
    //! this species, in terms of calculating thermodynamic functions
    size_t numMult_;

    //! Vector of species indices in the neutral molecule ThermoPhase
    std::vector<size_t> idNeutralMoleculeVec;

    //! Stoichiometric coefficient for this species using the Neutral Molecule Species
    //! in the vector idNeutralMoleculeVec
    std::vector<double> factorVec;

    //! Add 2RTln2 to the entropy and Gibbs free energies for this species
    /*!
     *  This is true if this species is not the special species
     */
    bool add2RTln2_;

    //! Vector of length equal to the number of species in the neutral molecule phase
    mutable std::vector<double> tmpNM;

    //! True if this species is the special species
    int specialSpecies_;
};
}

#endif



