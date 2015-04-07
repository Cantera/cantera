/**
 *  @file PDSS_HKFT.h
 *    Declarations for the class PDSS_HKFT (pressure dependent standard state)
 *    which handles calculations for a single species in a phase using the
 *    HKFT standard state
 *    (see \ref pdssthermo and class \link Cantera::PDSS_HKFT PDSS_HKFT\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_PDSS_HKFT_H
#define CT_PDSS_HKFT_H

class WaterPropsIAPWS;
#include "PDSS.h"

namespace Cantera
{
class XML_Node;
class VPStandardState;
class PDSS_Water;
class WaterProps;

//! Class for pressure dependent standard states corresponding to
//!  ionic solutes in electrolyte water.
/*!
 * @ingroup pdssthermo
 */
class PDSS_HKFT : public PDSS_Molar
{
public:

    

    //! @name  Constructors
    //! @{

    //! Constructor that initializes the object by examining the XML entries
    //! from the ThermoPhase object
    /*!
     *  This function calls the constructPDSS member function.
     *
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS_HKFT(VPStandardStateTP* tp, size_t spindex);

    //! Copy Constructor
    /*!
     * @param b object to be copied
     */
    PDSS_HKFT(const PDSS_HKFT& b);

    //! Assignment operator
    /*!
     * @param b Object to be copied
     */
    PDSS_HKFT& operator=(const PDSS_HKFT& b);

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSFile member function.
     *
     *  @param vptp_ptr  Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param inputFile String name of the input file
     *  @param id        String name of the phase in the input file. The default
     *                   is the empty string, in which case the first phase in the
     *                   file is used.
     */
    PDSS_HKFT(VPStandardStateTP* vptp_ptr, size_t spindex,
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
    PDSS_HKFT(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
              const XML_Node& phaseRef, bool spInstalled);

    //! Destructor for the phase
    virtual ~PDSS_HKFT();

    virtual PDSS* duplMyselfAsPDSS() const;

    //! @}
    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_mole() const;

#ifdef DEBUG_MODE
    //! Return the molar enthalpy in units of J kmol-1
    /*!
     * Returns the species standard state enthalpy in J kmol-1 at the
     * current temperature and pressure.
     *
     *  Note this is just an extra routine to check the arithmetic
     *
     * @return returns the species standard state enthalpy in  J kmol-1
     */
    doublereal enthalpy_mole2() const;
#endif

    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal molarVolume() const;
    virtual doublereal density() const;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    doublereal refPressure() const {
        return m_p0;
    }

    virtual doublereal gibbs_RT_ref() const;
    virtual doublereal enthalpy_RT_ref() const;
    virtual doublereal entropy_R_ref() const;
    virtual doublereal cp_R_ref() const;
    virtual doublereal molarVolume_ref() const;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    virtual void setState_TP(doublereal temp, doublereal pres);

    //! @}
    //! @name Initialization of the Object
    //! @{

    virtual void initThermo();

    //! Initialization of a PDSS object using an input XML file.
    /*!
     * This routine is a precursor to constructPDSSXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param vptp_ptr    Pointer to the Variable pressure %ThermoPhase object
     *                    This object must have already been malloced.
     * @param spindex     Species index within the phase
     * @param inputFile   XML file containing the description of the
     *                    phase
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
     *     - initThermo()                 (cascade)
     *     - getStuff from species Part of XML file
     *     - initThermoXML(phaseNode)      (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure %ThermoPhase object
     *                   This object must have already been malloced.
     * @param spindex    Species index within the phase
     * @param speciesNode XML Node containing the species information
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     * @param spInstalled  Boolean indicating whether the species is
     *                     already installed.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& speciesNode,
                          const XML_Node& phaseNode, bool spInstalled);

    virtual void initAllPtrs(VPStandardStateTP* vptp_ptr, VPSSMgr* vpssmgr_ptr,
                             SpeciesThermo* spthermo_ptr);

    //! This utility function reports back the type of parameterization and
    //! all of the parameters for the species, index.
    /*!
     * The following parameters are reported
     *
     * -   c[0] = m_deltaG_formation_tr_pr;
     * -   c[1] = m_deltaH_formation_tr_pr;
     * -   c[2] = m_Mu0_tr_pr;
     * -   c[3] = m_Entrop_tr_pr;
     * -   c[4] =  m_a1;
     * -   c[5] =  m_a2;
     * -   c[6] =  m_a3;
     * -   c[7] =  m_a4;
     * -   c[8] =  m_c1;
     * -   c[9] =  m_c2;
     * -   c[10] = m_omega_pr_tr;
     * .
     *
     * @param kindex     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(size_t& kindex, int& type, doublereal* const c,
                              doublereal& minTemp, doublereal& maxTemp,
                              doublereal& refPressure) const;

    //@}

private:
    //! Main routine that actually calculates the gibbs free energy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is eEqn. 59 in Johnson et al. (1992).
     */
    doublereal deltaG() const;

    //! Main routine that actually calculates the entropy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is Eqn. 61 in Johnson et al. (1992). Actually, there appears to
     *  be an error in the latter. This is a correction.
     */
    doublereal deltaS() const;

#ifdef DEBUG_MODE
    //! Routine that actually calculates the enthalpy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is an extra routine that was added to check the arithmetic
     */
    doublereal deltaH() const;
#endif

    //! Internal formula for the calculation of a_g()
    /*!
     * The output of this is in units of Angstroms
     *
     * @param temp Temperature (K)
     *
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal ag(const doublereal temp, const int ifunc = 0) const;

    //! Internal formula for the calculation of b_g()
    /*!
     * the output of this is unitless
     *
     * @param temp Temperature (K)
     *
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal bg(const doublereal temp, const int ifunc = 0) const;

    //!  function g appearing in the formulation
    /*!
     * Function g appearing in the Johnson et al formulation
     *
     * @param temp      Temperature kelvin
     * @param pres      Pressure (pascal)
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal g(const doublereal temp, const doublereal pres, const int ifunc = 0) const;

    //! Difference function f appearing in the formulation
    /*!
     * Function f appearing in the Johnson et al formulation of omega_j
     *   Eqn. 33 ref
     *
     * @param temp      Temperature kelvin
     * @param pres      Pressure (pascal)
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal f(const doublereal temp, const doublereal pres, const int ifunc = 0) const;

    //! Evaluate the Gstar value appearing in the HKFT formulation
    /*!
     * @param temp      Temperature kelvin
     * @param pres      Pressure (pascal)
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal gstar(const doublereal temp, const doublereal pres,
                     const int ifunc = 0) const;

    //!  Function to look up Element Free Energies
    /*!
     * This function looks up the argument string in the element database and
     * returns the associated 298 K Gibbs Free energy of the element in its
     * stable state.
     *
     * @param  elemName  String. Only the first 3 characters are significant
     *
     * @return value contains the Gibbs free energy for that element
     *
     * @exception CanteraError
     *    If a match is not found, a CanteraError is thrown as well
     */
    doublereal LookupGe(const std::string& elemName);

    //! Translate a Gibbs free energy of formation value to a NIST-based Chemical potential
    /*!
     * Internally, this function is used to translate the input value,
     * m_deltaG_formation_tr_pr, to the internally stored value,  m_Mu0_tr_pr.
     */
    void convertDGFormation();

private:
    //!  Water standard state calculator
    /*!
     *  derived from the equation of state for water.
     *  This object doesn't own the object. Just a shallow pointer.
     */
    PDSS_Water* m_waterSS;

    //! density of standard-state water
    /*!
     * internal temporary variable
     */
    mutable doublereal m_densWaterSS;

    //!  Pointer to the water property calculator
    WaterProps* m_waterProps;

    //! Born coefficient for the current ion or species
    doublereal m_born_coeff_j;

    //! Electrostatic radii
    doublereal m_r_e_j;

    //! Input value of deltaG of Formation at Tr and Pr    (cal gmol-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     *
     *  This is the delta G for the formation reaction of the
     *  ion from elements in their stable state at Tr, Pr.
     */
    doublereal m_deltaG_formation_tr_pr;

    //!  Input value of deltaH of Formation at Tr and Pr    (cal gmol-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     *
     *  This is the delta H for the formation reaction of the
     *  ion from elements in their stable state at Tr, Pr.
     */
    doublereal m_deltaH_formation_tr_pr;

    //! Value of the Absolute Gibbs Free Energy NIST scale at T_r and P_r
    /*!
     *  This is the NIST scale value of Gibbs free energy at T_r = 298.15
     *  and P_r = 1 atm.
     *
     *  J kmol-1
     */
    doublereal m_Mu0_tr_pr;

    //! Input value of S_j at Tr and Pr    (cal gmol-1 K-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     */
    doublereal m_Entrop_tr_pr;

    //! Input a1 coefficient (cal gmol-1 bar-1)
    doublereal m_a1;

    //!  Input a2 coefficient (cal gmol-1)
    doublereal m_a2;

    //!  Input a3 coefficient (cal K gmol-1 bar-1)
    doublereal m_a3;

    //!  Input a4 coefficient (cal K gmol-1)
    doublereal m_a4;

    //!  Input c1 coefficient (cal gmol-1 K-1)
    doublereal m_c1;

    //!  Input c2 coefficient (cal K gmol-1)
    doublereal m_c2;

    //! Input  omega_pr_tr coefficient(cal gmol-1)
    doublereal m_omega_pr_tr;

    //! y = dZdT = 1/(esp*esp) desp/dT at 298.15 and 1 bar
    doublereal m_Y_pr_tr;

    //! Z = -1 / relEpsilon at 298.15 and 1 bar
    doublereal m_Z_pr_tr;

    //! Reference pressure is 1 atm in units of bar= 1.0132
    doublereal m_presR_bar;

    //! small value that is not quite zero
    doublereal m_domega_jdT_prtr;

    //! Charge of the ion
    doublereal m_charge_j;

    //!  Static variable determining error exiting 
    /*!
     *   If true, then will error exit if there is an inconsistency in DG0, DH0, and DS0.
     *   If not, then will rewrite DH0 to be consistent with the other two.
     */
    static int s_InputInconsistencyErrorExit;
};

}

#endif
