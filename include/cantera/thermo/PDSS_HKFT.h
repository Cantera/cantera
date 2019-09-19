/**
 *  @file PDSS_HKFT.h
 *    Declarations for the class PDSS_HKFT (pressure dependent standard state)
 *    which handles calculations for a single species in a phase using the
 *    HKFT standard state
 *    (see \ref pdssthermo and class \link Cantera::PDSS_HKFT PDSS_HKFT\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_HKFT_H
#define CT_PDSS_HKFT_H

#include "PDSS.h"
#include "WaterProps.h"

namespace Cantera
{
class PDSS_Water;

//! Class for pressure dependent standard states corresponding to
//!  ionic solutes in electrolyte water.
/*!
 * @ingroup pdssthermo
 */
class PDSS_HKFT : public PDSS_Molar
{
public:
    //! Default Constructor
    PDSS_HKFT();

    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_mole() const;

    //! Return the molar enthalpy in units of J kmol-1
    /*!
     * Returns the species standard state enthalpy in J kmol-1 at the
     * current temperature and pressure.
     *
     *  Note this is just an extra routine to check the arithmetic
     *
     * @returns the species standard state enthalpy in J kmol-1
     */
    doublereal enthalpy_mole2() const;

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

    void setParent(VPStandardStateTP* phase, size_t k) {
        m_tp = phase;
        m_spindex = k;
    }

    virtual void initThermo();

     //! Set enthalpy of formation at Pr, Tr [J/kmol]
    void setDeltaH0(double dh0);

    //! Set Gibbs free energy of formation at Pr, Tr [J/kmol]
    void setDeltaG0(double dg0);

     //! Set entropy of formation at Pr, Tr [J/kmol/K]
    void setS0(double s0);

    //! Set "a" coefficients (array of 4 elements). Units of each coefficient
    //! are [J/kmol/Pa, J/kmol, J*K/kmol/Pa, J*K/kmol]
    void set_a(double* a);

    //! Set "c" coefficients (array of 2 elements). Units of each coefficient
    //! are [J/kmol/K, J*K/kmol]
    void set_c(double* c);
    void setOmega(double omega); //!< Set omega [J/kmol]

    void setParametersFromXML(const XML_Node& speciesNode);

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
     * @param kindex    Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the parameters for
     *                  the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(size_t& kindex, int& type, doublereal* const c,
                              doublereal& minTemp, doublereal& maxTemp,
                              doublereal& refPressure) const;
    //@}

private:
    VPStandardStateTP* m_tp; //!< Parent VPStandardStateTP (ThermoPhase) object
    size_t m_spindex; //!< Index of this species within the parent phase

    //! Main routine that actually calculates the Gibbs free energy difference
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

    //! Routine that actually calculates the enthalpy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is an extra routine that was added to check the arithmetic
     */
    doublereal deltaH() const;

    //! Internal formula for the calculation of a_g()
    /*!
     * The output of this is in units of Angstroms
     *
     * @param temp Temperature (K)
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
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    doublereal bg(const doublereal temp, const int ifunc = 0) const;

    //! function g appearing in the formulation
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

    //! Function to look up Element Free Energies
    /*!
     * This function looks up the argument string in the element database and
     * returns the associated 298 K Gibbs Free energy of the element in its
     * stable state.
     *
     * @param  elemName  String. Only the first 3 characters are significant
     * @return value contains the Gibbs free energy for that element
     *
     * @exception CanteraError
     *    If a match is not found, a CanteraError is thrown as well
     */
    doublereal LookupGe(const std::string& elemName);

    //! Translate a Gibbs free energy of formation value to a NIST-based Chemical potential
    /*!
     * Internally, this function is used to translate the input value,
     * m_deltaG_formation_tr_pr, to the internally stored value, m_Mu0_tr_pr.
     */
    void convertDGFormation();

private:
    //! Water standard state calculator
    /*!
     *  derived from the equation of state for water.
     *  This object doesn't own the object. Just a shallow pointer.
     */
    PDSS_Water* m_waterSS;

    //! density of standard-state water. internal temporary variable
    mutable doublereal m_densWaterSS;

    //!  Pointer to the water property calculator
    std::unique_ptr<WaterProps> m_waterProps;

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
