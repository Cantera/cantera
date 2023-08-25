/**
 *  @file PDSS_HKFT.h
 *    Declarations for the class PDSS_HKFT (pressure dependent standard state)
 *    which handles calculations for a single species in a phase using the
 *    HKFT standard state
 *    (see @ref pdssthermo and class @link Cantera::PDSS_HKFT PDSS_HKFT@endlink).
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

    double enthalpy_mole() const override;

    double intEnergy_mole() const override;
    double entropy_mole() const override;
    double gibbs_mole() const override;
    double cp_mole() const override;
    double molarVolume() const override;
    double density() const override;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    double refPressure() const {
        return m_p0;
    }

    double gibbs_RT_ref() const override;
    double enthalpy_RT_ref() const override;
    double entropy_R_ref() const override;
    double cp_R_ref() const override;
    double molarVolume_ref() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    void setState_TP(double temp, double pres) override;

    //! @}
    //! @name Initialization of the Object
    //! @{

    void setParent(VPStandardStateTP* phase, size_t k) override {
        m_tp = phase;
        m_spindex = k;
    }

    void initThermo() override;

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

    void getParameters(AnyMap& eosNode) const override;
    //! @}

private:
    VPStandardStateTP* m_tp; //!< Parent VPStandardStateTP (ThermoPhase) object
    size_t m_spindex; //!< Index of this species within the parent phase

    //! Main routine that actually calculates the Gibbs free energy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is Eqn. 59 in Johnson et al. @cite johnson1992.
     */
    double deltaG() const;

    //! Main routine that actually calculates the entropy difference
    //! between the reference state at Tr, Pr and T,P
    /*!
     *  This is Eqn. 61 in Johnson et al. @cite johnson1992. Actually, there appears to
     *  be an error in the latter. This is a correction.
     */
    double deltaS() const;

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
    double ag(const double temp, const int ifunc = 0) const;

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
    double bg(const double temp, const int ifunc = 0) const;

    //! function g appearing in the formulation
    /*!
     * Function @f$ g @f$ (Eqn. 49) appearing in the Johnson et al. @cite johnson1992
     * formulation.
     *
     * @param temp      Temperature kelvin
     * @param pres      Pressure (pascal)
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    double g(const double temp, const double pres, const int ifunc = 0) const;

    //! Difference function f appearing in the formulation
    /*!
     * Function @f$ f @f$ (Eqn. 52) appearing in the Johnson et al. @cite johnson1992
     * formulation of @f$ \omega_j @f$ (Eqn. 46).
     *
     * @param temp      Temperature kelvin
     * @param pres      Pressure (pascal)
     * @param ifunc     parameters specifying the desired information
     *                 - 0 function value
     *                 - 1 derivative wrt temperature
     *                 - 2 2nd derivative wrt temperature
     *                 - 3 derivative wrt pressure
     */
    double f(const double temp, const double pres, const int ifunc = 0) const;

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
    double gstar(const double temp, const double pres, const int ifunc = 0) const;

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
    double LookupGe(const string& elemName);

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
    PDSS_Water* m_waterSS = nullptr;

    UnitSystem m_units;

    //! density of standard-state water. internal temporary variable
    mutable double m_densWaterSS = -1.0;

    //!  Pointer to the water property calculator
    unique_ptr<WaterProps> m_waterProps;

    //! Input value of deltaG of Formation at Tr and Pr    (cal gmol-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     *
     *  This is the delta G for the formation reaction of the
     *  ion from elements in their stable state at Tr, Pr.
     */
    double m_deltaG_formation_tr_pr = NAN;

    //!  Input value of deltaH of Formation at Tr and Pr    (cal gmol-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     *
     *  This is the delta H for the formation reaction of the
     *  ion from elements in their stable state at Tr, Pr.
     */
    double m_deltaH_formation_tr_pr = NAN;

    //! Value of the Absolute Gibbs Free Energy NIST scale at T_r and P_r
    /*!
     *  This is the NIST scale value of Gibbs free energy at T_r = 298.15
     *  and P_r = 1 atm.
     *
     *  J kmol-1
     */
    double m_Mu0_tr_pr = 0.0;

    //! Input value of S_j at Tr and Pr    (cal gmol-1 K-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     */
    double m_Entrop_tr_pr = NAN;

    //! Input a1 coefficient (cal gmol-1 bar-1)
    double m_a1 = 0.0;

    //!  Input a2 coefficient (cal gmol-1)
    double m_a2 = 0.0;

    //!  Input a3 coefficient (cal K gmol-1 bar-1)
    double m_a3 = 0.0;

    //!  Input a4 coefficient (cal K gmol-1)
    double m_a4 = 0.0;

    //!  Input c1 coefficient (cal gmol-1 K-1)
    double m_c1 = 0.0;

    //!  Input c2 coefficient (cal K gmol-1)
    double m_c2 = 0.0;

    //! Input  omega_pr_tr coefficient(cal gmol-1)
    double m_omega_pr_tr = 0.0;

    //! y = dZdT = 1/(esp*esp) desp/dT at 298.15 and 1 bar
    double m_Y_pr_tr = 0.0;

    //! Z = -1 / relEpsilon at 298.15 and 1 bar
    double m_Z_pr_tr = 0.0;

    //! Reference pressure is 1 atm in units of bar= 1.0132
    double m_presR_bar = OneAtm * 1.0E-5;

    //! small value that is not quite zero
    double m_domega_jdT_prtr = 0.0;

    //! Charge of the ion
    double m_charge_j = 0.0;

    //!  Static variable determining error exiting
    /*!
     *   If true, then will error exit if there is an inconsistency in DG0, DH0, and DS0.
     *   If not, then will rewrite DH0 to be consistent with the other two.
     */
    static int s_InputInconsistencyErrorExit;
};

}

#endif
