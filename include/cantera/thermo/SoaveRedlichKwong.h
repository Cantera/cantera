//! @file SoaveRedlichKwong.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOAVEREDLICHKWONG_H
#define CT_SOAVEREDLICHKWONG_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Soave-Redlich-Kwong equation of state
 *
 * @ingroup thermoprops
 */
class SoaveRedlichKwong : public MixtureFugacityTP
{
public:
    explicit SoaveRedlichKwong(const std::string& infile="",
                               const std::string& id="");

    virtual std::string type() const {
        return "Soave-Redlich-Kwong";
    }

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     * Since the mass density, temperature, and mass fractions are stored,
     * this method uses these values to implement the
     * mechanical equation of state \f$ P(T, \rho, Y_1, \dots, Y_K) \f$.
     *
     * \f[
     *    P = \frac{RT}{v-b_{mix}}
     *        - \frac{\left(a\alpha\right)_{mix}}{v\left(v + b_{mix}\right)}
     * \f]
     *
     * where:
     *
     * \f[
     *    \alpha = \left[ 1 + \kappa \left(1-T_r^{0.5}\right)\right]^2
     * \f]
     *
     *  and
     *
     * \f[
     *    \kappa = \left(0.48508 + 1.55171\omega - 0.15613\omega^2\right)
     * \f]
     *
     * Coefficients \f$ a_{mix}, b_{mix} \f$ and \f$(a \alpha)_{mix}\f$ are calculated as
     *
     * \f[
     *    a_{mix} = \sum_i \sum_j X_i X_j a_{i, j} = \sum_i \sum_j X_i X_j \sqrt{a_i a_j}
     * \f]
     *
     * \f[
     *    b_{mix} = \sum_i X_i b_i
     * \f]
     *
     * \f[
     *   {a \alpha}_{mix} = \sum_i \sum_j X_i X_j {a \alpha}_{i, j}
     *       = \sum_i \sum_j X_i X_j \sqrt{a_i a_j} \sqrt{\alpha_i \alpha_j}
     * \f]
     */
    virtual double pressure() const;


    //! Calculate species-specific critical temperature
    /*!
     *  The temperature dependent parameter in SRK EoS is calculated as
     *       \f[ T_{crit} = (0.0778 a)/(0.4572 b R) \f] TODO
     *  Units: Kelvin
     *
     * @param a    species-specific coefficients used in P-R EoS
     * @param b    species-specific coefficients used in P-R EoS
     */
    double speciesCritTemperature(double a, double b) const;

    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void initThermo();
    virtual void getSpeciesParameters(const std::string& name,
                                      AnyMap& speciesNode) const;

protected:
    // Special functions inherited from MixtureFugacityTP

    virtual double liquidVolEst(double T, double& pres) const;
    virtual double densityCalc(double T, double pressure, int phase, double rhoguess);

    virtual double densSpinodalLiquid() const;
    virtual double densSpinodalGas() const;
    virtual double dpdVCalc(double T, double molarVol, double& presCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    //! Calculate temperature derivative \f$d(a \alpha)/dT\f$
    /*!
     *  These are stored internally.
     */
    double daAlpha_dT() const;

    //! Calculate second derivative \f$d^2(a \alpha)/dT^2\f$
    /*!
     *  These are stored internally.
     */
    double d2aAlpha_dT2() const;

public:
    //! Returns the isothermal compressibility. Units: 1/Pa.
    double isothermalCompressibility() const;

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    double thermalExpansionCoeff() const;

    //! Calculate \f$dp/dV\f$ and \f$dp/dT\f$ at the current conditions
    /*!
     *  These are stored internally.
     */
    void calculatePressureDerivatives() const;

    //! Update the \f$a\f$, \f$b\f$, and \f$\alpha\f$ parameters
    /*!
     *  The \f$a\f$ and the \f$b\f$ parameters depend on the mole fraction and the
     *  parameter \f$\alpha\f$ depends on the temperature. This function updates
     *  the internal numbers based on the state of the object.
     */
    virtual void updateMixingExpressions();

    //! Calculate the \f$a\f$, \f$b\f$, and \f$\a alpha\f$ parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     * @param aAlpha (output) Returns the (a*alpha) value.
     */
    void calculateAB(double& aCalc, double& bCalc, double& aAlphaCalc) const

    void calcCriticalConditions(double& pc, double& tc, double& vc) const;

    //! Prepare variables and call the function to solve the cubic equation of state
    int solveCubic(double T, double pres, double a, double b, double aAlpha,
                   double Vroot[3]) const;

protected:
    //! Value of \f$b\f$ in the equation of state
    /*!
     *  `m_b` is a function of the mole fractions and species-specific b values.
     */
    double m_b;

    //! Value of \f$a\f$ in the equation of state
    /*!
     *  `m_a` depends only on the mole fractions.
     */
    double m_a;

    //! Value of \f$\alpha\f$ in the equation of state
    /*!
     *  `m_aAlpha_mix` is a function of the temperature and the mole fractions.
     */
    double m_aAlpha_mix;

    // Vectors required to store species-specific a_coeff, b_coeff, alpha, kappa
    // and other derivatives. Length = m_kk.
    vector_fp m_b_coeffs;
    vector_fp m_kappa;
    vector_fp m_acentric; //!< acentric factor for each species, length #m_kk
    mutable vector_fp m_dalphadT;
    mutable vector_fp m_d2alphadT2;
    vector_fp m_alpha;

    // Matrices for Binary coefficients a_{i,j} and {a*alpha}_{i.j} are saved in an
    // array form. Size = (m_kk, m_kk).
    Array2D m_a_coeffs;
    Array2D m_aAlpha_binary;

    //! Explicitly-specified binary interaction parameters, to enable serialization
    std::map<std::string, std::map<std::string, double>> m_binaryParameters;

    int m_NSolns;

    double m_Vroot[3];

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_pp;

    // Partial molar volumes of the species
    mutable vector_fp m_partialMolarVolumes;

    //! The derivative of the pressure with respect to the volume
    /*!
     * Calculated at the current conditions. temperature and mole number kept
     * constant
     */
    mutable double m_dpdV;

    //! The derivative of the pressure with respect to the temperature
    /*!
     *  Calculated at the current conditions. Total volume and mole number kept
     *  constant
     */
    mutable double m_dpdT;

    //! Vector of derivatives of pressure with respect to mole number
    /*!
     *  Calculated at the current conditions. Total volume, temperature and
     *  other mole number kept constant
     */
    mutable vector_fp m_dpdni;

    enum class CoeffSource { EoS, CritProps, Database };
    //! For each species, specifies the source of the a, b, and omega coefficients
    std::vector<CoeffSource> m_coeffSource;

private:
    //! Omega constant: omega_a used in Soave-Redlich-Kwong equation of state
    /*!
     *  This value is calculated by solving SRK cubic equation at the critical point.
     */
    static const double omega_a;

    //! Omega constant: omega_b used in Soave-Redlich-Kwong equation of state
    /*!
     *  This value is calculated by solving SRK cubic equation at the critical point.
     */
    static const double omega_b;

    //! Omega constant for the critical molar volume
    static const double omega_vc;
};
}

#endif
