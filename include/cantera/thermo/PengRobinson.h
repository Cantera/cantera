//! @file PengRobinson.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PENGROBINSON_H
#define CT_PENGROBINSON_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Peng-Robinson equation of state
 *
 * @ingroup thermoprops
 */
class PengRobinson : public MixtureFugacityTP
{
public:

    //! Construct and initialize a PengRobinson object directly from an
    //! input file
    /*!
     * @param infile    Name of the input file containing the phase YAML data.
     *                  If blank, an empty phase will be created.
     * @param id        ID of the phase in the input file. If empty, the
     *                  first phase definition in the input file will be used.
     */
    explicit PengRobinson(const std::string& infile="",
                              const std::string& id="");

    virtual std::string type() const {
        return "PengRobinson";
    }

    //! @name Molar Thermodynamic properties
    //! @{

    virtual double cp_mole() const;
    virtual double cv_mole() const;

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
     *    P = \frac{RT}{v-b_{mix}} - \frac{\left(\alpha a\right)_{mix}}{v^2 + 2b_{mix}v - b_{mix}^2}
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
     *    \kappa = \left(0.37464 + 1.54226\omega - 0.26992\omega^2\right), \qquad \qquad \text{For } \omega <= 0.491 \\
     *    \kappa = \left(0.379642 + 1.487503\omega - 0.164423\omega^2 +  0.016667\omega^3 \right), \qquad \text{For } \omega > 0.491
     * \f]
     *
     * Coefficients \f$ a_mix, b_mix \f$ and \f$(a \alpha)_{mix}\f$ are calculated as
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
     *   {a \alpha}_{mix} = \sum_i \sum_j X_i X_j {a \alpha}_{i, j} = \sum_i \sum_j X_i X_j \sqrt{a_i a_j} \sqrt{\alpha_i \alpha_j}
     * \f]
     *
     *
     */
    
    virtual double pressure() const;

    // @}

public:

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase.
     * The ideal gas mixture is considered as the standard or reference state here.
     * Since the activity for an ideal gas mixture is simply the mole fraction,
     * for an ideal gas,  \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default is to
     *          assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m^3 / kmol.
     */
    virtual double standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(double* ac) const;

    /// @name  Partial Molar Properties of the Solution
    //@{

    virtual void getChemPotentials(double* mu) const;
    virtual void getPartialMolarEnthalpies(double* hbar) const;
    virtual void getPartialMolarEntropies(double* sbar) const;
    virtual void getPartialMolarIntEnergies(double* ubar) const;
    virtual void getPartialMolarCp(double* cpbar) const;
    virtual void getPartialMolarVolumes(double* vbar) const;

    //! Calculate the temperature dependent interaction parameter \f$\alpha\f$ needed for P-R EoS
    /*!
     *  The temperature dependent parameter in P-R EoS is calculated as
     *       \f$ \alpha = [1 + \kappa(1 - \sqrt{T/T_{crit}})]^2  \f$
     *  kappa is a function calulated based on the acentric factor.
     *  Units: unitless
     * 
     * @param a    species-specific coefficients used in P-R EoS
     * @param b    species-specific coefficients used in P-R EoS
     * @param w    the acentric factor
     */

    virtual void calculateAlpha(const std::string& species, double a, double b, double w);
    
    virtual double speciesCritTemperature(double a, double b) const;

public:
    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void initThermo();

    //! Retrieve a and b coefficients by looking up tabulated critical parameters
    /*!
    *  If pureFluidParameters are not provided for any species in the phase,
    *  consult the critical properties tabulated in data/inputs/critProperties.xml.
    *  If the species is found there, calculate pure fluid parameters a_k and b_k as:
    *  \f[ a_k = 0.4278 R^2 T_c^2 / P_c \f]
    *
    *  and:
    *  \f[ b_k = 0.08664 R T_c/ P_c \f]
    *
    *  @param iName    Name of the species
    */
    virtual vector_fp getCoeff(const std::string& iName);

    //! Set the pure fluid interaction parameters for a species
    /*!
     *
     *  @param species   Name of the species
     *  @param a         "a" parameter in the Peng-Robinson model [Pa-m^6/kmol^2]
     *  @param b         "b" parameter in the Peng-Robinson model [m^3/kmol]
     *  @param w         acentric factor
     */
    void setSpeciesCoeffs(const std::string& species, double a, double b,
                          double w);

    //! Set values for the interaction parameter between two species
    /*!
     *  The "a" parameter for interactions between species *i* and *j* is
     *  assumed by default to be computed as:
     *  \f[ a_{ij} = \sqrt{a_{i, 0} a_{j, 0}} + \sqrt{a_{i, 1} a_{j, 1}} T \f]
     *
     *
     *  @param species_i   Name of one species
     *  @param species_j   Name of the other species
     *  @param a           constant term in the "a" expression [Pa-m^6/kmol^2]
     *  @param alpha       dimensionless function of T_r and \omega
     */
    void setBinaryCoeffs(const std::string& species_i,
                         const std::string& species_j, double a, double alpha);

protected:
    // Special functions inherited from MixtureFugacityTP
    virtual double sresid() const;
    virtual double hresid() const;

    virtual double liquidVolEst(double TKelvin, double& pres) const;
    virtual double densityCalc(double TKelvin, double pressure, int phase, double rhoguess);

    virtual double densSpinodalLiquid() const;
    virtual double densSpinodalGas() const;
    virtual double pressureCalc(double TKelvin, double molarVol) const;
    virtual double dpdVCalc(double TKelvin, double molarVol, double& presCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    double daAlpha_dT() const;
    double d2aAlpha_dT2() const;

public:    

    //! Calculate dpdV and dpdT at the current conditions
    /*!
     *  These are stored internally.
     */
    void calculatePressureDerivatives() const;

    //! Update the a, b and alpha parameters
    /*!
     *  The a and the b parameters depend on the mole fraction and the
     *  parameter alpha depends on the temperature. This function updates 
     *  the internal numbers based on the state of the object.
     */
    virtual void updateMixingExpressions();

    //! Calculate the a, b and the alpha parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     * @param aAlpha (output) Returns the (a*alpha) value.
     */
    void calculateAB(double& aCalc, double& bCalc, double& aAlpha) const;

    void calcCriticalConditions(double& pc, double& tc, double& vc) const;

    //! Prepare variables and call the function to solve the cubic equation of state
    int solveCubic(double T, double pres, double a, double b, double aAlpha,
                   double Vroot[3]) const;
protected:
    //! Value of b in the equation of state
    /*!
     *  m_b_current is a function the mole fractions and species-specific b values.
     */
    double m_b_current;

    //! Value of a and alpha in the equation of state
    /*!
     *  m_aAlpha_current is a function of the temperature and the mole fractions. m_a_current depends only on the mole fractions.
     */
    double m_a_current;
    double m_aAlpha_current;

    // Vectors required to store a_coeff, b_coeff, alpha, kappa and other values for every species. Length = m_kk
    vector_fp m_b_vec_Curr;
    vector_fp m_kappa_vec;
    mutable vector_fp m_dalphadT_vec_Curr;
    mutable vector_fp m_d2alphadT2;
    vector_fp m_alpha_vec_Curr;    

    //Matrices for Binary coefficients a_{i,j} and {a*alpha}_{i.j} are saved in an Array form. Length =  (m_kk, m_kk)
    Array2D m_a_vec_Curr; 
    Array2D m_aAlpha_vec_Curr; 

    int m_NSolns;

    double m_Vroot[3];

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

public:
    //! Omega constants: a0 (= omega_a) and b0 (= omega_b) values used in Peng-Robinson equation of state
    /*!
     *  These values are calculated by solving P-R cubic equation at the critical point.
     */
    static const double omega_a;

    //! Omega constant for b
    static const double omega_b;

    //! Omega constant for the critical molar volume
    static const double omega_vc;
};
}

#endif
