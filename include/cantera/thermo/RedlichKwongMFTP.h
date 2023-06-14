//! @file RedlichKwongMFTP.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REDLICHKWONGMFTP_H
#define CT_REDLICHKWONGMFTP_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Redlich-Kwong equation of state
 *
 * @ingroup thermoprops
 */
class RedlichKwongMFTP : public MixtureFugacityTP
{
public:
    //! Construct a RedlichKwongMFTP object from an input file
    /*!
     * @param infile Name of the input file containing the phase definition.
     *               If blank, an empty phase will be created.
     * @param id     name (ID) of the phase in the input file. If empty, the
     *               first phase definition in the input file will be used.
     */
    explicit RedlichKwongMFTP(const std::string& infile="",
                              const std::string& id="");

    virtual std::string type() const {
        return "Redlich-Kwong";
    }

    //! @name Molar Thermodynamic properties
    //! @{
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;
    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  Since the mass density, temperature, and mass fractions are stored,
     *  this method uses these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots, Y_K) \f$.
     *
     * \f[
     *    P = \frac{RT}{v-b_{mix}} - \frac{a_{mix}}{T^{0.5} v \left( v + b_{mix} \right) }
     * \f]
     */
    virtual doublereal pressure() const;

    //! @}

public:

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase. Since the activity
     * for an ideal gas mixture is simply the mole fraction, for an ideal gas
     * \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default is to
     *          assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    //! @name  Partial Molar Properties of the Solution
    //! @{

    //! Get the array of non-dimensional species chemical potentials.
    //! These are partial molar Gibbs free energies.
    /*!
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function, here, calling getChemPotentials() and
     * then dividing by RT. No need for child classes to handle.
     *
     * @param mu    Output vector of non-dimensional species chemical potentials
     *              Length: m_kk.
     * @deprecated To be removed after Cantera 3.0. Use getChemPotentials() instead.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;
    virtual void getPartialMolarCp(double* cpbar) const {
        throw NotImplementedError("RedlichKwongMFTP::getPartialMolarCp");
    }
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    //! @}

public:
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

    //! Set the pure fluid interaction parameters for a species
    /*!
     *  The "a" parameter for species *i* in the Redlich-Kwong model is assumed
     *  to be a linear function of temperature:
     *  \f[ a = a_0 + a_1 T \f]
     *
     *  @param species   Name of the species
     *  @param a0        constant term in the expression for the "a" parameter
     *      of the specified species [Pa-m^6/kmol^2]
     *  @param a1        temperature-proportional term in the expression for the
     *      "a" parameter of the specified species [Pa-m^6/kmol^2/K]
     *  @param b         "b" parameter in the Redlich-Kwong model [m^3/kmol]
     */
    void setSpeciesCoeffs(const std::string& species, double a0, double a1,
                              double b);

    //! Set values for the interaction parameter between two species
    /*!
     *  The "a" parameter for interactions between species *i* and *j* is
     *  assumed by default to be computed as:
     *  \f[ a_{ij} = \sqrt(a_{i,0} a_{j,0}) + \sqrt(a_{i,1} a_{j,1}) T \f]
     *
     *  This function overrides the defaults with the specified parameters:
     *  \f[ a_{ij} = a_{ij,0} + a_{ij,1} T \f]
     *
     *  @param species_i   Name of one species
     *  @param species_j   Name of the other species
     *  @param a0          constant term in the "a" expression [Pa-m^6/kmol^2]
     *  @param a1          temperature-proportional term in the "a" expression
     *      [Pa-m^6/kmol^2/K]
     */
    void setBinaryCoeffs(const std::string& species_i,
                         const std::string& species_j, double a0, double a1);
    //! @}

protected:
    // Special functions inherited from MixtureFugacityTP
    virtual doublereal sresid() const;
    virtual doublereal hresid() const;

public:
    virtual doublereal liquidVolEst(doublereal TKelvin, doublereal& pres) const;
    virtual doublereal densityCalc(doublereal T, doublereal pressure, int phase, doublereal rhoguess);

    virtual doublereal densSpinodalLiquid() const;
    virtual doublereal densSpinodalGas() const;
    virtual doublereal dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const;

    virtual double isothermalCompressibility() const;
    virtual double thermalExpansionCoeff() const;
    virtual double soundSpeed() const;

    //! Calculate dpdV and dpdT at the current conditions
    /*!
     *  These are stored internally.
     */
    void pressureDerivatives() const;

    //! Update the a and b parameters
    /*!
     *  The a and the b parameters depend on the mole fraction and the
     *  temperature. This function updates the internal numbers based on the
     *  state of the object.
     */
    virtual void updateMixingExpressions();

    //! Calculate the a and the b parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param temp  Temperature (TKelvin)
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     */
    void calculateAB(doublereal temp, doublereal& aCalc, doublereal& bCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    doublereal da_dt() const;

    void calcCriticalConditions(doublereal& pc, doublereal& tc, doublereal& vc) const;

    //! Prepare variables and call the function to solve the cubic equation of state
    int solveCubic(double T, double pres, double a, double b, double Vroot[3]) const;

protected:
    //! Form of the temperature parameterization
    /*!
     *  - 0 = There is no temperature parameterization of a or b
     *  - 1 = The a_ij parameter is a linear function of the temperature
     */
    int m_formTempParam = 0;

    //! Value of b in the equation of state
    /*!
     *  m_b is a function of the temperature and the mole fraction.
     */
    double m_b_current = 0.0;

    //! Value of a in the equation of state
    /*!
     *  a_b is a function of the temperature and the mole fraction.
     */
    double m_a_current = 0.0;

    vector_fp a_vec_Curr_;
    vector_fp b_vec_Curr_;

    Array2D a_coeff_vec;

    //! Explicitly-specified binary interaction parameters
    std::map<std::string, std::map<std::string, std::pair<double, double>>> m_binaryParameters;

    enum class CoeffSource { EoS, CritProps, Database };
    //! For each species, specifies the source of the a and b coefficients
    std::vector<CoeffSource> m_coeffSource;

    int NSolns_ = 0;

    double Vroot_[3] = {0.0, 0.0, 0.0};

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_pp;

    // Partial molar volumes of the species
    mutable vector_fp m_partialMolarVolumes;

    //! The derivative of the pressure wrt the volume
    /*!
     * Calculated at the current conditions. temperature and mole number kept
     * constant
     */
    mutable double dpdV_ = 0.0;

    //! The derivative of the pressure wrt the temperature
    /*!
     *  Calculated at the current conditions. Total volume and mole number kept
     *  constant
     */
    mutable double dpdT_ = 0.0;

    //! Vector of derivatives of pressure wrt mole number
    /*!
     *  Calculated at the current conditions. Total volume, temperature and
     *  other mole number kept constant
     */
    mutable vector_fp dpdni_;

private:
    //! Omega constant for a -> value of a in terms of critical properties
    /*!
     *  this was calculated from a small nonlinear solve
     */
    static const doublereal omega_a;

    //! Omega constant for b
    static const doublereal omega_b;

    //! Omega constant for the critical molar volume
    static const doublereal omega_vc;
};
}

#endif
