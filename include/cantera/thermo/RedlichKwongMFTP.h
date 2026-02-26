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
 * # Implementation of a multi-species Redlich-Kwong equation of state
 *
 * The equation of state is
 * @f[
 * P(T,v,\mathbf{X}) = \frac{RT}{v-b} - \frac{a}{\sqrt{T} v (v+b)}
 * @f]
 * with mole-fraction mixing rules
 * @f[
 * a = \sum_i \sum_j X_i X_j a_{ij}(T),\qquad
 * b = \sum_i X_i b_i
 * @f]
 * and (in the present model) linear temperature dependence
 * @f[
 * a_{ij}(T)=a_{ij,0}+a_{ij,1}T.
 * @f]
 *
 * ## Notation used in property derivations {#redlich-Kwong-notation}
 *
 * Several equation of state derivatives are used in calculating the thermodynamic
 * properties of the mixture. These include:
 * @f[
 * P_T = \frac{R}{v-b}
 *   - \frac{1}{\sqrt{T} v(v+b)}
 *   \left(\frac{\partial a}{\partial T} - \frac{a}{2T}\right)
 * @f]
 * @f[
 * P_v = -\frac{RT}{(v-b)^2} + \frac{a(2v+b)}{\sqrt{T} v^2(v+b)^2}
 * @f]
 * @f[
 * v_T = -\frac{P_T}{P_v},\qquad
 * v_{TT} = -\frac{P_{TT} + 2 P_{Tv} v_T + P_{vv} v_T^2}{P_v}.
 * @f]
 *
 * where subscript notation such as @f$P_T@f$, @f$P_{Tv}@f$, and @f$v_{TT}@f$ denote
 * partial derivatives with composition and the remaining state variable held fixed.
 *
 * The following shorthand is used in some property calculations:
 * @f[
 * A_k = \sum_i X_i a_{ki},\qquad
 * A'_k = \sum_i X_i \frac{d a_{ki}}{dT},\qquad
 * S_k = 2T A'_k - 3A_k
 * @f]
 * @f[
 * F = T\frac{da}{dT} - \frac{3}{2}a,\qquad
 * L = \ln \left(\frac{v+b}{v}\right),\qquad
 * C = -\frac{L}{b} + \frac{1}{v+b}.
 * @f]
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
    explicit RedlichKwongMFTP(const string& infile="", const string& id="");

    string type() const override {
        return "Redlich-Kwong";
    }

    //! @name Molar Thermodynamic properties
    //! @{
    double cp_mole() const override;
    double cv_mole() const override;
    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  Since the mass density, temperature, and mass fractions are stored,
     *  this method uses these values to implement the
     *  mechanical equation of state @f$ P(T, \rho, Y_1, \dots, Y_K) @f$.
     *
     * @f[
     *    P = \frac{RT}{v-b_{mix}} - \frac{a_{mix}}{T^{0.5} v \left( v + b_{mix} \right) }
     * @f]
     */
    double pressure() const override;

    //! @}

public:

    //! Returns the standard concentration @f$ C^0_k @f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase. Since the activity
     * for an ideal gas mixture is simply the mole fraction, for an ideal gas
     * @f$ C^0_k = P/\hat R T @f$.
     *
     * @param k Optional parameter indicating the species. The default is to
     *          assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    double standardConcentration(size_t k=0) const override;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    void getActivityCoefficients(span<double> ac) const override;

    //! @name  Partial Molar Properties of the Solution
    //! @{

    void getChemPotentials(span<double> mu) const override;

    //! Return partial molar enthalpies @f$ \bar{h}_k @f$ (J/kmol).
    /*!
     * The partial molar enthalpies for the Redlich-Kwong equation of state are:
     * @f[
     * \bar{h}_k = h_k^\t{ref} + h^E_k
     *   - \left(v - T \left.\frac{\partial v}{\partial T}\right|_{P,\mathbf{X}}\right)
     * \left(\frac{\partial P}{\partial n_k}\right)_{T,V,n_j}
     * @f]
     * where the implementation starts from the ideal-gas reference values
     * @f$ h_k^\t{ref}(T) @f$ and adds a constant volume departure term
     * @f$ h^E_k @f$ and a final term to adjust to constant pressure.
     *
     * Here,
     * @f[
     * \left(\frac{\partial P}{\partial n_k}\right)_{T,V,n_j} =
     *     \frac{RT}{v-b} + \frac{RT b_k}{(v-b)^2}
     *     - \frac{2 A_k}{v(v+b) \sqrt{T}} + \frac{a b_k}{v(v+b)^2 \sqrt{T}},
     * @f]
     * and the intermediate departure term is
     * @f[
     * h^E_k = v\left(\frac{\partial P}{\partial n_k}\right)_{T,V,n_j}
     *     - RT - \frac{b_k}{b^2\sqrt{T}} L F + \frac{1}{b\sqrt{T}} L S_k
     *     + \frac{b_k}{(v+b)b\sqrt{T}} F
     * @f]
     * where
     * @f[
     * L = \ln \left(\frac{v+b}{v}\right),\quad
     * F = T\left.\frac{\partial a}{\partial T}\right|_{\mathbf{X}} - \frac{3}{2}a,\quad
     * S_k = 2T A'_k - 3A_k.
     * @f]
     */
    void getPartialMolarEnthalpies(span<double> hbar) const override;
    void getPartialMolarEntropies(span<double> sbar) const override;
    void getPartialMolarIntEnergies(span<double> ubar) const override;

    //! Get the species molar internal energies associated with the derivatives
    //! of total internal energy at constant-volume [J/kmol].
    /*!
     * This method computes
     * @f[
     * \tilde{u}_k = \frac{\partial U}{\partial n_k}\Bigg|_{T,V,n_{j\ne k}}
     *     = u_k^\t{ref} + \frac{1}{b\sqrt{T}} \left(L S_k + b_k F C\right)
     * @f]
     * Where @f$ L, S_k, F, @f$ and @f$ C @f$ are defined in the
     * [notation section](#redlich-Kwong-notation) and @f$ u_k^\t{ref} @f$ is the
     * ideal-gas reference internal energy.
     *
     * For non-ideal phases like Redlich-Kwong, these are distinct from the
     * partial molar internal energies
     * @f$ \bar{u}_k = \left(\partial U/\partial n_k\right)_{T,P,n_{j\ne k}} @f$
     * calculated by the getPartialMolarIntEnergies() method.
     */
    void getPartialMolarIntEnergies_TV(span<double> ubar) const override;

    //! Get the partial molar heat capacities at constant pressure [J/kmol/K].
    /*!
     * The partial molar heat capacities at constant pressure are defined as
     * @f[
     * \bar{c}_{p,k} = \left(\partial \bar{h}_k / \partial T\right)_{P,\mathbf{X}}
     * @f]
     * Expanding this expression via the chain rule at constant pressure and composition
     * gives:
     * @f[
     * \bar{c}_{p,k} = c_{p,k}^\t{ref}
     * + \frac{d h^{E,v}_k}{dT}\Bigg|_{P,\mathbf{X}}
     * + T v_{TT} \Pi_k
     * - \left(v - T v_T\right)\frac{d\Pi_k}{dT}\Bigg|_{P,\mathbf{X}},
     * @f]
     * where
     * @f[
     * \Pi_k \equiv \left(\frac{\partial P}{\partial n_k}\right)_{T,V,n_j}
     * @f]
     * and @f$ v_T @f$ and @f$ v_{TT} @f$ are the partial derivatives of the molar
     * volume as given in the [notation section](#redlich-Kwong-notation).
     *
     * The derivative of the intermediate departure term can be expanded as
     * @f{eqnarray*}{
     * \frac{d h^{E,v}_k}{dT}\Bigg|_{P,\mathbf{X}}
     * &=& v_T \Pi_k + v \frac{d\Pi_k}{dT}\Bigg|_{P,\mathbf{X}} - R
     * - \frac{b_k}{b^2}\frac{d}{dT}\left(\frac{L F}{\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}} \\
     * &+& \frac{1}{b}\frac{d}{dT}\left(\frac{L S_k}{\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}}
     * + \frac{b_k}{b}\frac{d}{dT}\left(\frac{F}{(v+b)\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}}.
     * @f}
     *
     * The temperature derivatives appearing here are path derivatives at constant
     * pressure and composition, and include both explicit temperature dependence and
     * implicit dependence through @f$ v(T)\vert_{P,\mathbf{X}} @f$. These terms are
     * computed as:
     * @f[
     * \frac{d}{dT}\left(\frac{L F}{\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}}
     * = \frac{1}{\sqrt{T}}\left(L_T F + L F_T - \frac{L F}{2T}\right),
     * \qquad
     *
     * @f]
     * @f[
     * \frac{d}{dT}\left(\frac{L S_k}{\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}}
     * = \frac{1}{\sqrt{T}}\left(L_T S_k + L S_{k,T} - \frac{L S_k}{2T}\right),
     * @f]
     * @f[
     * \frac{d}{dT}\left(\frac{F}{(v+b)\sqrt{T}}\right)\Bigg|_{P,\mathbf{X}}
     * = \frac{F_T}{(v+b)\sqrt{T}}
     * - \frac{F v_T}{(v+b)^2\sqrt{T}}
     * - \frac{F}{2T(v+b)\sqrt{T}}
     * @f]
     * where
     * @f[
     * L_T \equiv \frac{dL}{dT}\Bigg|_{P,\mathbf{X}} = -\frac{b v_T}{v(v+b)},
     * @f]
     * @f[
     * F_T \equiv \frac{dF}{dT}\Bigg|_{P,\mathbf{X}}
     *     =-\frac{1}{2}\frac{\partial a}{\partial T}\Bigg|_{\mathbf{X}},\t{ and}
     * @f]
     * @f[
     * S_{k,T} \equiv \frac{dS_k}{dT}\Bigg|_{P,\mathbf{X}} = -A'_k.
     * @f]
     */
    void getPartialMolarCp(span<double> cpbar) const override;

    //! Get the species molar heat capacities associated with the constant volume
    //! partial molar internal energies [J/kmol/K].
    /*!
     * This method computes
     * @f[
     * \tilde{c}_{v,k}
     *   = \left(\partial \tilde{u}_k/\partial T\right)_{V,\mathbf{n}}.
     * @f]
     * where @f$ \tilde{u}_k @f$ are the species molar internal energies computed at
     * constant volume as in getPartialMolarIntEnergies_TV().
     *
     * Using the [previously-introduced notation](#redlich-Kwong-notation),
     * @f[
     * \tilde{u}_k = u_k^\t{ref}(T) + \frac{1}{b\sqrt{T}} \left(L S_k + b_k F C\right),
     * @f]
     * At constant volume and composition, @f$ L @f$ and @f$ C @f$ are constant and only
     * @f$ F @f$ and @f$ S_k @f$ contribute temperature derivatives. Defining
     * @f[
     * \tilde{u}_k^E \equiv \frac{L S_k + b_k F C}{b\sqrt{T}},
     * @f]
     * Then
     * @f[
     * \tilde{c}_{v,k}
     * = c_{v,k}^\t{ref}(T) + \frac{d\tilde{u}_k^E}{dT}\Bigg|_{V,\mathbf{n}}
     * @f]
     * where @f$ c_{v,k}^\t{ref} @f$ are the ideal gas specific heat capacities and
     * @f[
     * \frac{d\tilde{u}_k^E}{dT}\Bigg|_{V,\mathbf{n}} =
     * \frac{1}{b\sqrt{T}} \left( L \frac{dS_k}{dT}\Bigg|_{\mathbf{X}}
     *                            + b_k C \frac{dF}{dT}\Bigg|_{\mathbf{X}} \right)
     * - \frac{\tilde{u}_k^E}{2T}.
     * @f]
     * For the temperature-dependent form @f$ a_{ij}(T)=a_{ij,0}+a_{ij,1}T @f$
     * used here:
     * @f[
     * \frac{dA_k}{dT}\Bigg|_{\mathbf{X}} = A'_k,\qquad
     * \frac{dS_k}{dT}\Bigg|_{\mathbf{X}} = -A'_k,\qquad
     * \frac{dF}{dT}\Bigg|_{\mathbf{X}} = -\frac{1}{2}
     *     \left.\frac{\partial a}{\partial T}\right|_{\mathbf{X}}.
     * @f]
     */
    void getPartialMolarCv_TV(span<double> cvbar) const override;
    void getPartialMolarVolumes(span<double> vbar) const override;
    //! @}

public:
    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{

    bool addSpecies(shared_ptr<Species> spec) override;
    void initThermo() override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;

    //! Set the pure fluid interaction parameters for a species
    /*!
     *  The "a" parameter for species *i* in the Redlich-Kwong model is assumed
     *  to be a linear function of temperature:
     *  @f[ a = a_0 + a_1 T @f]
     *
     *  @param species   Name of the species
     *  @param a0        constant term in the expression for the "a" parameter
     *      of the specified species [Pa-m^6/kmol^2]
     *  @param a1        temperature-proportional term in the expression for the
     *      "a" parameter of the specified species [Pa-m^6/kmol^2/K]
     *  @param b         "b" parameter in the Redlich-Kwong model [m^3/kmol]
     */
    void setSpeciesCoeffs(const string& species, double a0, double a1, double b);

    //! Set values for the interaction parameter between two species
    /*!
     *  The "a" parameter for interactions between species *i* and *j* is
     *  assumed by default to be computed as:
     *  @f[ a_{ij} = \sqrt(a_{i,0} a_{j,0}) + \sqrt(a_{i,1} a_{j,1}) T @f]
     *
     *  This function overrides the defaults with the specified parameters:
     *  @f[ a_{ij} = a_{ij,0} + a_{ij,1} T @f]
     *
     *  @param species_i   Name of one species
     *  @param species_j   Name of the other species
     *  @param a0          constant term in the "a" expression [Pa-m^6/kmol^2]
     *  @param a1          temperature-proportional term in the "a" expression
     *      [Pa-m^6/kmol^2/K]
     */
    void setBinaryCoeffs(const string& species_i,
                         const string& species_j, double a0, double a1);
    //! @}

protected:
    // Special functions inherited from MixtureFugacityTP
    double sresid() const override;
    double hresid() const override;

public:
    double liquidVolEst(double TKelvin, double& pres) const override;
    double densityCalc(double T, double pressure, int phase, double rhoguess) override;
    double dpdVCalc(double TKelvin, double molarVol, double& presCalc) const override;

    double isothermalCompressibility() const override;
    double thermalExpansionCoeff() const override;
    double internalPressure() const override;
    double soundSpeed() const override;

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
    void updateMixingExpressions() override;

    //! Calculate the a and the b parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param temp  Temperature (TKelvin)
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     */
    void calculateAB(double temp, double& aCalc, double& bCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    double da_dt() const;

    void calcCriticalConditions(double& pc, double& tc, double& vc) const override;

    //! Prepare variables and call the function to solve the cubic equation of state
    int solveCubic(double T, double pres, double a, double b, span<double> Vroot) const;

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

    vector<double> a_vec_Curr_;
    vector<double> b_vec_Curr_;

    Array2D a_coeff_vec;

    //! Explicitly-specified binary interaction parameters
    map<string, map<string, pair<double, double>>> m_binaryParameters;

    enum class CoeffSource { EoS, CritProps, Database };
    //! For each species, specifies the source of the a and b coefficients
    vector<CoeffSource> m_coeffSource;

    int NSolns_ = 0;

    double Vroot_[3] = {0.0, 0.0, 0.0};

    //! Temporary storage - length = m_kk.
    mutable vector<double> m_pp;

    // Partial molar volumes of the species
    mutable vector<double> m_partialMolarVolumes;

    mutable vector<double> m_dAkdT; //!< Temporary storage for dA_k/dT; length #m_kk.

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
    mutable vector<double> dpdni_;

private:
    //! Omega constant for a -> value of a in terms of critical properties
    /*!
     *  this was calculated from a small nonlinear solve
     */
    static const double omega_a;

    //! Omega constant for b
    static const double omega_b;

    //! Omega constant for the critical molar volume
    static const double omega_vc;
};
}

#endif
