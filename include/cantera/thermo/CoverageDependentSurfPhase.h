/**
 *  @file CoverageDependentSurfPhase.h
 *  Header for a thermodynamics model of a coverage-dependent surface
 *  phase derived from SurfPhase, applying adsorbate lateral interaction
 *  correction factors to the SurfPhase thermodynamic properties.
 *  (see @ref thermoprops and class
 *  @link Cantera::CoverageDependentSurfPhase CoverageDependentSurfPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_COVERAGEDEPENDENTSURFPHASE_H
#define CT_COVERAGEDEPENDENTSURFPHASE_H

#include "SurfPhase.h"
#include "cantera/base/Array.h"

namespace Cantera
{

//! A thermodynamic model for a coverage-dependent surface phase, applying surface
//! species lateral interaction correction factors to the ideal surface phase
//! properties.
/*!
 * The ideal surface phase assumes no lateral interaction among surface species.
 * This coverage-dependent surface phase allows adding lateral interaction
 * correction terms to the ideal surface phase (SurfPhase) thermodynamic properties
 * so that more accurate surface species thermochemistry can be achieved.
 *
 * ## Coverage-dependent Thermodynamic Properties Formulations
 *
 * At a low-coverage limit, a surface species thermochemistry is the same as
 * that of ideal surface species since there are no adsorbates in the vicinity
 * to cause lateral interaction. Therefore, it is logical to set ideal surface
 * species properties as the low-coverage limit and add lateral interaction terms
 * to them as excess properties. Accordingly, standard state coverage-dependent
 * enthalpy, entropy, and heat capacity of a surface species @f$ k @f$ can be
 * formulated as follows.
 *
 * @f[
 *  h_k^o(T,\theta)
 *      = \underbrace{h_k^{o,ideal}(T)
 *        + \int_{298}^{T}c_{p,k}^{o,ideal}(T)dT}_{\text{low-coverage limit}}
 *        + \underbrace{h_k^{o,cov}(T,\theta)
 *        + \int_{298}^{T}c_{p,k}^{o,cov}(T,\theta)dT}_{\text{coverage dependence}}
 *
 * @f]
 *
 * @f[
 *  s_k^o(T,\theta)
 *     = \underbrace{s_k^{o,ideal}(T)
 *       + \int_{298}^{T}\frac{c_{p,k}^{o,ideal}(T)}{T}dT}_{\text{low-coverage limit}}
 *       + \underbrace{s_k^{o,cov}(T,\theta)
 *       + \int_{298}^{T}\frac{c_{p,k}^{o,cov}(T,\theta)}{T}dT}_{\text{coverage
 *         dependence}}
 * @f]
 *
 * @f[
 *  c_{p,k}^o(T,\theta)
 *      = \underbrace{c_{p,k}^{o,ideal}(T)}_{\text{low-coverage limit}}
 *        + \underbrace{c_{p,k}^{o,cov}(T,\theta)}_{\text{coverage dependence}}
 * @f]
 *
 * ## Mathematical Models for Coverage-dependent Correction Terms
 *
 * Coverage-dependent correction terms for enthalpy and entropy can be calculated
 * with one of the four algebraic models: linear dependency model, polynomial
 * dependency model, piecewise-linear, and interpolative dependency model.
 * In the dependency model equations, a coverage-dependent correction term is denoted
 * by @f$ f^{cov} @f$ where @f$ f @f$ can be either enthalpy (@f$ h^{cov} @f$) or
 * entropy (@f$ s^{cov} @f$). Because lateral interaction can compose of both
 * self- and cross- interactions, the total correction term of species @f$ k @f$
 * is a sum of all interacting species @f$ j @f$ which can include itself.
 * Coefficients @f$ c^{(1)}_{k,j}-c^{(6)}_{k,j} @f$ are user-provided parameters
 * that can be given in a input yaml.
 *
 * Linear dependency model:
 * @f[
 *  f^{cov}_k(\theta) = \sum_j c^{(1)}_{k,j} \theta_j
 * @f]
 *
 * Polynomial dependency model:
 * @f[
 *  f^{cov}_k(\theta) =
 *   \sum_j \left[c^{(1)}_{k,j}\theta_j + c^{(2)}_{k,j}\theta_j^2
 *           + c^{(3)}_{k,j}\theta_j^3 + c^{(4)}_{k,j}\theta_j^4\right]
 * @f]
 *
 * Piecewise-linear dependency model:
 * @f[
 * f^{cov}_k(\theta) = \sum_j \left\{
 *  \begin{array}{ll}
 *  c^{(5)}_{k,j}\theta_j & \text{, } \theta_j \leq \theta^\text{change}_{k,j} \\
 *  \left[c^{(6)}_{k,j}(\theta_j - \theta^\text{change}_{k,j})
 *  + (c^{(5)}_{k,j}\theta^\text{change}_{k,j})\right]
 *  & \text{, } \theta_j > \theta^\text{change}_{k,j} \\
 *  \end{array}
 *  \right.
 * @f]
 *
 * Interpolative dependency model:
 * @f[
 *  f^{cov}_k(\theta) =
 *   \sum_j \left[\frac{f^{cov}_k(\theta^{higher}_j) - f^{cov}_k(\theta^{lower}_j)}
 *   {\theta^{higher}_j - \theta^{lower}_j}(\theta_j - \theta^{lower}_j)
 *   + f^{cov}_k (\theta^{lower}_j)\right] \\
 *   \text{where } \theta^{lower}_j \leq \theta_j < \theta^{higher}_j
 * @f]
 *
 * Coverage-dependent heat capacity is calculated using an equation with a
 * quadratic dependence on coverages and a logarithmic dependence on temperature.
 * Temperature is nondimensionalized with a reference temperature of 1 K.
 * The coverage-dependent heat capacity of species @f$ k @f$ is a sum of
 * all quantities dependent on coverage of species @f$ j @f$. Coefficients
 * @f$ c^{(a)}_{k,j} \text{ and } c^{(b)}_{k,j} @f$ are user-provided parameters
 * that can be given in an input yaml.
 *
 * Coverage-dependent heat capacity model:
 * @f[
 *  c^{cov}_{p,k}(\theta) =
 *   \sum_j \left(c^{(a)}_{k,j} \ln\left(\frac{T}{1\text{ K}}\right)
 *   + c^{(b)}_{k,j}\right) \theta_j^2
 * @f]
 */
class CoverageDependentSurfPhase : public SurfPhase
{
public:
    //! A struct to store sets of parameters used in coverage-dependent enthalpy
    //! and entropy calculations by a polynomial equation or a linear equation
    //! in CoverageDependentSurfPhase.
    struct PolynomialDependency
    {
        /*!
         * @param k index of a target species whose enthalpy and entropy are calculated
         * @param j index of a species whose coverage affects enthalpy and entropy of
         *          a target species
         * @param dep_map map of coverage-dependency parameters
         */
        PolynomialDependency(size_t k, size_t j, const AnyMap& dep_map);

        //! index of a target species whose enthalpy and entropy is calculated
        size_t k;
        //! index of a species whose coverage affects enthalpy and entropy of
        //! a target species
        size_t j;
        //! array of polynomial coefficients describing coverage-dependent enthalpy
        //! [J/kmol] in order of 1st-order, 2nd-order, 3rd-order, and 4th-order
        //! coefficients (@f$ c^{(1)}, c^{(2)}, c^{(3)}, \text{ and } c^{(4)} @f$
        //! in the linear or the polynomial dependency model)
        vector<double> enthalpy_coeffs;
        //! array of polynomial coefficients describing coverage-dependent entropy
        //! [J/kmol/K] in order of 1st-order, 2nd-order, 3rd-order, and 4th-order
        //! coefficients (@f$ c^{(1)}, c^{(2)}, c^{(3)}, \text{ and } c^{(4)} @f$
        //! in the linear or the polynomial dependency model)
        vector<double> entropy_coeffs;
        //! boolean indicating whether the dependency is linear
        bool isLinear;
    };

    //! A struct to store sets of parameters used in coverage-dependent enthalpy
    //! and entropy calculations by a interpolative equation or a piecewise-linear
    //! equation in CoverageDependentSurfPhase.
    struct InterpolativeDependency
    {
        /*!
         * @param k index of a target species whose enthalpy and entropy are calculated
         * @param j index of a species whose coverage affects enthalpy and entropy of
         *          a target species
         * @param dep_map map of coverage-dependency parameters
         * @param node species node of a target species
         */
        InterpolativeDependency(size_t k, size_t j,
                                const AnyMap& dep_map, const AnyBase& node);

        //! index of a target species whose enthalpy and entropy are calculated
        size_t k;
        //! index of a species whose coverage affects enthalpy and entropy of
        //! a target species
        size_t j;
        //! map of <coverage[dimensionless], enthalpy[J/kmol]> pairs
        map<double, double> enthalpy_map;
        //! map of <coverage[dimensionless], entropy[J/kmol/K]> pairs
        map<double, double> entropy_map;
        //! boolean indicating whether the dependency is piecewise-linear
        bool isPiecewise;
    };

    //! A struct to store sets of parameters used in coverage-dependent heat capacity
    //! calculations by a log-quadratic equation in CoverageDependentSurfPhase.
    struct HeatCapacityDependency
    {
        /*!
         * @param k index of a target species whose heat capacity is calculated
         * @param j index of a species whose coverage affects heat capacity of
         *          a target species
         */
        HeatCapacityDependency(size_t k, size_t j):
                               k(k), j(j), coeff_a(0.0), coeff_b(0.0) {}
        //! index of a target species whose heat capacity is calculated
        size_t k;
        //! index of a species whose coverage affects heat capacity of
        //! a target species
        size_t j;
        //! coefficient @f$ c^{(a)} @f$ [J/kmol/K] in the coverage-dependent
        //! heat capacity model
        double coeff_a;
        //! coefficient @f$ c^{(b)} @f$ [J/kmol/K] in the coverage-dependent
        //! heat capacity model
        double coeff_b;
    };

    //! Construct and initialize a CoverageDependentSurfPhase ThermoPhase object
    //! directly from an ASCII input file.
    /*!
     * @param infile name of the input file. If blank, an empty phase will be created.
     * @param id     name of the phase id in the file. If blank, the first phase
     *               in the file is used.
     */
    explicit CoverageDependentSurfPhase(const string& infile="",
                                        const string& id="");

    string type() const override {
        return "coverage-dependent-surface";
    }

    //! Add interpolative coverage dependence parameters for a species
    /*!
     *  @param int_deps  list of parameters as an InterpolativeDependency object
     */
    void addInterpolativeDependency(const InterpolativeDependency& int_deps);

    void initThermo() override;
    bool addSpecies(shared_ptr<Species> spec) override;
    void getParameters(AnyMap& phaseNode) const override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;

    //! @name Methods calculating reference state thermodynamic properties
    //! Reference state properties are evaluated at @f$ T \text{ and }
    //! \theta^{ref} @f$. With coverage fixed at a reference value,
    //! reference state properties are effectively only dependent on temperature.
    //! @{
    void getEnthalpy_RT_ref(double* hrt) const override;
    void getEntropy_R_ref(double* sr) const override;
    void getCp_R_ref(double* cpr) const override;
    void getGibbs_RT_ref(double* grt) const override;
    //! @}

    //! @name Methods calculating standard state thermodynamic properties
    //! Standard state properties are evaluated at @f$ T \text{ and } \theta @f$,
    //! and thus are dependent both on temperature and coverage.
    //! @{

    //! Get the nondimensionalized standard state enthalpy vector.
    /*!
     * @f[
     *      \frac{h^o_k(T,\theta)}{RT}
     *          = \frac{h^{ref}_k(T) + h^{cov}_k(T,\theta)
     *            + \int_{298}^{T} c^{cov}_{p,k}(T,\theta)dT}{RT}
     * @f]
     */
    void getEnthalpy_RT(double* hrt) const override;

    //! Get the nondimensionalized standard state entropy vector.
    /*!
     * @f[
     *      \frac{s^o_k(T,\theta)}{R}
     *          = \frac{s^{ref}_k(T) + s^{cov}_k(T,\theta)
     *            + \int_{298}^{T}\frac{c^{cov}_{p,k}(T,\theta)}{T}dT}{R}
     *            - \ln\left(\frac{1}{\theta_{ref}}\right)
     * @f]
     */
    void getEntropy_R(double* sr) const override;

    //! Get the nondimensionalized standard state heat capacity vector.
    /*!
     * @f[
     *      \frac{c^o_{p,k}(T,\theta)}{RT}
     *          = \frac{c^{ref}_{p,k}(T) + c^{cov}_{p,k}(T,\theta)}{RT}
     * @f]
     */
    void getCp_R(double* cpr) const override;

    //! Get the nondimensionalized standard state gibbs free energy vector.
    /*!
     * @f[
     *      \frac{g^o_k(T,\theta)}{RT}
     *          = \frac{h^o_k(T,\theta)}{RT} + \frac{s^o_k(T,\theta)}{R}
     * @f]
     */
    void getGibbs_RT(double* grt) const override;

    //! Get the standard state gibbs free energy vector. Units: J/kmol.
    /*!
     * @f[
     *      g^o_k(T,\theta) = h^o_k(T,\theta) + Ts^o_k(T,\theta)
     * @f]
     */
    void getPureGibbs(double* g) const override;

    //! Get the standard state chemical potential vector. Units: J/kmol.
    /*!
     * @f[
     *      \mu^o_k(T,\theta) = h^o_k(T,\theta) + Ts^o_k(T,\theta)
     * @f]
     */
    void getStandardChemPotentials(double* mu0) const override;
    //! @}

    //! @name Methods calculating partial molar thermodynamic properties
    //! Partial molar properties are evaluated at @f$ T \text{ and } \theta @f$,
    //! and thus are dependent both on temperature and coverage.
    //! @{

    //! Get the partial molar enthalpy vector. Units: J/kmol.
    /*!
     * @f[
     *      \tilde{h}_k(T,\theta) = h^o_k(T,\theta)
     * @f]
     */
    void getPartialMolarEnthalpies(double* hbar) const override;

    //! Get the partial molar entropy vector. Units: J/kmol/K.
    /*!
     * @f[
     *      \tilde{s}_k(T,\theta) = s^o_k(T,\theta) - R\ln(\theta_k)
     * @f]
     */
    void getPartialMolarEntropies(double* sbar) const override;

    //! Get the partial molar heat capacity vector. Units: J/kmol/K.
    /*!
     * @f[
     *      \tilde{c}_{p,k}(T,\theta) = c^o_{p,k}(T,\theta)
     * @f]
     */
    void getPartialMolarCp(double* cpbar) const override;

    //! Get the chemical potential vector. Units: J/kmol.
    /*!
     * @f[
     *      \mu_k(T,\theta) = \mu^o_k(T,\theta) + RT\ln(\theta_k)
     * @f]
     */
    void getChemPotentials(double* mu) const override;
    //! @}

    //! @name Methods calculating Phase thermodynamic properties
    //! Phase properties are evaluated at @f$ T \text{ and } \theta @f$,
    //! and thus are dependent both on temperature and coverage.

    //! @{

    //! Return the solution's molar enthalpy. Units: J/kmol
    /*!
     * @f[
     *      \hat h(T,\theta) = \sum_k \theta_k \tilde{h}_k(T,\theta)
     * @f]
     */
    double enthalpy_mole() const override;

    //! Return the solution's molar entropy. Units: J/kmol/K
    /*!
     * @f[
     *      \hat s(T,\theta) = \sum_k \theta_k \tilde{s}_k(T,\theta)
     * @f]
     */
    double entropy_mole() const override;

    //! Return the solution's molar heat capacity. Units: J/kmol/K
    /*!
     * @f[
     *      \hat{c_p}(T,\theta) = \sum_k \theta_k \tilde{c_p}_k(T,\theta)
     * @f]
     */
    double cp_mole() const override;
    //! @}

protected:
    //! Temporary storage for the coverages.
    mutable vector<double> m_cov;

    //! Temporary storage for the coverage-dependent enthalpies.
    mutable vector<double> m_h_cov;

    //! Temporary storage for the coverage-dependent entropies.
    mutable vector<double> m_s_cov;

    //! Temporary storage for the coverage-dependent heat capacities.
    mutable vector<double> m_cp_cov;

    //! Temporary storage for the coverage-dependent chemical potentials.
    mutable vector<double> m_mu_cov;

    //! Temporary storage for the sum of reference state enthalpies and
    //! coverage-dependent enthalpies.
    mutable vector<double> m_enthalpy;

    //! Temporary storage for the sum of reference state entropies and
    //! coverage-dependent entropies.
    mutable vector<double> m_entropy;

    //! Temporary storage for the sum of reference state heat capacities and
    //! coverage-dependent heat capacities.
    mutable vector<double> m_heatcapacity;

    //! Temporary storage for the sum of reference state chemical potentials
    //! and coverage-dependent chemical potentials.
    mutable vector<double> m_chempot;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the linear and polynomial dependency equations.
    vector<PolynomialDependency> m_PolynomialDependency;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the piecewise-linear and interpolative dependency equations.
    vector<InterpolativeDependency> m_InterpolativeDependency;

    //! Array of heat capacity coverage dependency parameters.
    vector<HeatCapacityDependency> m_HeatCapacityDependency;

private:
    //! Storage for the user-defined reference state coverage which has to be
    //! greater than 0.0 and less than or equal to 1.0. default = 1.0.
    double m_theta_ref;

    //! Last value of the state number processed.
    mutable int m_stateNumlast;

    //! Update the species coverage-dependent thermodynamic functions.
    /*!
     * The coverage-dependent enthalpy and entropy are only re-evaluated
     * if the coverage has changed. The coverage-dependent heat capacity
     * is only re-evaluated if the coverage or temperature has changed.
     */
    void _updateCovDepThermo() const;

    //! Update the total (reference state + coverage-dependent)
    //! thermodynamic functions.
    /*!
     * Calls subroutines for both ideal species thermodynamic update and
     * coverage-dependent species thermodynamic update.
     */
    void _updateTotalThermo() const;

};

}

#endif
