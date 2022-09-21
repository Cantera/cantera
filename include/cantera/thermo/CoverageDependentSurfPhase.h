/**
 *  @file CoverageDependentSurfPhase.h
 *  Header for a thermodynamics model of a coverage-dependent surface
 *  phase derived from SurfPhase, applying adsorbate lateral interaction
 *  correction factors to the SurfPhase thermodynamic properties.
 *  (see \ref thermoprops and class
 *  \link Cantera::CoverageDependentSurfPhase CoverageDependentSurfPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_COVERAGEDEPENDENTSURFPHASE_H
#define CT_COVERAGEDEPENDENTSURFPHASE_H

#include "SurfPhase.h"
#include "cantera/base/Array.h"

namespace Cantera
{

//! A thermodynamic model for a coverage-dependent surface phase, applying adsorbate
//! lateral interaction correction factors to the ideal surface phase properties.
/*!
 * The ideal surface phase assumes no lateral interaction among adsorbates.
 * This coverage-dependent surface phase allows adding lateral interaction
 * correction terms to the ideal surface phase (SurfPhase) thermodynamic properties
 * so that more accurate adsorbate thermochemistry can be achieved.
 *
 * ## Coverage-dependent Thermodynamic Properties Formulations
 *
 * At a low-coverage limit, an adsorbate thermochemistry is the same as
 * that of ideal species since there are no adsorbates in the vicinity to cause
 * lateral interaction. Therefore, it is logical to set ideal species properties
 * as the low-coverage limit properties and add lateral interaction terms, i.e.
 * lateral interaction correction terms, to them as excess properties.
 * Accordingly, standard state coverage-dependent enthalpy, entropy,
 * and heat capacity of a species \f$ k \f$ can be formulated as follows.
 *
 * \f[
 *  h_k^o(T,\theta)
 *      = \underbrace{h_k^{o,ideal}(T)
 *        + \int_{298}^{T}cp_k^{o,ideal}(T)dT}_{\text{low-coverage limit}}
 *        + \underbrace{h_k^{o,cov}(T,\theta)
 *        + \int_{298}^{T}cp_k^{o,cov}(T,\theta)dT}_{\text{coverage dependence}}
 *
 * \f]
 *
 * \f[
 *  s_k^o(T,\theta)
 *     = \underbrace{s_k^{o,ideal}(T)
 *       + \int_{298}^{T}\frac{cp_k^{o,ideal}(T)}{T}dT}_{\text{low-coverage limit}}
 *       + \underbrace{s_k^{o,cov}(T,\theta)
 *       + \int_{298}^{T}\frac{cp_k^{o,cov}(T,\theta)}{T}dT}_{\text{coverage
 *         dependence}}
 * \f]
 *
 * \f[
 *  cp_k^o(T,\theta)
 *      = cp_k^{o,ideal}(T) + cp_k^{o,cov}(T,\theta)
 * \f]
 *
 * ## Mathematical Models for Coverage-dependent Correction Terms
 *
 * Coverage-dependent correction terms for enthalpy and entropy can be calculated
 * with one of the four algebraic models: linear dependecy model, piecewise-linear
 * dependency model, polynomial dependency model, and interpolative dependency model.
 * In the dependency model equations, a coverage-dependent correction term is denoted
 * by \f$ f^{cov} \f$ where \f$ f \f$ can be either enthalpy or entropy. Because
 * lateral interaction can compose of both self- and cross- interactions, the total
 * correction term of species \f$ k \f$ is a sum of all interacting species \f$ j \f$
 * which can include itself. Coefficients \f$ c1_{k,j}-c7_{k,j} \f$ are user-provided
 * parameters that can be given in input mechanism.
 *
 * Linear dependency model:
 * \f[
 *  f^{cov}_k(\theta) = \sum_j c1_{k,j} \theta_j
 * \f]
 *
 * Piecewise-linear dependency model:
 * \f[
 * f^{cov}_k(\theta) = \sum_j \left\{
 *  \begin{array}{ll}
 *  c2_{k,j}\theta_j & \text{, } \theta_j \leq \theta^\text{change}_{k,j} \\
 *  \left[c3_{k,j}(\theta_j - \theta^\text{change}_{k,j})
 *  + (c2_{k,j}\theta^\text{change}_{k,j})\right]
 *  & \text{, } \theta_j > \theta^\text{change}_{k,j} \\
 *  \end{array}
 *  \right.
 * \f]
 *
 * Polynomial dependency model:
 * \f[
 *  f^{cov}_k(\theta) =
 *   \sum_j \left[c4_{k,j}\theta_j + c5_{k,j}\theta_j^2
 *           + c6_{k,j}\theta_j^3 + c7_{k,j}\theta_j^4\right]
 * \f]
 *
 * Interpolative dependency model:
 * \f[
 *  f^{cov}_k(\theta) =
 *   \sum_j \left[\frac{f^{cov}_k(\theta^{higher}_j) - f^{cov}_k(\theta^{lower}_j)}
 *   {\theta^{higher}_j - \theta^{lower}_j}(\theta_j - \theta^{lower}_j)
 *   + f^{cov}_k (\theta^{lower}_j)\right] \text{, }
 *   \theta^{lower}_j \leq \theta_j < \theta^{higher}_j
 * \f]
 *
 * Coverage-dependent heat capacity is calculated using an equation with a
 * quadratic dependence on coverages and a logarithmic dependence on temperature.
 * Temperature is nondimesionalized with a reference temperature of 1 K.
 * The coverage-dependent heat capacity of species \f$ k \f$ is a sum of
 * all quantities dependent on coverage of species \f$ j \f$. Coefficients
 * \f$ c8_{k,j} \text{ and } c9_{k,j} \f$ are user-provided parameters that can be
 * given in input mechanism.
 *
 * \f[
 *  cp^{cov}_k(\theta) =
 *   \sum_j \left(c8_{k,j} ln\left(\frac{T}{1\text{ K}}\right)
 *   + c9_{k,j}\right) \theta_j^2
 * \f]
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
         * @param enthalpy_coeffs array of polynomial coefficients describing
         *                        coverage-depdendent enthalpy [J/kmol] in order of
         *                        1st-order, 2nd-order, 3rd-order, and 4th-order
         *                        coefficients
         * @param entropy_coeffs array of polynomial coefficients describing
         *                       coverage-dependent entropy [J/kmol/K] in order of
         *                       1st-order, 2nd-order, 3rd-order, and 4th-order
         *                       coefficients
         */
        PolynomialDependency(size_t k, size_t j,
                             vector_fp enthalpy_coeffs={0.0, 0.0, 0.0, 0.0, 0.0},
                             vector_fp entropy_coeffs={0.0, 0.0, 0.0, 0.0, 0.0}):
                             k(k), j(j), enthalpy_coeffs(enthalpy_coeffs),
                             entropy_coeffs(entropy_coeffs) {}
        //! index of a target species whose enthalpy and entropy is calculated
        size_t k;
        //! index of a species whose coverage affects enthalpy and entropy of
        //! a target species
        size_t j;
        //! array of polynomial coefficients describing coverage-depdendent enthalpy
        //! [J/kmol] in order of 1st-order, 2nd-order, 3rd-order, and 4th-order
        //! coefficients
        vector_fp enthalpy_coeffs;
        //! array of polynomial coefficients describing coverage-depdendent entropy
        //! [J/kmol/K] in order of 1st-order, 2nd-order, 3rd-order, and 4th-order
        //! coefficients
        vector_fp entropy_coeffs;
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
         * @param enthalpy_map map of <coverage[dimensionless], enthalpy[J/kmol]> pairs
         * @param entropy_map map of <coverage[dimensionless], entropy[J/kmol/K]> pairs
         */
        InterpolativeDependency(size_t k, size_t j,
                                std::map<double, double> enthalpy_map={{0.0, 0.0},
                                                                       {1.0, 0.0}},
                                std::map<double, double> entropy_map={{0.0, 0.0},
                                                                      {1.0, 0.0}}):
                                k(k), j(j), enthalpy_map(enthalpy_map),
                                entropy_map(entropy_map) {}
        //! index of a target species whose enthalpy and entropy are calculated
        size_t k;
        //! index of a species whose coverage affects enthalpy and entropy of
        //! a target species
        size_t j;
        //! map of <coverage[dimensionless], enthalpy[J/kmol]> pairs
        std::map<double, double> enthalpy_map;
        //! map of <coverage[dimensionless], entropy[J/kmol/K]> pairs
        std::map<double, double> entropy_map;
    };

    //! A struct to store sets of parameters used in coverage-dependent heat capacity
    //! calculations by a log-quadratic equation in CoverageDependentSurfPhase.
    struct HeatCapacityDependency
    {
        /*!
         * @param k index of a target species whose heat capacity is calculated
         * @param j index of a species whose coverage affects heat capacity of
         *          a target species
         * @param coeff_a coefficient a [J/kmol/K]
         * @param coeff_b coefficient b [J/kmol/K]
         */
        HeatCapacityDependency(size_t k, size_t j,
                               double coeff_a=0.0, double coeff_b=0.0):
                               k(k), j(j), coeff_a(coeff_a), coeff_b(coeff_b) {}
        //! index of a target species whose heat capacity is calculated
        size_t k;
        //! index of a species whose coverage affects heat capacity of
        //! a target species
        size_t j;
        //! coefficient a [J/kmol/K]
        double coeff_a;
        //! coefficient b [J/kmol/K]
        double coeff_b;
    };

    //! Construct and initialize a CoverageDependentSurfPhase ThermoPhase object
    //! directly from an ASCII input file.
    /*!
     * @param infile name of the input file. If blank, an empty phase will be created.
     * @param id     name of the phase id in the file. If blank, the first phase
     *               in the file is used.
     */
    explicit CoverageDependentSurfPhase(const std::string& infile="",
                                        const std::string& id="");

    //! Add interpolative coverage dependece parameters for a species
    /*!
     *  @param int_deps  list of parameters as an InterpolativeDependency object
     */
    void addInterpolativeDependency(const InterpolativeDependency& int_deps);

    virtual void initThermo();
    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void getParameters(AnyMap& phaseNode) const;
    virtual void getSpeciesParameters(const std::string& name,
                                      AnyMap& speciesNode) const;

    //! @name Methods calculating reference state thermodynamic properties
    //! Reference state properties are evaluated at \f$ T \text{ and }
    //! \theta^{ref} \f$. With coverage fixed at a reference value,
    //! reference state properties are effectively only dependent on temperature.
    //! @{
    virtual void getEnthalpy_RT_ref(double* hrt) const;
    virtual void getEntropy_R_ref(double* sr) const;
    virtual void getCp_R_ref(double* cpr) const;
    virtual void getGibbs_RT_ref(double* grt) const;
    //! @}

    //! @name Methods calculating standard state thermodynamic properties
    //! Standard state properties are evaluated at \f$ T \text{ and } \theta \f$,
    //! and thus are dependent both on temperature and coverage.
    //! @{

    //! Get the nondimensionalized standard state enthalpy vector.
    /*!
     * \f[
     *      \frac{h^o_k(T,\theta)}{RT}
     *          = \frac{h^{ref}_k(T) + h^{cov}_k(T,\theta)
     *            + \int_{298}^{T} cp^{cov}_k(T,\theta)dT}{RT}
     * \f]
     */
    virtual void getEnthalpy_RT(double* hrt) const;

    //! Get the nondimensionalized standard state entropy vector.
    /*!
     * \f[
     *      \frac{s^o_k(T,\theta)}{R}
     *          = \frac{s^{ref}_k(T) + s^{cov}_k(T,\theta)
     *            + \int_{298}^{T}\frac{cp^{cov}_k(T,\theta)}{T}dT}{R}
     *            - ln\left(\frac{1}{\theta_{ref}}\right)
     * \f]
     */
    virtual void getEntropy_R(double* sr) const;

    //! Get the nondimensionalized standard state heat capacity vector.
    /*!
     * \f[
     *      \frac{cp^o_k(T,\theta)}{RT}
     *          = \frac{cp^{ref}_k(T) + cp^{cov}_k(T,\theta)}{RT}
     * \f]
     */
    virtual void getCp_R(double* cpr) const;

    //! Get the nondimensionalized standard state gibbs free energy vector.
    /*!
     * \f[
     *      \frac{g^o_k(T,\theta)}{RT}
     *          = \frac{h^o_k(T,\theta)}{RT} + \frac{s^o_k(T,\theta)}{R}
     * \f]
     */
    virtual void getGibbs_RT(double* grt) const;

    //! Get the standard state gibbs free energy vector. Units: J/kmol.
    /*!
     * \f[
     *      g^o_k(T,\theta) = h^o_k(T,\theta) + Ts^o_k(T,\theta)
     * \f]
     */
    virtual void getPureGibbs(double* g) const;

    //! Get the standard state chemical potential vector. Units: J/kmol.
    /*!
     * \f[
     *      \mu^o_k(T,\theta) = h^o_k(T,\theta) + Ts^o_k(T,\theta)
     * \f]
     */
    virtual void getStandardChemPotentials(double* mu0) const;
    //! @}

    //! @name Methods calculating partial molar thermodynamic properties
    //! Partial molar properties are evaluated at \f$ T \text{ and } \theta \f$,
    //! and thus are dependent both on temperature and coverage.
    //! @{

    //! Get the partial molar enthalpy vector. Units: J/kmol.
    /*!
     * \f[
     *      \overline{h}_k(T,\theta) = h^o_k(T,\theta)
     * \f]
     */
    virtual void getPartialMolarEnthalpies(double* hbar) const;

    //! Get the partial molar entropy vector. Units: J/kmol/K.
    /*!
     * \f[
     *      \overline{s}_k(T,\theta) = s^o_k(T,\theta) - Rln(\theta_k)
     * \f]
     */
    virtual void getPartialMolarEntropies(double* sbar) const;

    //! Get the partial molar heat capacity vector. Units: J/kmol/K.
    /*!
     * \f[
     *      \overline{cp}_k(T,\theta) = cp^o_k(T,\theta)
     * \f]
     */
    virtual void getPartialMolarCp(double* cpbar) const;

    //! Get the chemical potential vector. Units: J/kmol.
    /*!
     * \f[
     *      \mu_k(T,\theta) = \mu^o_k(T,\theta) + RTln(\theta_k)
     * \f]
     */
    virtual void getChemPotentials(double* mu) const;
    //! @}

    //! @name Methods calculating Phase thermodynamic properties
    //! Phase properties are evaluated at \f$ T \text{ and } \theta \f$,
    //! and thus are dependent both on temperature and coverage.

    //! @{

    //! Return the solution's molar enthalpy. Units: J/kmol
    /*!
     * \f[
     *      \hat h(T,\theta) = \sum_k \theta_k \overline{h}_k(T,\theta)
     * \f]
     */
    virtual double enthalpy_mole() const;

    //! Return the solution's molar entropy. Units: J/kmol/K
    /*!
     * \f[
     *      \hat s(T,\theta) = \sum_k \theta_k \overline{s}_k(T,\theta)
     * \f]
     */
    virtual double entropy_mole() const;

    //! Return the solution's molar cp. Units: J/kmol/K
    /*!
     * \f[
     *      \hat{cp}(T,\theta) = \sum_k \theta_k \overline{cp}_k(T,\theta)
     * \f]
     */
    virtual double cp_mole() const;
    //! @}

protected:
    //! Temporary storage for the coverages.
    mutable vector_fp m_cov;

    //! Temporary storage for the coverage-dependent enthalpies.
    mutable vector_fp m_h_cov;

    //! Temporary storage for the coverage-dependent entropies.
    mutable vector_fp m_s_cov;

    //! Temporary storage for the coverage-dependent heat capacities.
    mutable vector_fp m_cp_cov;

    //! Temporary storage for the coverage-dependent chemical potentials.
    mutable vector_fp m_mu_cov;

    //! Temporary storage for the sum of reference state enthalpies and
    //! coverage-dependent enthalpies.
    mutable vector_fp m_enthalpy;

    //! Temporary storage for the sum of reference state entropies and
    //! coverage-dependent entropies.
    mutable vector_fp m_entropy;

    //! Temporary storage for the sum of reference state heat capacities and
    //! coverage-dependent heat capacities.
    mutable vector_fp m_heatcapacity;

    //! Temporary storage for the sum of reference state chemical potentials
    //! and coverage-dependent chemical potentials.
    mutable vector_fp m_chempot;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the linear and polynomial dependency equations.
    std::vector<PolynomialDependency> m_PolynomialDependency;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the piecewise-linear and interpolative dependency equations.
    std::vector<InterpolativeDependency> m_InterpolativeDependency;

    //! Array of heat capacity coverage dependency parameters.
    std::vector<HeatCapacityDependency> m_HeatCapacityDependency;

    //! Explicitly-specified k-j interaction parameters, to enable serialization.
    //! For linear model. <name_k>: <name_j: PolynomialDependency_index>>.
    std::map<std::string, std::map<std::string, size_t>> m_indexmap_lin;

    //! Explicitly-specified k-j interaction parameters, to enable serialization.
    //! For polynomial model. <name_k>: <name_j: PolynomialDependency_index>>.
    std::map<std::string, std::map<std::string, size_t>> m_indexmap_poly;

    //! Explicitly-specified k-j interaction parameters, to enable serialization.
    //! For piecewise-linear model. <name_k>: <name_j: InterpolativeDependency_index>>.
    std::map<std::string, std::map<std::string, size_t>> m_indexmap_pwlin;

    //! Explicitly-specified k-j interaction parameters, to enable serialization.
    //! For interpolative model. <name_k>: <name_j: InterpolativeDependency_index>>.
    std::map<std::string, std::map<std::string, size_t>> m_indexmap_int;

    //! Explicitly-specified k-j interaction parameters, to enable serialization.
    //! <name_k>: <name_j: HeatCapacityDependency_index>>.
    std::map<std::string, std::map<std::string, size_t>> m_indexmap_cp;

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
