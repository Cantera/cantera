/**
 *  @file CoverageDependentSurfPhase.h
 *  Header for a simple thermodynamics model of a coverage-dependent surface
 *  phase derived from SurfPhase, assuming an ideal solution model
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

//! Set of parameters modifying SurfPhase enthalpy and entropy based on fractional surface coverages
//! using a polynomial model. Linear model is a subset of the polynomial model.
struct PolynomialDependency
{
    //! Constructor
    //! @param name_k_  name of a species whose thermodynamics are calculated
    //! @param name_j_  name of a species whose coverage affects thermodynamics of the target species
    //! @param enthalpy_coeffs_  array of polynomial coefficients describing enthalpy change
    //!                          containing 1st-order, 2nd-order, 3rd-order, and 4th-order
    //!                          coefficients as a function of coverage [J/kmol]
    //! @param entropy_coeffs_  array of polynomial coefficients describing entropy change
    //!                         containing 1st-order, 2nd-order, 3rd-order, and 4th-order
    //!                         coefficients as a function of coverage [J/kmol/K]
    PolynomialDependency(std::string name_k_, std::string name_j_,
                         vector_fp enthalpy_coeffs_, vector_fp entropy_coeffs_):
                         name_k(name_k_), name_j(name_j_),
                         enthalpy_coeffs(enthalpy_coeffs_),
                         entropy_coeffs(entropy_coeffs_) {}
    PolynomialDependency() {}
    //! name of a species whose thermodynamics are calculated
    std::string name_k;
    //! name of a species whose coverage affects thermodynamics of the target species
    std::string name_j;
    //! array of polynomial coefficients describing enthalpy change containing 1st-order, 2nd-order,
    //! 3rd-order, and 4th-order coefficients as a function of coverage [J/kmol]
    vector_fp enthalpy_coeffs;
    //! array of polynomial coefficients describing entropy change containing 1st-order, 2nd-order,
    //! 3rd-order, and 4th-order coefficients as a function of coverage [J/kmol/K]
    vector_fp entropy_coeffs;
};

//! Set of parameters modifying SurfPhase enthalpy and entropy based on fractional surface coverages
//! using a piecewise linear model.
struct PiecewiseDependency
{
    //! Constructor
    //! @param name_k_  name of a species whose thermodynamics are calculated
    //! @param name_j_  name of a species whose coverage affects thermodynamics of the target species
    //! @param enthalpy_params_  array of three parameters to calculate coverage-dependent enthalpy:
    //!                          slope of enthalpy change in the first region [J/kmol], slope of
    //!                          enthalpy change in the second region [J/kmol], and coverage
    //!                          dividing first and second region [dimensionless]
    //! @param entropy_params_  array of three parameters to calculate coverage-dependent entropy:
    //!                         slope of entropy change in the first region [J/kmol/K], slope of
    //!                         entropy change in the second region [J/kmol/K], and coverage
    //!                         dividing first and second region [dimensionless]
    PiecewiseDependency(std::string name_k_, std::string name_j_,
                        vector_fp enthalpy_params_, vector_fp entropy_params_):
                        name_k(name_k_), name_j(name_j_),
                        enthalpy_params(enthalpy_params_),
                        entropy_params(entropy_params_) {}
    PiecewiseDependency() {}
    //! name of a species whose thermodynamics are calculated
    std::string name_k;
    //! name of a species whose coverage affects thermodynamics of the target species
    std::string name_j;
    //! array of three parameters to calculate coverage-dependent enthalpy: slope of enthalpy change
    //! in the first region [J/kmol], slope of enthalpy change in the second region [J/kmol],
    //! and coverage dividing first and second region [dimensionless]
    vector_fp enthalpy_params;
    //! array of three parameters to calculate coverage-dependent entropy: slope of entropy change
    //! in the first region [J/kmol/K], slope of entropy change in the second region [J/kmol/K],
    //! and coverage dividing first and second region [dimensionless]
    vector_fp entropy_params;
};

//! Set of parameters modifying SurfPhase enthalpy and entropy based on fractional surface coverages
//! using a interpolative model.
struct InterpolativeDependency
{
    //! Constructor
    //! @param name_k_  name of a species whose thermodynamics are calculated
    //! @param name_j_  name of a species whose coverage affects thermodynamics of the target species
    //! @param enthalpy_coverages_  array of coverages for coverage-dependent enthalpy interpolation
    //!                             [dimensionless]
    //! @param enthalpies_ array of enthalpies at corresponding coverages in enthalpy-coverages
    //!                    [J/kmol]
    //! @param entropy_coverages_  array of coverages for coverage-dependent entropy interpolation
    //!                            [dimensionless]
    //! @param entropies_  array of entropies at corresponding coverages in entropy-coverages
    //!                    [J/kmol/K]
    InterpolativeDependency(std::string name_k_, std::string name_j_,
                            vector_fp enthalpy_coverages_, vector_fp enthalpies_,
                            vector_fp entropy_coverages_, vector_fp entropies_):
                            name_k(name_k_), name_j(name_j_),
                            enthalpy_coverages(enthalpy_coverages_),
                            enthalpies(enthalpies_),
                            entropy_coverages(entropy_coverages_),
                            entropies(entropies_) {}
    InterpolativeDependency() {}
    //! name of a species whose thermodynamics are calculated
    std::string name_k;
    //! name of a species whose coverage affects thermodynamics of the target species
    std::string name_j;
    //! array of coverages for coverage-dependent enthalpy interpolation [dimensionless]
    vector_fp enthalpy_coverages;
    //! array of enthalpies at corresponding coverages in enthalpy-coverages [J/kmol]
    vector_fp enthalpies;
    //! array of coverages for coverage-dependent entropy interpolation [dimensionless]
    vector_fp entropy_coverages;
    //! array of entropies at corresponding coverages in entropy-coverages [J/kmol/K]
    vector_fp entropies;
};

//! Set of parameters modifying SurfPhase heat capacity based on fractional surface coverages
//! using a quadratic model.
struct HeatCapacityDependency
{
    //! Constructor
    //! @param name_k_  name of a species whose thermodynamics are calculated
    //! @param name_j_  name of a species whose coverage affects thermodynamics of the target species
    //! @param cpcov_a_  log model coefficient a [J/kmol/K]
    //! @param cpcov_b_  log model coefficient b [J/kmol/K]
    HeatCapacityDependency(std::string name_k_, std::string name_j_,
                           double cpcov_a_, double cpcov_b_):
                           name_k(name_k_), name_j(name_j_),
                           cpcov_a(cpcov_a_), cpcov_b(cpcov_b_) {}
    HeatCapacityDependency() {}
    //! name of a species whose thermodynamics are calculated
    std::string name_k;
    //! name of a species whose coverage affects thermodynamics of the target species
    std::string name_j;
    double cpcov_a; //! log model coefficient a [J/kmol/K]
    double cpcov_b; //! log model coefficient b [J/kmol/K]
};

class CoverageDependentSurfPhase : public SurfPhase
{
public:
    //! Constructor.
    //! @param n0 Site Density of the Surface Phase [kmol m^-2]
    CoverageDependentSurfPhase(double n0 = 1.0);

    //! Construct and initialize a CoverageDependentSurfPhase ThermoPhase object directly
    //! from an input file
    //! @param infile name of the input file
    //! @param id     name of the phase id in the file. If this is blank, the first phase
    //!               in the file is used
    explicit CoverageDependentSurfPhase(const std::string& infile, const std::string& id);

    //! Set the polynomial coverage dependece for species
    /*!
     *  enthalpy and entropy are sum of ideal surface species enthalpy and entropy
     *  and coverage-dependent enthalpy and entropy which are caclulated with
     *  a polynomial function of coverages:
     *
     *  For the linear depedence,
     *  \f[ h^{cov}_k(\theta) = \sum_j h^{slope}_{k,j} \theta_j \f]
     *  \f[ s^{cov}_k(\theta) = \sum_j s^{slope}_{k,j} \theta_j \f]
     *
     *  For the polynomial dependence,
     *  \f[ h^{cov}_k(\theta) =
     *      \sum_j (h^{1st-order}_{k,j} \theta_j + h^{2nd-order}_{k,j} \theta_j^2
     *              + h^{3rd-order}_{k,j} \theta_j^3 + h^{4th-order}_{k,j} \theta_j^4)
     *  \f]
     *  \f[ s^{cov}_k(\theta) =
     *      \sum_j (s^{1st-order}_{k,j} \theta_j + s^{2nd-order}_{k,j} \theta_j^2
     *              + s^{3rd-order}_{k,j} \theta_j^3 + s^{4th-order}_{k,j} \theta_j^4)
     *  \f]
     *  @param poly_deps  list of parameters as a PolynomialDependency object
     */
    void setPolynomialDependency(const PolynomialDependency& poly_deps);

    //! Set the piecewise linear coverage dependece for species
    /*!
     *  enthalpy and entropy are sum of ideal surface species enthalpy and entropy
     *  and coverage-dependent enthalpy and entropy which are caclulated with
     *  a piecewise linear function of coverages:
     *
     *  \f[ h^{cov}_k(\theta) = \sum_j h^{low}_{k,j} \theta_j
     *      \text{, for } \theta_j \leq \theta^{change}_{k,j}
     *  \f]
     *  \f[ h^{cov}_k(\theta) = \sum_j h^{high}_{k,j}
     *      (\theta_j - \theta^{change}_{k,j}) +
     *      (h^{low}_{k,j} \theta^{change}_{k,j}) \text{, for }
     *      \theta_j > \theta^{change}_{k,j}
     *  \f]
     *
     *  \f[ s^{cov}_k(\theta) = \sum_j s^{low}_{k,j} \theta_j
     *      \text{, for } \theta_j \leq \theta^\text{change}_{k,j}
     *  \f]
     *  \f[ s^{cov}_k(\theta) = \sum_j s^{high}_{k,j}
     *      (\theta_j - \theta^{change}_{k,j}) +
     *      (s^{low}_{k,j} \theta^{change}_{k,j}) \text{, for }
     *      \theta_j > \theta^{change}_{k,j}
     *  \f]
     *  @param plin_deps  list of parameters as a PiecewiseDependency object
     */
    void setPiecewiseDependency(const PiecewiseDependency& plin_deps);

    //! Set the interpolative coverage dependece for species
    /*!
     *  enthalpy and entropy are sum of ideal surface species enthalpy and entropy
     *  and coverage-dependent enthalpy and entropy which are caclulated with
     *  a linearly interpolated function of coverages:
     *
     *  \f[ h^{cov}_k(\theta) = \sum_j \frac{h^{right}_{k,j}
     *      - h^{left}_{k,j}}{\theta^{right}_{k,j} - \theta^{left}_{k,j}}
     *      (\theta_j - \theta^{left}_{k,j}) + h^{left}_{k,j}
     *      \text{, where } \theta^{left}_{k,j} \leq \theta_j < \theta^{right}_{k,j}
     *  \f]
     *
     *  \f[ s^{cov}_k(\theta) = \sum_j \frac{s^{right}_{k,j}
     *      - s^{left}_{k,j}}{\theta^{right}_{k,j} - \theta^{left}_{k,j}}
     *      (\theta_j - \theta^{left}_{k,j}) + s^{left}_{k,j}
     *      \text{, where } \theta^{left}_{k,j} \leq \theta_j < \theta^{right}_{k,j}
     *  \f]
     *
     *  @param int_deps  list of parameters as an InterpolativeDependency object
     */
    void setInterpolativeDependency(const InterpolativeDependency& int_deps);

    //! Set the heat capacity coverage dependece for species
    /*!
     *  heat capacity is sum of ideal surface species heat capacity and coverage-dependent
     *  heat capacity which is caculated using a function with quadratic depedence on coverages
     *  and a logarithimic dependence on temperature:

     *  \f[ cp^{cov}_k(\theta) = \sum_j (a_{k,j} ln(T) + b_{k,j}) \theta_j^2 \f]
     *
     *  @param cpcov_deps  list of parameters as a HeatCapacityDependency object
     */
    void setHeatCapacityDependency(const HeatCapacityDependency& cpcov_deps);

    virtual void initThermo();
    virtual bool addSpecies(shared_ptr<Species> spec);

    //! @copydoc SurfPhase::setCoverages
    //! After setting coverages, covstateNum is incremented by 1.
    virtual void setCoverages(const double* theta);

    //! @copydoc SurfPhase::setCoveragesNoNorm
    //! After setting coverages, covstateNum is incremented by 1.
    virtual void setCoveragesNoNorm(const double* theta);

    //! Return the covstateNum
    int statecovNumber() const {
        return m_covstateNum;
    }

    //! Convert given eneregy to the quantity with its default unit of J/kmol.
    double convertEnergy(double value, const std::string& src) const;

    //! Convert given energy per temperature to the quantity with its default unit of J/kmol/K.
    double convertEnergy_T(double value, const std::string& src) const;

    // Functions calculating reference state thermodyanmic properties--------------

    //! Return the nondimensionalized reference state enthalpy.
    /*!
     * Nondimensionalized reference state enthalpy is evaluated at T, $P_{ref}$ and $\theta_{ref}$:
     *
     * \f[
     *      h^{ref}_k / RT = h_k(T, P_{ref}, \theta_{ref}) / RT
     * \f]
     */
    virtual void getEnthalpy_RT_ref(double* hrt) const;
    //! Return the nondimensionalized reference state entropy.
    /*!
     * Nondimensionalized reference state entropy is evaluated at T, P_{ref} and $\theta_{ref}$:
     *
     * \f[
     *      s^{ref}_k / R = s_k(T, P_{ref}, \theta_\text{ref}) / R
     * \f]
     */
    virtual void getEntropy_R_ref(double* sr) const;
    //! Return the nondimensionalized reference state cp.
    /*!
     * Nondimensionalized reference state cp is evaluated at T, P_{ref} and $\theta_{ref}$:
     *
     * \f[
     *      cp^{ref}_k / R = cp_k(T, P_{ref}, \theta_\text{ref}) / R
     * \f]
     */
    virtual void getCp_R_ref(double* cpr) const;
    //! Return the nondimensionalized reference state Gibbs free energy.
    /*!
     * \f[
     *      g^{ref}_k / RT = h^{ref}_k /RT - s^{ref}_k /R
     * \f]
     */
    virtual void getGibbs_RT_ref(double* grt) const;

    // Functions calculating standard state thermodyanmic properties---------------

    //! Return the nondimensionalized standard state enthalpy.
    /*!
     * Nondimensionalized standard state enthalpy is evaluated at T, P and $\theta$:
     *
     * \f[
     *      h^o_k / RT = (h^{ref}_k + h^{cov}_k(T, P, \theta)) / RT
     *                   + \int_{T_{ref}}^{T} cp^{cov}_k(T, P, \theta) dT / RT
     * \f]
     */
    virtual void getEnthalpy_RT(double* hrt) const;
    //! Return the nondimensionalized standard state entropy.
    /*!
     * Nondimensionalized standard state entropy is evaluated at T, P and $\theta$:
     *
     * \f[
     *      s^o_k / R = (s^{ref}_k + s^{cov}_k(T, P, \theta)) / R - ln(1 / \theta_{ref})
     *                   + \int_{T_{ref}}^{T} cp^{cov}_k(T, P, \theta) / T dT / R
     * \f]
     */
    virtual void getEntropy_R(double* sr) const;
    //! Return the nondimensionalized standard state cp.
    /*!
     * Nondimensionalized standard state cp is evaluated at T, P and $\theta$:
     *
     * \f[
     *      cp^o_k / RT = (cp^{ref}_k + cp^{cov}_k(T, P, \theta)) / RT
     * \f]
     */
    virtual void getCp_R(double* cpr) const;
    //! Return the nondimensionalized standard state Gibbs free energy.
    /*!
     * \f[
     *      g^o_k / RT = h^o_k / RT - s^o_k / R
     * \f]
     */
    virtual void getGibbs_RT(double* grt) const;
    //! Return the standard state Gibbs free energy. Units: J/kmol
    /*!
     * \f[
     *      g^o_k = h^o_k - T s^o_k
     * \f]
     */
    virtual void getPureGibbs(double* g) const;
    //! Return the standard state chemical potential. Units: J/kmol
    /*!
     * \f[
     *      \mu^o_k = h^o_k - T s^o_k
     * \f]
     */
    virtual void getStandardChemPotentials(double* mu0) const;

    // Functions calculating partial molar thermodyanmic properties----------------

    //! Return the partial molar enthalpy. Units: J/kmol
    /*!
     * Partial molar enthalpy is evaluated at T, P and $\theta$:
     *
     * \f[
     *      \overline{h}_k = h^o_k(T, P, \theta)
     * \f]
     */
    virtual void getPartialMolarEnthalpies(double* hbar) const;
    //! Return the partial molar entropy. Units: J/kmol/K
    /*!
     * Partial molar entropy is evaluated at T, P and $\theta$:
     *
     * \f[
     *      \overline{s}_k = s^o_k(T, P, \theta) - R ln(\theta)
     * \f]
     */
    virtual void getPartialMolarEntropies(double* sbar) const;
    //! Return the partial molar cp. Units: J/kmol/K
    /*!
     * Partial molar cp is evaluated at T, P and $\theta$:
     *
     * \f[
     *      \overline{cp}_k = cp^o_k(T, P, \theta)
     * \f]
     */
    virtual void getPartialMolarCp(double* cpbar) const;

    //! Return the chemical potential. Units: J/kmol
    /*!
     * \f[
     *      \mu_k = mu^o_k(T, P, \theta) + RT ln(\theta)
     * \f]
     */
    virtual void getChemPotentials(double* mu) const;

    // Functions calculating mixture thermodyanmic properties----------------------

    //! Return the solution's molar enthalpy. Units: J/kmol
    /*!
     * Asumming an ideal solution,
     * \f[
     * \hat h(T, P, \theta) = \sum_k \theta_k \overline{h}_k(T, \theta)
     *                      = \sum_k \theta_k h^o_k(T, \theta)
     * \f]
     *
     * \see MultiSpeciesThermo
     */
    virtual double enthalpy_mole() const;

    //! Return the solution's molar entropy. Units: J/kmol/K
    /*!
     * \f[
     * \hat s(T, P, \theta) = \sum_k \theta_k \overline{s}_k(T, \theta)
     *                      = \sum_k \theta_k s^o_k(T, \theta) - R ln(\theta)
     * \f]
     */
    virtual double entropy_mole() const;

    //! Return the solution's molar cp. Units: J/kmol/K
    /*!
     * Asumming an ideal solution,
     * \f[
     * \hat{cp} (T, P, \theta) = \sum_k \theta_k \overline{cp}_k(T, \theta)
     *                         = \sum_k \theta_k cp^o_k(T, \theta)
     * \f]
     *
     * \see MultiSpeciesThermo
     */
    virtual double cp_mole() const;

protected:
    //! Temporary storage for the reference state enthalpies.
    mutable vector_fp m_h_ref;

    //! Temporary storage for the reference state entropies.
    mutable vector_fp m_s_ref;

    //! Temporary storage for the reference state heat capacities.
    mutable vector_fp m_cp_ref;

    //! Temporary storage for the reference state chemical potentials.
    mutable vector_fp m_mu_ref;

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
    //! coverage-dependent heat capapcities.
    mutable vector_fp m_heatcapacity;

    //! Temporary storage for the sum of reference state chemical potentials
    //! and coverage-dependent chemical potentials.
    mutable vector_fp m_chempot;

    /*! Indicate the surface coverages have changed by incrementing the
     *  coverage state number by one. This function is called whenever
     *  the setCoverages() or setCoveragesNoNorm() was used to update the
     *  coverage array.
     */
    void coverageChanged();

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the linear and polynomial models.
    mutable std::vector<PolynomialDependency> m_PolynomialDependency;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the piecewise linear model.
    mutable std::vector<PiecewiseDependency> m_PiecewiseDependency;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the interpolative model.
    mutable std::vector<InterpolativeDependency> m_InterpolativeDependency;

    //! Array of heat capacity coverage dependency parameters.
    mutable std::vector<HeatCapacityDependency> m_HeatCapacityDependency;

    //! Storage for the linear and polynomial model coverage dependency
    //! parameters for enthalpy.
    mutable vector_fp m_polynomial_h;

    //! Storage for the linear and polynomial model coverage dependency
    //! parameters for entropy.
    mutable vector_fp m_polynomial_s;

    //! Storage for the piecewise linear model coverage dependency parameters
    //! for enthalpy.
    mutable vector_fp m_piecewise_h;

    //! Storage for the piecewise linear model coverage dependency parameters
    //! for entropy.
    mutable vector_fp m_piecewise_s;

    //! Storage for the interpolative model coverage array for enthalpy.
    mutable vector_fp m_interpolative_hcov;

    //! Storage for the interpolative model enthalpy array.
    mutable vector_fp m_interpolative_h;

    //! Storage for the interpolative model coverage array for entropy.
    mutable vector_fp m_interpolative_scov;

    //! Storage for the interpolative model entropy array.
    mutable vector_fp m_interpolative_s;

    //! Storage for the coverage-dependent heat capacity log model coefficients.
    mutable double m_cpcov_a, m_cpcov_b;

private:
    //! Coverage state change variable. Whenever the coverage array changes,
    //! this number is incremented.
    int m_covstateNum;

    //! Last value of the coverage state number processed.
    mutable int m_covstateNumlast;

    //! Flag indicating polynomial coverage-dependent model is being used.
    mutable bool m_has_polynomial_dependency;

    //! Flag indicating piecewise linear coverage-dependent model is being used.
    mutable bool m_has_piecewise_dependency;

    //! Flag indicating interpolative coverage-dependent model is being used.
    mutable bool m_has_interpolative_dependency;

    //! Flag indicating coverage-dependent heat capacity is being calculated.
    mutable bool m_has_heatcapacity_dependency;

    //! Flag indicating user has defined the reference coverage.
    mutable bool m_has_ref_coverage;

    //! Storage for the user-defined reference state coverage which has to be
    //! greater than 0.0 and less than or equal to 1.0.
    mutable double m_theta_ref;

    //! Update the species ideal reference state thermodynamic functions
    /*!
     * The polynomials for the ideal surface reference functions are only
     * reevaluated if the temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevaluation of the thermo
     *               polynomials. default = false.
     */
    void _updateReferenceThermo(bool force=false) const;

    //! Update the species coverage-dependent thermodynamic functions
    /*!
     * The coverage-dependent enthalpy and entropy are only reevaluated
     * if the coverage has changed. The coverage-dependent heat capacity
     * is only reevaluated if the coverage or temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevaluation of the
     *               coverage-dependent functions. default = false.
     */
    void _updateCovDepThermo(bool force=false) const;

    //! Update the total (reference state + coverage-dependent) thermodynamic functions
    /*!
     * Calls subroutines for ideal species thermodynamic update as well as
     * coverage-dependent species thermodynamic update
     *
     * @param force  Boolean, which if true, forces a reevaluation of the
     *               reference state as well as coverage-dependent functions.
     *               default = false.
     */
    void _updateThermo(bool force=false) const;

};

}

#endif
