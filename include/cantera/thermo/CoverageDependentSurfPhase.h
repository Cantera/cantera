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
 *  The ideal surface phase assumes no lateral interaction among adsorbates.
 *  This coverage-dependent surface phase allows adding lateral interaction
 *  correction terms to the ideal surface phase thermodynamic properties so that
 *  more accurate adsorbate thermochemistry can be achieved.
 *
 *  ## Coverage-dependent Thermodynamic Properties Formulations
 *
 *  Naturally, at a coverage limit an adsorbate thermochemistry is same as the
 *  that of ideal species since there are no near adsorbates to render lateral
 *  interactions. Therefore, it is logical to set ideal species properties
 *  as the low coverage quantities and apply lateral interaction terms as excess
 *  quantities. For example, a species coverage-dependent enthalpy of formation
 *  can be formulated as follows.
 *  \f[
 *       \Delta_f h_k^{\circ,298} = \underbrace{ \Delta_f h_k^{\circ,298,
 *           \text{ideal species}}}_{\text{low-coverage limit}}
 *           + \underbrace{{\Delta_f h_k^{\circ,298}}^E
 *           \left( \theta \right)}_{\text{excess enthalpy of formation}} \nonumber
 *  \f]
 *
 *  ## Mathematical Models for Coverage-dependent Correction Terms
 *
 *  The coverage-dependent correction terms can be given as one of the four
 *  functional models: linear dependecy model, piecewise-linear dependency model,
 *  polynomial dependency model, and interpolative dependency model.
 *
 */
class CoverageDependentSurfPhase : public SurfPhase
{
public:
    //! Default constructor
    CoverageDependentSurfPhase();

    //! Construct and initialize a CoverageDependentSurfPhase ThermoPhase object
    //! directly from an input file
    //! @param infile name of the input file
    //! @param id     name of the phase id in the file. If this is blank, the first
    //!               phase in the file is used
    explicit CoverageDependentSurfPhase(const std::string& infile,
                                        const std::string& id="");

    //! Set of parameters modifying SurfPhase enthalpy and entropy based on
    //! surface coverages using a polynomial model.
    //! Linear model is a subset of the polynomial model.
    struct PolynomialDependency
    {
        //! Constructor
        //! @param k_  index of a target species whose thermodynamics are calculated
        //! @param j_  index of a species whose coverage affects thermodynamics of
        //!            a target species
        //! @param enthalpy_coeffs_  array of polynomial coefficients describing
        //!                          enthalpy change [J/kmol] containing 1st-order,
        //!                          2nd-order, 3rd-order, and 4th-order coefficients
        //!                          as a function of coverage
        //! @param entropy_coeffs_  array of polynomial coefficients describing entropy
        //!                         change [J/kmol/K] containing 1st-order, 2nd-order,
        //!                         3rd-order, and 4th-order coefficients as a function
        //!                         of coverage
        PolynomialDependency(size_t k_, size_t j_,
                            vector_fp enthalpy_coeffs_, vector_fp entropy_coeffs_):
                            k(k_), j(j_),
                            enthalpy_coeffs(enthalpy_coeffs_),
                            entropy_coeffs(entropy_coeffs_) {}
        PolynomialDependency() {}
        //! index of a target species whose thermodynamics are calculated
        size_t k;
        //! index of a species whose coverage affects thermodynamics of a target species
        size_t j;
        //! array of polynomial coefficients describing enthalpy change [J/kmol] containing
        //! 1st-order, 2nd-order, 3rd-order, and 4th-order coefficients as a function
        //! of coverage
        vector_fp enthalpy_coeffs;
        //! array of polynomial coefficients describing entropy change [J/kmol/K]
        //! containing 1st-order, 2nd-order, 3rd-order, and 4th-order coefficients as
        //! a function of coverage
        vector_fp entropy_coeffs;
    };

    //! Set of parameters modifying SurfPhase enthalpy and entropy based on
    //! surface coverages using a interpolative model.
    //! Piecewise model is a subset of interpolative with three data points.
    struct InterpolativeDependency
    {
        //! Constructor
        //! @param k_  index of a target species whose thermodynamics are calculated
        //! @param j_  index of a species whose coverage affects thermodynamics of
        //!            a target species
        //! @param enthalpy_map_  map of <coverage, enthalpy> as a key-value pair of
        //!                       coverage [dimensionless] and enthalpy [J/kmol]
        //! @param entropy_map_  map of <coverage, entropy> as a key-value pair of
        //!                      coverage [dimensionless] and entropy [J/kmol/K]
        InterpolativeDependency(size_t k_, size_t j_,
                                std::map<double, double> enthalpy_map_,
                                std::map<double, double> entropy_map_):
                                k(k_), j(j_),
                                enthalpy_map(enthalpy_map_),
                                entropy_map(entropy_map_){}
        InterpolativeDependency() {}
        //! index of a target species whose thermodynamics are calculated
        size_t k;
        //! index of a species whose coverage affects thermodynamics of a target species
        size_t j;
        //! map of coverage-enthalpy pairs for coverage-dependent enthalpy interpolation
        std::map<double, double> enthalpy_map;
        //! map of coverage-entropy pairs for coverage-dependent entropy interpolation
        std::map<double, double> entropy_map;
    };

    //! Set of parameters modifying SurfPhase heat capacity based on surface coverages
    //! using a quadratic model.
    struct HeatCapacityDependency
    {
        //! Constructor
        //! @param k_  index of a target species whose thermodynamics are calculated
        //! @param j_  index of a species whose coverage affects thermodynamics of
        //!            a target species
        //! @param cpcov_a_  log model coefficient a [J/kmol/K]
        //! @param cpcov_b_  log model coefficient b [J/kmol/K]
        HeatCapacityDependency(size_t k_, size_t j_,
                            double cpcov_a_, double cpcov_b_):
                            k(k_), j(j_),
                            cpcov_a(cpcov_a_), cpcov_b(cpcov_b_) {}
        HeatCapacityDependency() {}
        //! index of a target species whose thermodynamics are calculated
        size_t k;
        //! index of a species whose coverage affects thermodynamics of a target
        //! species
        size_t j;
        double cpcov_a; //! log model coefficient a [J/kmol/K]
        double cpcov_b; //! log model coefficient b [J/kmol/K]
    };

    //! Set the polynomial coverage dependece for species
    /*!
     *  enthalpy and entropy are sum of ideal surface species enthalpy and entropy
     *  and coverage-dependent enthalpy and entropy which are calculated with
     *  a polynomial function of coverages:
     *
     *  For the linear dependence,
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

    //! Set the interpolative coverage dependece for species
    /*!
     *  enthalpy and entropy are sum of ideal surface species enthalpy and entropy
     *  and coverage-dependent enthalpy and entropy which are calculated with
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
     *  heat capacity is sum of ideal surface species heat capacity and
     *  coverage-dependent heat capacity which is calculated using a function with
     *  quadratic dependence on coverages and a logarithmic dependence on temperature:

     *  \f[ cp^{cov}_k(\theta) = \sum_j (a_{k,j} ln(T) + b_{k,j}) \theta_j^2 \f]
     *
     *  @param cpcov_deps  list of parameters as a HeatCapacityDependency object
     */
    void setHeatCapacityDependency(const HeatCapacityDependency& cpcov_deps);

    virtual void initThermo();
    virtual bool addSpecies(shared_ptr<Species> spec);

    //! @name Functions calculating reference state thermodynamic properties
    //! @{
    //! Get the nondimensionalized reference state enthalpy.
    /*!
     * Nondimensionalized reference state enthalpy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      h^{ref}_k / RT = h_k(T, P_{ref}, \theta_{ref}) / RT
     * \f]
     */
    virtual void getEnthalpy_RT_ref(double* hrt) const;

    //! Get the nondimensionalized reference state entropy.
    /*!
     * Nondimensionalized reference state entropy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      s^{ref}_k / R = s_k(T, P_{ref}, \theta_\text{ref}) / R
     * \f]
     */
    virtual void getEntropy_R_ref(double* sr) const;

    //! Get the nondimensionalized reference state cp.
    /*!
     * Nondimensionalized reference state cp is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      cp^{ref}_k / R = cp_k(T, P_{ref}, \theta_\text{ref}) / R
     * \f]
     */
    virtual void getCp_R_ref(double* cpr) const;

    //! Return the nondimensionalized reference state Gibbs free energy.
    /*!
     * Nondimensionalized reference state Gibbs free energy is evaluated
     * at \f$ T, P_{ref} \f$ and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      g^{ref}_k / RT = h^{ref}_k /RT - s^{ref}_k /R
     * \f]
     */
    virtual void getGibbs_RT_ref(double* grt) const;
    //! @}

    //! @name Functions calculating standard state thermodynamic properties
    //! @{
    //! Get the nondimensionalized standard state enthalpy.
    /*!
     * Nondimensionalized standard state enthalpy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      h^o_k / RT = (h^{ref}_k + h^{cov}_k(T, P, \theta)) / RT
     *                   + \int_{T_{ref}}^{T} cp^{cov}_k(T, P, \theta) dT / RT
     * \f]
     */
    virtual void getEnthalpy_RT(double* hrt) const;

    //! Get the nondimensionalized standard state entropy.
    /*!
     * Nondimensionalized standard state entropy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      s^o_k / R = (s^{ref}_k + s^{cov}_k(T, P, \theta))
     *                  / R - ln(1 / \theta_{ref})
     *                  + \int_{T_{ref}}^{T} cp^{cov}_k(T, P, \theta) / T dT / R
     * \f]
     */
    virtual void getEntropy_R(double* sr) const;

    //! Get the nondimensionalized standard state cp.
    /*!
     * Nondimensionalized standard state cp is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      cp^o_k / RT = (cp^{ref}_k + cp^{cov}_k(T, P, \theta)) / RT
     * \f]
     */
    virtual void getCp_R(double* cpr) const;

    //! Get the nondimensionalized standard state Gibbs free energy.
    /*!
     * Nondimensionalized standard state Gibbs free energy evaluated
     * at \f$ T, P_{ref} \f$ and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      g^o_k / RT = h^o_k / RT - s^o_k / R
     * \f]
     */
    virtual void getGibbs_RT(double* grt) const;

    //! Get the standard state Gibbs free energy. Units: J/kmol
    /*!
     * Standard state Gibbs free energy evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      g^o_k = h^o_k - T s^o_k
     * \f]
     */
    virtual void getPureGibbs(double* g) const;

    //! Get the standard state chemical potential. Units: J/kmol
    /*!
     * Standard state chemical potential evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      \mu^o_k = h^o_k - T s^o_k
     * \f]
     */
    virtual void getStandardChemPotentials(double* mu0) const;
    //! @}

    //! @name Functions calculating partial molar thermodynamic properties
    //! @{
    //! Get the partial molar enthalpy. Units: J/kmol
    /*!
     * Partial molar enthalpy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      \overline{h}_k = h^o_k(T, P, \theta)
     * \f]
     */
    virtual void getPartialMolarEnthalpies(double* hbar) const;
    //! Get the partial molar entropy. Units: J/kmol/K
    /*!
     * Partial molar entropy is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      \overline{s}_k = s^o_k(T, P, \theta) - R ln(\theta)
     * \f]
     */
    virtual void getPartialMolarEntropies(double* sbar) const;

    //! Get the partial molar cp. Units: J/kmol/K
    /*!
     * Partial molar cp is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      \overline{cp}_k = cp^o_k(T, P, \theta)
     * \f]
     */
    virtual void getPartialMolarCp(double* cpbar) const;

    //! Get the chemical potential. Units: J/kmol
    /*!
     * Chemical potential is evaluated at \f$ T, P_{ref} \f$
     * and \f$ \theta_{ref} \f$:
     *
     * \f[
     *      \mu_k = mu^o_k(T, P, \theta) + RT ln(\theta)
     * \f]
     */
    virtual void getChemPotentials(double* mu) const;
    //! @}

    //! @name Functions calculating mixture thermodynamic properties
    //! @{
    //! Return the solution's molar enthalpy. Units: J/kmol
    /*!
     * Assuming an ideal mixing,
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
     * Assuming an ideal mixing,
     * \f[
     * \hat{cp} (T, P, \theta) = \sum_k \theta_k \overline{cp}_k(T, \theta)
     *                         = \sum_k \theta_k cp^o_k(T, \theta)
     * \f]
     *
     * \see MultiSpeciesThermo
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
    //! the linear and polynomial models.
    std::vector<PolynomialDependency> m_PolynomialDependency;

    //! Array of enthalpy and entropy coverage dependency parameters used in
    //! the interpolative model.
    std::vector<InterpolativeDependency> m_InterpolativeDependency;

    //! Array of heat capacity coverage dependency parameters.
    std::vector<HeatCapacityDependency> m_HeatCapacityDependency;

private:
    //! Last value of the state number processed.
    mutable int m_stateNumlast;

    //! Storage for the user-defined reference state coverage which has to be
    //! greater than 0.0 and less than or equal to 1.0. default = 1.0.
    mutable double m_theta_ref;

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
    void _updateTotalThermo(bool force=false) const;

};

}

#endif
