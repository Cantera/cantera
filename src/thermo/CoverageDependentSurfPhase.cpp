/**
 *  @file CoverageDependentSurfPhase.cpp
 *  Definitions for a simple thermodynamic model of a coverage-dependent
 *  surface phase derived from SurfPhase, assuming an ideal solution model
 *  (see \ref thermoprops and class
 *  \link Cantera::CoverageDependentSurfPhase CoverageDependentSurfPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/CoverageDependentSurfPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/Species.h"

using namespace std;

namespace Cantera
{
CoverageDependentSurfPhase::CoverageDependentSurfPhase(double n0):
    m_covstateNum(-1),
    m_covstateNumlast(0),
    m_has_polynomial_dependency(false),
    m_has_piecewise_dependency(false),
    m_has_interpolative_dependency(false),
    m_has_heatcapacity_dependency(false),
    m_has_ref_coverage(false)
{
    setSiteDensity(n0);
    setNDim(2);
}

CoverageDependentSurfPhase::CoverageDependentSurfPhase(const std::string& infile,
                                                       const std::string& id_):
    m_covstateNum(-1),
    m_covstateNumlast(0),
    m_has_polynomial_dependency(false),
    m_has_piecewise_dependency(false),
    m_has_interpolative_dependency(false),
    m_has_heatcapacity_dependency(false),
    m_has_ref_coverage(false)
{
    initThermoFile(infile, id_);
}

void CoverageDependentSurfPhase::setPolynomialDependency(const PolynomialDependency&
                                                         poly_deps)
{
    size_t k = speciesIndex(poly_deps.name_k);
    size_t j = speciesIndex(poly_deps.name_j);

    if (k == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setPolynomialDependency",
            "Unknown species '{}'.", poly_deps.name_k);
    }

    if (j == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setPolynomialDependency",
            "Unknown species '{}'.", poly_deps.name_k);
    }

    m_PolynomialDependency.push_back(poly_deps);
    m_has_polynomial_dependency = true;
}

void CoverageDependentSurfPhase::setPiecewiseDependency(const PiecewiseDependency&
                                                        plin_deps)
{
    size_t k = speciesIndex(plin_deps.name_k);
    size_t j = speciesIndex(plin_deps.name_j);
    double hcov_change = plin_deps.enthalpy_params[2];
    double scov_change = plin_deps.entropy_params[2];

    if (k == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setPiecewiseDependency",
            "Unknown species '{}'.", plin_deps.name_k);
    }

    if (j == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setPiecewiseDependency",
            "Unknown species '{}'.", plin_deps.name_j);
    }

    if (hcov_change <= 0.0 || hcov_change > 1.0 || scov_change <= 0.0 || scov_change > 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::setPiecewiseDependency",
            "Coverage where slope changes must be greater than 0.0 and less than or equal to 1.0.");
    }

    m_PiecewiseDependency.push_back(plin_deps);
    m_has_piecewise_dependency = true;
}

void CoverageDependentSurfPhase::setInterpolativeDependency(const
                                                            InterpolativeDependency&
                                                            int_deps)
{
    size_t k = speciesIndex(int_deps.name_k);
    size_t j = speciesIndex(int_deps.name_j);

    if (k == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Unknown species '{}'.", int_deps.name_k);
    }
    if (int_deps.enthalpy_coverages.size() != int_deps.enthalpies.size()) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Sizes of coverages array and enthalpies array are not equal.");
    }
    double hcov_last = 0.0;
    for (const auto& hcov_now : int_deps.enthalpy_coverages) {
        if (hcov_now < hcov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
                "Coverages are not in ascending order.");
        }
        hcov_last = hcov_now;
    }
    if (int_deps.enthalpy_coverages.front() != 0.0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The first element of enthalpy-coverages array must be 0.0.");
    }
    if (int_deps.enthalpy_coverages.back() != 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The last element of enthalpy-coverages array must be 1.0.");
    }

    if (j == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Unknown species '{}'.", int_deps.name_j);
    }
    if (int_deps.entropy_coverages.size() != int_deps.entropies.size()) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Sizes of coverages array and entropies array are not equal.");
    }
    double scov_last = 0.0;
    for (const auto& scov_now : int_deps.entropy_coverages) {
        if (scov_now < scov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
                "Coverages are not in ascending order.");
        }
        scov_last = scov_now;
    }
    if (int_deps.entropy_coverages.front() != 0.0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The first element of entropy-coverages array must be 0.0.");
    }
    if (int_deps.entropy_coverages.back() != 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The last element of entropy-coverages array must be 1.0.");
    }

    m_InterpolativeDependency.push_back(int_deps);
    m_has_interpolative_dependency = true;
}

void CoverageDependentSurfPhase::setHeatCapacityDependency(const HeatCapacityDependency&
                                                           cpcov_deps)
{
    size_t k = speciesIndex(cpcov_deps.name_k);
    size_t j = speciesIndex(cpcov_deps.name_j);

    if (k == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setHeatCapacityDependency",
            "Unknown species '{}'.", cpcov_deps.name_k);
    }

    if (j == npos) {
        throw CanteraError("CoverageDependentSurfPhase::setHeatCapacityDependency",
            "Unknown species '{}'.", cpcov_deps.name_j);
    }

    m_HeatCapacityDependency.push_back(cpcov_deps);
    m_has_heatcapacity_dependency = true;
}

void CoverageDependentSurfPhase::initThermo()
{
    if (m_input.hasKey("site-density")) {
        // Units are kmol/m^2 for surface phases or kmol/m for edge phases
        setSiteDensity(m_input.convert("site-density",
            Units(1.0, 0, -static_cast<double>(m_ndim), 0, 0, 0, 1)));
    }
    if (m_input.hasKey("reference-state-coverage")) {
        m_theta_ref = m_input["reference-state-coverage"].as<double>();
        if (m_theta_ref <= 0.0 || m_theta_ref > 1.0) {
            throw CanteraError("CoverageDependentSurfPhase::initThermo",
               "Reference state coverage must be greater than 0.0 and less than or equal to 1.0.");
        }
        m_has_ref_coverage = true;
    }
    for (auto& item : m_species) {
        // Read enthalpy and entropy dependencies from species 'input' information
        // (i.e. as specified in a YAML input file) for both self- and cross-
        // interactions.
        if (item.second->input.hasKey("coverage-dependencies")) {
            auto cov_map = item.second->input["coverage-dependencies"];
            for (const auto& item2 : cov_map) {
                auto cov_map2 = item2.second.as<AnyMap>();
                // For linear and polynomial model
                if (cov_map2["model"] == "Linear" || cov_map2["model"] == "Polynomial") {
                    m_polynomial_h = {0.0, 0.0, 0.0, 0.0};
                    m_polynomial_s = {0.0, 0.0, 0.0, 0.0};
                    if (cov_map2.hasKey("enthalpy")) {
                        auto enthalpy = cov_map2["enthalpy"].as<double>();
                        m_polynomial_h[0] = enthalpy;
                    } else {
                        if (cov_map2.hasKey("enthalpy-1st-order")) {
                            auto h_a = cov_map2["enthalpy-1st-order"].as<double>();
                            m_polynomial_h[0] = h_a;
                        }
                        if (cov_map2.hasKey("enthalpy-2nd-order")) {
                            auto h_b = cov_map2["enthalpy-2nd-order"].as<double>();
                            m_polynomial_h[1] = h_b;
                        }
                        if (cov_map2.hasKey("enthalpy-3rd-order")) {
                            auto h_c = cov_map2["enthalpy-3rd-order"].as<double>();
                            m_polynomial_h[2] = h_c;
                        }
                        if (cov_map2.hasKey("enthalpy-4th-order")) {
                            auto h_d = cov_map2["enthalpy-4th-order"].as<double>();
                            m_polynomial_h[3] = h_d;
                        }
                    }
                    if (cov_map2.hasKey("enthalpy-unit")) {
                        auto enthalpy_unit = cov_map2["enthalpy-unit"].as<string>();
                        for (size_t i = 0; i < 4; i++) {
                            m_polynomial_h[i] = convertEnergy(m_polynomial_h[i], enthalpy_unit);
                        }
                    }

                    if (cov_map2.hasKey("entropy")) {
                        auto entropy = cov_map2["entropy"].as<double>();
                        m_polynomial_s[0] = entropy;
                    } else {
                        if (cov_map2.hasKey("entropy-1st-order")) {
                            auto s_a = cov_map2["entropy-1st-order"].as<double>();
                            m_polynomial_s[0] = s_a;
                        }
                        if (cov_map2.hasKey("entropy-2nd-order")) {
                            auto s_b = cov_map2["entropy-2nd-order"].as<double>();
                            m_polynomial_s[1] = s_b;
                        }
                        if (cov_map2.hasKey("entropy-3rd-order")) {
                            auto s_c = cov_map2["entropy-3rd-order"].as<double>();
                            m_polynomial_s[2] = s_c;
                        }
                        if (cov_map2.hasKey("entropy-4th-order")) {
                            auto s_d = cov_map2["entropy-4th-order"].as<double>();
                            m_polynomial_s[3] = s_d;
                        }
                    }
                    if (cov_map2.hasKey("entropy-unit")) {
                        auto entropy_unit = cov_map2["entropy-unit"].as<string>();
                        for (size_t i = 0; i < 4; i++) {
                            m_polynomial_s[i] = convertEnergy_T(m_polynomial_s[i], entropy_unit);
                        }
                    }

                    PolynomialDependency poly_deps(item.first, item2.first,
                                                   m_polynomial_h, m_polynomial_s);
                    setPolynomialDependency(poly_deps);
                // For piecewise linear model
                } else if (cov_map2["model"] == "Piecewise-Linear") {
                    m_piecewise_h = {0.0, 0.0, 0.5};
                    m_piecewise_s = {0.0, 0.0, 0.5};
                    if (cov_map2.hasKey("enthalpy-low")) {
                        auto enthalpy_low = cov_map2["enthalpy-low"].as<double>();
                        auto enthalpy_high = cov_map2["enthalpy-high"].as<double>();
                        auto enthalpy_change = cov_map2["enthalpy-change"].as<double>();
                        if (cov_map2.hasKey("enthalpy-unit")) {
                            auto enthalpy_unit = cov_map2["enthalpy-unit"].as<string>();
                            enthalpy_low = convertEnergy(enthalpy_low, enthalpy_unit);
                            enthalpy_high = convertEnergy(enthalpy_high, enthalpy_unit);
                        }
                        m_piecewise_h[0] = enthalpy_low;
                        m_piecewise_h[1] = enthalpy_high;
                        m_piecewise_h[2] = enthalpy_change;
                    }

                    if (cov_map2.hasKey("entropy-low")) {
                        auto entropy_low = cov_map2["entropy-low"].as<double>();
                        auto entropy_high = cov_map2["entropy-high"].as<double>();
                        auto entropy_change = cov_map2["entropy-change"].as<double>();
                        if (cov_map2.hasKey("entropy-unit")) {
                            auto entropy_unit = cov_map2["entropy-unit"].as<string>();
                            entropy_low = convertEnergy_T(entropy_low, entropy_unit);
                            entropy_high = convertEnergy_T(entropy_high, entropy_unit);
                        }
                        m_piecewise_s[0] = entropy_low;
                        m_piecewise_s[1] = entropy_high;
                        m_piecewise_s[2] = entropy_change;
                    }

                    PiecewiseDependency plin_deps(item.first, item2.first,
                                                  m_piecewise_h,
                                                  m_piecewise_s);
                    setPiecewiseDependency(plin_deps);
                // For interpolative model
                } else if (cov_map2["model"] == "Interpolative") {
                    m_interpolative_hcov = {0.0, 1.0};
                    m_interpolative_h = {0.0, 0.0};
                    m_interpolative_scov = {0.0, 1.0};
                    m_interpolative_s = {0.0, 0.0};
                    if (cov_map2.hasKey("enthalpy-coverages")) {
                        m_interpolative_hcov = cov_map2["enthalpy-coverages"].as<vector_fp>();
                        auto enthalpies = cov_map2["enthalpies"].as<vector_fp>();
                        if (cov_map2.hasKey("enthalpy-unit")) {
                            m_interpolative_h.resize(enthalpies.size());
                            auto enthalpy_unit = cov_map2["enthalpy-unit"].as<string>();
                            for (size_t k = 0; k < enthalpies.size(); k++) {
                                m_interpolative_h[k] = convertEnergy(enthalpies[k], enthalpy_unit);
                            }
                        } else {
                            m_interpolative_h = enthalpies;
                        }
                    }

                    if (cov_map2.hasKey("entropy-coverages")) {
                        m_interpolative_scov = cov_map2["entropy-coverages"].as<vector_fp>();
                        auto entropies = cov_map2["entropies"].as<vector_fp>();
                        if (cov_map2.hasKey("entropy-unit")) {
                            m_interpolative_s.resize(entropies.size());
                            auto entropy_unit = cov_map2["entropy-unit"].as<string>();
                            for (size_t k = 0; k < entropies.size(); k++) {
                                m_interpolative_s[k] = convertEnergy_T(entropies[k], entropy_unit);
                            }
                        } else {
                            m_interpolative_s = entropies;
                        }
                    }

                    InterpolativeDependency int_deps(item.first, item2.first,
                                                     m_interpolative_hcov,
                                                     m_interpolative_h,
                                                     m_interpolative_scov,
                                                     m_interpolative_s);
                    setInterpolativeDependency(int_deps);
                } else {
                    throw CanteraError("CoverageDependentSurfPhase::initThermo",
                        "Unrecognized coverage-dependencies model between '{}' and '{}'.",
                        item.first, item2.first);
                }
                // For coverage-dependent heat capacity parameters, if present
                if (cov_map2.hasKey("heat-capacity-a")) {
                    auto cpcov_a = cov_map2["heat-capacity-a"].as<double>();
                    auto cpcov_b = cov_map2["heat-capacity-b"].as<double>();
                    if (cov_map2.hasKey("heat-capacity-unit")){
                        auto cpcov_unit = cov_map2["heat-capacity-unit"].as<string>();
                        m_cpcov_a = convertEnergy_T(cpcov_a, cpcov_unit);
                        m_cpcov_b = convertEnergy_T(cpcov_b, cpcov_unit);
                    } else {
                        m_cpcov_a = cpcov_a;
                        m_cpcov_b = cpcov_b;
                    }
                    HeatCapacityDependency cpcov_deps(item.first, item2.first,
                                                      m_cpcov_a, m_cpcov_b);
                    setHeatCapacityDependency(cpcov_deps);
                }
            }
        }
    }
}

bool CoverageDependentSurfPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = SurfPhase::addSpecies(spec);
    if (added) {
        m_h_ref.push_back(0.0);
        m_s_ref.push_back(0.0);
        m_cp_ref.push_back(0.0);
        m_mu_ref.push_back(0.0);
        m_h_cov.push_back(0.0);
        m_s_cov.push_back(0.0);
        m_cp_cov.push_back(0.0);
        m_mu_cov.push_back(0.0);
        m_enthalpy.push_back(0.0);
        m_entropy.push_back(0.0);
        m_heatcapacity.push_back(0.0);
        m_chempot.push_back(0.0);
    }
    return added;
}

void CoverageDependentSurfPhase::setCoverages(const double* theta)
{
    SurfPhase::setCoverages(theta);
    coverageChanged();
}

void CoverageDependentSurfPhase::setCoveragesNoNorm(const double* theta)
{
    SurfPhase::setCoveragesNoNorm(theta);
    coverageChanged();
}

double CoverageDependentSurfPhase::convertEnergy(double value,
                                                 const std::string& src) const
{
    // Convert to J/kmol
    Units usrc(src);
    if (usrc.convertible(Units("J/kmol"))) {
        value *= usrc.factor();
    } else if (src == "eV" || src == "meV" ) {
        value *= Avogadro * usrc.factor();
    } else {
        throw CanteraError("CoverageDependentSurfPhase::convertEnergy",
            "Don't understand units '{}' as energy", src);
    }
    return value;
}

double CoverageDependentSurfPhase::convertEnergy_T(double value,
                                                  const std::string& src) const
{
    // Convert to J/kmol/K
    Units usrc(src);
    if (usrc.convertible(Units("J/kmol/K"))) {
        value *= usrc.factor();
    } else if (src == "eV/K" || src == "meV/K") {
        value *= Avogadro * usrc.factor();
    } else {
        throw CanteraError("CoverageDependentSurfPhase::convertEnergy_T",
            "Don't understand units '{}' as energy over temperature", src);
    }
    return value;
}

void CoverageDependentSurfPhase::coverageChanged() {
    m_covstateNum++;
}

// Functions calculating reference state thermodyanmic properties--------------

void CoverageDependentSurfPhase::getGibbs_RT_ref(double* grt) const
{
    _updateReferenceThermo();
    scale(m_mu_ref.begin(), m_mu_ref.end(), grt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEnthalpy_RT_ref(double* hrt) const
{
    _updateReferenceThermo();
    scale(m_h_ref.begin(), m_h_ref.end(), hrt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R_ref(double* sr) const
{
    _updateReferenceThermo();
    scale(m_s_ref.begin(), m_s_ref.end(), sr, 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getCp_R_ref(double* cpr) const
{
    _updateReferenceThermo();
    scale(m_cp_ref.begin(), m_cp_ref.end(), cpr, 1.0/GasConstant);
}

// Functions calculating standard state thermodyanmic properties---------------

void CoverageDependentSurfPhase::getEnthalpy_RT(double* hrt) const
{
    _updateThermo();
    scale(m_enthalpy.begin(), m_enthalpy.end(), hrt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R(double* sr) const
{
    _updateThermo();
    scale(m_entropy.begin(), m_entropy.end(), sr, 1.0/GasConstant);
    if (m_has_ref_coverage) {
        double tmp = -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            sr[k] -= tmp;
        }
    }
}

void CoverageDependentSurfPhase::getCp_R(double* cpr) const
{
    _updateThermo();
    scale(m_heatcapacity.begin(), m_heatcapacity.end(), cpr, 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getGibbs_RT(double* grt) const
{
    _updateThermo();
    scale(m_chempot.begin(), m_chempot.end(), grt, 1.0/RT());
    if (m_has_ref_coverage) {
        double tmp = -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            grt[k] += tmp;
        }
    }
}

void CoverageDependentSurfPhase::getPureGibbs(double* g) const
{
    getGibbs_RT(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void CoverageDependentSurfPhase::getStandardChemPotentials(double* mu0) const
{
    _updateThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu0);
    if (m_has_ref_coverage) {
        double tmp = RT() * -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            mu0[k] += tmp;
        }
    }
}

// Functions calling partial molar thermodyanmic properties----------------

void CoverageDependentSurfPhase::getPartialMolarEnthalpies(double* hbar) const
{
    _updateThermo();
    copy(m_enthalpy.begin(), m_enthalpy.end(), hbar);
}

void CoverageDependentSurfPhase::getPartialMolarEntropies(double* sbar) const
{
    _updateThermo();
    copy(m_entropy.begin(), m_entropy.end(), sbar);
    if (m_has_ref_coverage) {
        for (size_t k = 0; k < m_kk; k++) {
            sbar[k] -= GasConstant * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber) / m_theta_ref);
        }
    } else {
        for (size_t k = 0; k < m_kk; k++) {
            sbar[k] -= GasConstant * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber));
        }
    }
}

void CoverageDependentSurfPhase::getPartialMolarCp(double* cpbar) const
{
    _updateThermo();
    copy(m_heatcapacity.begin(), m_heatcapacity.end(), cpbar);
}


void CoverageDependentSurfPhase::getChemPotentials(double* mu) const
{
    _updateThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu);
    if (m_has_ref_coverage) {
        for (size_t k = 0; k < m_kk; k++) {
            mu[k] += RT() * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber) / m_theta_ref);
        }
    } else {
        for (size_t k = 0; k < m_kk; k++) {
            mu[k] += RT() * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber));
        }
    }
}

// Functions calculating mixture thermodyanmic properties--------------------------

double CoverageDependentSurfPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(m_enthalpy);
}

double CoverageDependentSurfPhase::entropy_mole() const
{
    _updateThermo();
    double entropy = 0.0;
    if (m_has_ref_coverage) {
        for (size_t k = 0; k < m_kk; k++) {
            entropy += moleFraction(k) * (m_entropy[k] -
                GasConstant * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber) / m_theta_ref));
        }
    } else {
        for (size_t k = 0; k < m_kk; k++) {
            entropy += moleFraction(k) * (m_entropy[k] -
                GasConstant * log(std::max(concentration(k) * size(k)/m_n0,
                SmallNumber)));
        }
    }
    return entropy;
}

double CoverageDependentSurfPhase::cp_mole() const
{
    _updateThermo();
    return mean_X(m_heatcapacity);
}

void CoverageDependentSurfPhase::_updateReferenceThermo(bool force) const
{
    double tnow = temperature();
    if (m_tlast != tnow || force) {
        m_spthermo.update(tnow, m_cp_ref.data(), m_h_ref.data(), m_s_ref.data());
        for (size_t k = 0; k < m_kk; k++) {
            m_h_ref[k] *= GasConstant * tnow;
            m_s_ref[k] *= GasConstant;
            m_cp_ref[k] *= GasConstant;
            m_mu_ref[k] = m_h_ref[k] - tnow * m_s_ref[k];
        }
    }
}

void CoverageDependentSurfPhase::_updateCovDepThermo(bool force) const
{
    int covstateNumnow = statecovNumber();
    double tnow = temperature();
    if (m_covstateNumlast != covstateNumnow || force) {
        for (size_t k = 0; k < m_kk; k++) {
            m_h_cov[k] = 0.0;
            m_s_cov[k] = 0.0;
        }
        vector_fp cov(m_kk, 0.0);
        SurfPhase::getCoverages(cov.data());

        // For linear and polynomial model
        if (m_has_polynomial_dependency) {
            for (auto& item : m_PolynomialDependency) {
                size_t k = speciesIndex(item.name_k);
                size_t j = speciesIndex(item.name_j);
                for (size_t i = 0; i < 4; i++) {
                    double h_coeff = item.enthalpy_coeffs[i];
                    if (h_coeff != 0.0) {
                        double expo = double (i+1);
                        m_h_cov[k] += h_coeff * pow(cov[j], expo);
                    }
                    double s_coeff = item.entropy_coeffs[i];
                    if (s_coeff != 0.0) {
                        double expo = double (i+1);
                        m_s_cov[k] += s_coeff * pow(cov[j], expo);
                    }
                }
            }
        }

        // For piecewise linear model
        if (m_has_piecewise_dependency) {
            for (auto& item : m_PiecewiseDependency) {
                size_t k = speciesIndex(item.name_k);
                size_t j = speciesIndex(item.name_j);
                double h_slope_low = item.enthalpy_params[0];
                double h_slope_high = item.enthalpy_params[1];
                double h_cov_change = item.enthalpy_params[2];
                if (cov[j] <= h_cov_change) {
                    m_h_cov[k] += h_slope_low * cov[j];
                } else {
                    m_h_cov[k] += h_slope_high
                        * (cov[j] - h_cov_change)
                        + (h_cov_change * h_slope_low);
                }
                double s_slope_low = item.entropy_params[0];
                double s_slope_high = item.entropy_params[1];
                double s_cov_change = item.entropy_params[2];
                if (cov[j] <= s_cov_change) {
                    m_s_cov[k] += s_slope_low * cov[j];
                } else {
                    m_s_cov[k] += s_slope_high
                        * (cov[j] - s_cov_change)
                        + (s_cov_change * s_slope_low);
                }
            }
        }

        // For interpolative model
        if (m_has_interpolative_dependency) {
            for (auto& item : m_InterpolativeDependency) {
                size_t k = speciesIndex(item.name_k);
                size_t j = speciesIndex(item.name_j);
                vector_fp h_covs = item.enthalpy_coverages;
                size_t i_inter = 0;
                for (size_t i = 0; i < (h_covs.size()-1); i++) {
                    if (h_covs[i] <= cov[j] && cov[j] <= h_covs[i+1]) {
                        i_inter = i;
                        break;
                    }
                }
                double h_left = item.enthalpies[i_inter];
                double h_right = item.enthalpies[i_inter+1];
                m_h_cov[k] += (h_right - h_left) /
                    (h_covs[i_inter+1] - h_covs[i_inter])
                    * (cov[j]  - h_covs[i_inter])
                    + h_left;

                vector_fp s_covs = item.entropy_coverages;
                i_inter = 0;
                for (size_t i = 0; i < (s_covs.size()-1); i++) {
                    if (s_covs[i] <= cov[j] && cov[j] <= s_covs[i+1]) {
                        i_inter = i;
                        break;
                    }
                }
                double s_left = item.entropies[i_inter];
                double s_right = item.entropies[i_inter+1];
                m_s_cov[k] += (s_right - s_left) /
                    (s_covs[i_inter+1] - s_covs[i_inter])
                    * (cov[j] - s_covs[i_inter])
                    + s_left;
            }
        }
    }

    if (m_covstateNumlast != covstateNumnow || m_tlast != tnow || force) {
        // For coverage-depedent heat capacity
        if (m_has_heatcapacity_dependency) {
            for (size_t k = 0; k < m_kk; k++) {
            m_cp_cov[k] = 0.0;
            }
            vector_fp cov(m_kk, 0.0);
            SurfPhase::getCoverages(cov.data());

            for (auto& item : m_HeatCapacityDependency) {
                size_t k = speciesIndex(item.name_k);
                size_t j = speciesIndex(item.name_j);
                double a = item.cpcov_a;
                double b = item.cpcov_b;
                m_cp_cov[k] += (a * log(tnow) + b) * cov[j] * cov[j];
                double int_cp_tnow = tnow * (a * log(tnow) - a + b);
                double int_cp_298 = 298.15 * (a * log(298.15) - a + b);
                m_h_cov[k] += (int_cp_tnow - int_cp_298) * cov[j] * cov[j];
                double int_cp_T_tnow = log(tnow) * (a * log(tnow) + 2 * b);
                double int_cp_T_298 = log(298.15) * (a * log(298.15) + 2 * b);
                m_s_cov[k] += 0.5 * (int_cp_T_tnow - int_cp_T_298) * cov[j] * cov[j];
            }
        }
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_mu_cov[k] = m_h_cov[k] - tnow * m_s_cov[k];
    }
    m_covstateNumlast = covstateNumnow;
    m_tlast = tnow;
}

void CoverageDependentSurfPhase::_updateThermo(bool force) const
{
    _updateReferenceThermo(force);
    _updateCovDepThermo(force);

    for (size_t k = 0; k < m_kk; k++) {
        m_enthalpy[k] = m_h_ref[k] + m_h_cov[k];
        m_entropy[k] = m_s_ref[k] + m_s_cov[k];
        m_heatcapacity[k] = m_cp_ref[k] + m_cp_cov[k];
        m_chempot[k] = m_mu_ref[k] + m_mu_cov[k];
    }
}

}
