/**
 *  @file CoverageDependentSurfPhase.cpp
 *  Definitions for a thermodynamics model of a coverage-dependent surface
 *  phase derived from SurfPhase, applying adsorbate lateral interaction
 *  correction factors to the SurfPhase thermodynamic properties.
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
#include "cantera/thermo/Species.h"

using namespace std;

namespace Cantera
{

CoverageDependentSurfPhase::CoverageDependentSurfPhase(const std::string& infile,
                                                       const std::string& id_):
    m_stateNumlast(-2),
    m_theta_ref(1.0)
{
    setNDim(2);
    initThermoFile(infile, id_);
}

void CoverageDependentSurfPhase::setInterpolativeDependency(const
                                                            InterpolativeDependency&
                                                            int_deps)
{
    double hcov_last = 0.0;
    for (auto iter=int_deps.enthalpy_map.begin();iter!=int_deps.enthalpy_map.end();
        ++iter){
        if (iter->first < hcov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency"
            , "Coverages are not in ascending order.");
        }
        hcov_last = iter->first;
    }
    if (int_deps.enthalpy_map.count(0.0) == 0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The first element of enthalpy-coverages array must be 0.0.");
    }
    if (int_deps.enthalpy_map.count(1.0) == 0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The last element of enthalpy-coverages array must be 1.0.");
    }

    double scov_last = 0.0;
    for (auto iter=int_deps.entropy_map.begin();iter!=int_deps.entropy_map.end();
        ++iter){
        if (iter->first < scov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency"
            , "Coverages are not in ascending order.");
        }
        scov_last = iter->first;
    }
    if (int_deps.entropy_map.count(0.0) == 0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The first element of entropy-coverages array must be 0.0.");
    }
    if (int_deps.entropy_map.count(1.0) == 0) {
        throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "The last element of entropy-coverages array must be 1.0.");
    }

    m_InterpolativeDependency.push_back(int_deps);
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
            throw InputFileError("CoverageDependentSurfPhase::initThermo",
                m_input, "Reference state coverage must be greater than 0.0 and \
                less than or equal to 1.0.");
        }
    }
    for (auto& item : m_species) {
        // Read enthalpy and entropy dependencies from species 'input' information
        // (i.e. as specified in a YAML input file) for both self- and cross-
        // interactions.
        if (item.second->input.hasKey("coverage-dependencies")) {
            auto& cov_map = item.second->input["coverage-dependencies"];
            for (const auto& item2 : cov_map) {
                size_t k = speciesIndex(item.first);
                size_t j = speciesIndex(item2.first);
                if (k == npos) {
                   throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unknown species '{}'.", item.first);
                }
                if (j == npos) {
                    throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unknown species '{}'.", item2.first);
                }
                auto& cov_map2 = item2.second.as<AnyMap>();
                // For linear model
                if (cov_map2["model"] == "linear") {
                    PolynomialDependency poly_deps(k, j);
                    if (cov_map2.hasKey("enthalpy")) {
                        poly_deps.enthalpy_coeffs[1] = cov_map2.convert("enthalpy",
                                                                        "J/kmol");
                    }
                    if (cov_map2.hasKey("entropy")) {
                        poly_deps.entropy_coeffs[1] = cov_map2.convert("entropy",
                                                                       "J/kmol/K");
                    }

                    m_PolynomialDependency.push_back(poly_deps);
                // For polynomial(4th) model
                } else if (cov_map2["model"] == "polynomial") {
                    PolynomialDependency poly_deps(k, j);
                    if (cov_map2.hasKey("enthalpy-coefficients")) {
                        poly_deps.enthalpy_coeffs = cov_map2.convertVector(
                            "enthalpy-coefficients", "J/kmol");
                        poly_deps.enthalpy_coeffs.insert(
                            poly_deps.enthalpy_coeffs.begin(), 0.0);
                    }
                    if (cov_map2.hasKey("entropy-coefficients")) {
                        poly_deps.entropy_coeffs = cov_map2.convertVector(
                            "entropy-coefficients", "J/kmol/K");
                        poly_deps.entropy_coeffs.insert(
                            poly_deps.entropy_coeffs.begin(), 0.0);
                    }

                    m_PolynomialDependency.push_back(poly_deps);
                // For piecewise linear model
                } else if (cov_map2["model"] == "piecewise-linear") {
                    InterpolativeDependency int_deps(k, j);
                    if (cov_map2.hasKey("enthalpy-low") ||
                        cov_map2.hasKey("enthalpy-change") ||
                        cov_map2.hasKey("enthalpy-high")) {
                        auto cov_change = cov_map2["enthalpy-change"].as<double>();
                        int_deps.enthalpy_map[cov_change] =
                            cov_map2.convert("enthalpy-low", "J/kmol") * cov_change;
                        int_deps.enthalpy_map[1.0] = (1.0 - cov_change) *
                            cov_map2.convert("enthalpy-high", "J/kmol")
                            + int_deps.enthalpy_map[cov_change];
                    }
                    if (cov_map2.hasKey("entropy-low") ||
                        cov_map2.hasKey("entropy-change") ||
                        cov_map2.hasKey("entropy-high")) {
                        auto cov_change = cov_map2["entropy-change"].as<double>();
                        int_deps.entropy_map[cov_change] =
                            cov_map2.convert("entropy-low", "J/kmol/K") * cov_change;
                        int_deps.entropy_map[1.0] = (1.0 - cov_change) *
                            cov_map2.convert("entropy-high", "J/kmol/K")
                            + int_deps.entropy_map[cov_change];
                    }

                    setInterpolativeDependency(int_deps);
                // For interpolative model
                } else if (cov_map2["model"] == "interpolative") {
                    InterpolativeDependency int_deps(k, j);
                    if (cov_map2.hasKey("enthalpy-coverages") ||
                        cov_map2.hasKey("enthalpies")) {
                        auto hcovs = cov_map2["enthalpy-coverages"].as<vector_fp>();
                        vector_fp enthalpies = cov_map2.convertVector("enthalpies",
                                                                      "J/kmol");
                        if (hcovs.size() != enthalpies.size()) {
                            throw InputFileError("CoverageDependentSurfPhase::\
                            setInterpolativeDependency", item.second->input,
                            "Sizes of coverages array and enthalpies array are \
                            not equal.");
                        }
                        for (size_t i = 0; i < hcovs.size(); i++) {
                            int_deps.enthalpy_map[hcovs[i]] = enthalpies[i];
                        }
                    }
                    if (cov_map2.hasKey("entropy-coverages") ||
                        cov_map2.hasKey("entropies")) {
                        auto scovs = cov_map2["entropy-coverages"].as<vector_fp>();
                        vector_fp entropies = cov_map2.convertVector("entropies",
                                                                     "J/kmol/K");
                        if (scovs.size() != entropies.size()) {
                            throw InputFileError("CoverageDependentSurfPhase::\
                            setInterpolativeDependency", item.second->input,
                            "Sizes of coverages array and entropies array are \
                            not equal.");
                        }
                        for (size_t i = 0; i < scovs.size(); i++) {
                            int_deps.entropy_map[scovs[i]] = entropies[i];
                        }
                    }

                    setInterpolativeDependency(int_deps);
                } else {
                    throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unrecognized coverage dependency model. \
                        Model must be 'linear', 'piecewise-linear', 'polynomial', \
                        or 'interpolative'.");
                }
                // For coverage-dependent heat capacity parameters, if present
                if (cov_map2.hasKey("heat-capacity-a")) {
                    HeatCapacityDependency cpcov_deps(k, j);
                    cpcov_deps.cpcov_a = cov_map2.convert("heat-capacity-a", "J/kmol/K");
                    cpcov_deps.cpcov_b = cov_map2.convert("heat-capacity-b", "J/kmol/K");

                    m_HeatCapacityDependency.push_back(cpcov_deps);
                }
            }
        }
    }
}

bool CoverageDependentSurfPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = SurfPhase::addSpecies(spec);
    if (added) {
        m_cov.push_back(0.0);
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

void CoverageDependentSurfPhase::getGibbs_RT_ref(double* grt) const
{
    SurfPhase::_updateThermo();
    scale(m_mu0.begin(), m_mu0.end(), grt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEnthalpy_RT_ref(double* hrt) const
{
    SurfPhase::_updateThermo();
    scale(m_h0.begin(), m_h0.end(), hrt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R_ref(double* sr) const
{
    SurfPhase::_updateThermo();
    scale(m_s0.begin(), m_s0.end(), sr, 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getCp_R_ref(double* cpr) const
{
    SurfPhase::_updateThermo();
    scale(m_cp0.begin(), m_cp0.end(), cpr, 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getEnthalpy_RT(double* hrt) const
{
    _updateTotalThermo();
    scale(m_enthalpy.begin(), m_enthalpy.end(), hrt, 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R(double* sr) const
{
    _updateTotalThermo();
    scale(m_entropy.begin(), m_entropy.end(), sr, 1.0/GasConstant);
    if (m_theta_ref != 1.0) {
        double tmp = -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            sr[k] -= tmp;
        }
    }
}

void CoverageDependentSurfPhase::getCp_R(double* cpr) const
{
    _updateTotalThermo();
    scale(m_heatcapacity.begin(), m_heatcapacity.end(), cpr, 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getGibbs_RT(double* grt) const
{
    _updateTotalThermo();
    scale(m_chempot.begin(), m_chempot.end(), grt, 1.0/RT());
    if (m_theta_ref != 1.0) {
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
    _updateTotalThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu0);
    if (m_theta_ref != 1.0) {
        double tmp = RT() * -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            mu0[k] += tmp;
        }
    }
}

void CoverageDependentSurfPhase::getPartialMolarEnthalpies(double* hbar) const
{
    _updateTotalThermo();
    copy(m_enthalpy.begin(), m_enthalpy.end(), hbar);
}

void CoverageDependentSurfPhase::getPartialMolarEntropies(double* sbar) const
{
    _updateTotalThermo();
    copy(m_entropy.begin(), m_entropy.end(), sbar);
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= GasConstant * log(std::max(m_cov[k], SmallNumber)
            / m_theta_ref);
    }
}

void CoverageDependentSurfPhase::getPartialMolarCp(double* cpbar) const
{
    _updateTotalThermo();
    copy(m_heatcapacity.begin(), m_heatcapacity.end(), cpbar);
}

void CoverageDependentSurfPhase::getChemPotentials(double* mu) const
{
    _updateTotalThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu);
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += RT() * log(std::max(m_cov[k], SmallNumber) / m_theta_ref);
    }
}

double CoverageDependentSurfPhase::enthalpy_mole() const
{
    _updateTotalThermo();
    return mean_X(m_enthalpy);
}

double CoverageDependentSurfPhase::entropy_mole() const
{
    _updateTotalThermo();
    double entropy = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        entropy += moleFraction(k) * (m_entropy[k] -
            GasConstant * log(std::max(m_cov[k], SmallNumber) / m_theta_ref));
    }
    return entropy;
}

double CoverageDependentSurfPhase::cp_mole() const
{
    _updateTotalThermo();
    return mean_X(m_heatcapacity);
}

void CoverageDependentSurfPhase::_updateCovDepThermo() const
{
    int stateNumnow = stateMFNumber();
    double tnow = temperature();
    if (m_stateNumlast != stateNumnow || m_tlast != tnow) {
        for (size_t k = 0; k < m_kk; k++) {
            m_h_cov[k] = 0.0;
            m_s_cov[k] = 0.0;
            m_cp_cov[k] = 0.0;
        }
        getCoverages(m_cov.data());

        // For linear and polynomial model
        for (auto& item : m_PolynomialDependency) {
            m_h_cov[item.k] += poly4(m_cov[item.j], item.enthalpy_coeffs.data());
            m_s_cov[item.k] += poly4(m_cov[item.j], item.entropy_coeffs.data());
        }

        // For piecewise linear and interpolative model
        for (auto& item : m_InterpolativeDependency) {
            auto h_iter = item.enthalpy_map.upper_bound(m_cov[item.j]);
            auto s_iter = item.entropy_map.upper_bound(m_cov[item.j]);
            AssertThrowMsg(h_iter != item.enthalpy_map.end(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[item.j]);
            AssertThrowMsg(h_iter != item.enthalpy_map.begin(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[item.j]);
            AssertThrowMsg(s_iter != item.entropy_map.end(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[item.j]);
            AssertThrowMsg(s_iter != item.entropy_map.begin(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[item.j]);

            double highHcov = h_iter->first;
            double highH = h_iter->second;
            double lowHcov = (--h_iter)->first;
            double lowH = h_iter->second;

            double highScov = s_iter->first;
            double highS = s_iter->second;
            double lowScov = (--s_iter)->first;
            double lowS = s_iter->second;

            m_h_cov[item.k] += (highH - lowH) / (highHcov - lowHcov)
                * (m_cov[item.j] - lowHcov) + lowH;

            m_s_cov[item.k] += (highS - lowS) / (highScov - lowScov)
                * (m_cov[item.j] - lowScov) + lowS;
        }

        // For coverage-dependent heat capacity
        for (auto& item : m_HeatCapacityDependency) {
            double a = item.cpcov_a;
            double b = item.cpcov_b;
            m_cp_cov[item.k] += (a * log(tnow) + b) * m_cov[item.j] * m_cov[item.j];
            double int_cp_tnow = tnow * (a * log(tnow) - a + b);
            double int_cp_298 = 298.15 * (a * log(298.15) - a + b);
            m_h_cov[item.k] += (int_cp_tnow - int_cp_298) * m_cov[item.j]
                               * m_cov[item.j];
            double int_cp_T_tnow = log(tnow) * (a * log(tnow) + 2 * b);
            double int_cp_T_298 = log(298.15) * (a * log(298.15) + 2 * b);
            m_s_cov[item.k] += 0.5 * (int_cp_T_tnow - int_cp_T_298) * m_cov[item.j]
                               * m_cov[item.j];
        }

        for (size_t k = 0; k < m_kk; k++) {
            m_mu_cov[k] = m_h_cov[k] - tnow * m_s_cov[k];
        }
        m_stateNumlast = stateNumnow;
    }
}

void CoverageDependentSurfPhase::_updateTotalThermo() const
{
    _updateCovDepThermo();
    SurfPhase::_updateThermo();

    for (size_t k = 0; k < m_kk; k++) {
        m_enthalpy[k] = m_h0[k] + m_h_cov[k];
        m_entropy[k] = m_s0[k] + m_s_cov[k];
        m_heatcapacity[k] = m_cp0[k] + m_cp_cov[k];
        m_chempot[k] = m_mu0[k] + m_mu_cov[k];
    }
}

}
