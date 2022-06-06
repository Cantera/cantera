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
CoverageDependentSurfPhase::CoverageDependentSurfPhase():
    m_stateNumlast(-2),
    m_theta_ref(1.0)
{
    setSiteDensity(1.0);
    setNDim(2);
}

CoverageDependentSurfPhase::CoverageDependentSurfPhase(const std::string& infile,
                                                       const std::string& id_)
{
    CoverageDependentSurfPhase();
    initThermoFile(infile, id_);
}

void CoverageDependentSurfPhase::setPolynomialDependency(const PolynomialDependency&
                                                         poly_deps)
{
    m_PolynomialDependency.push_back(poly_deps);
}

void CoverageDependentSurfPhase::setPiecewiseDependency(const PiecewiseDependency&
                                                        plin_deps)
{
    double hcov_change = plin_deps.enthalpy_params[2];
    double scov_change = plin_deps.entropy_params[2];

    if (hcov_change <= 0.0 || hcov_change > 1.0 || scov_change <= 0.0 || scov_change > 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::setPiecewiseDependency",
            "Coverage where slope changes must be greater than 0.0 and less than or equal to 1.0.");
    }

    m_PiecewiseDependency.push_back(plin_deps);
}

void CoverageDependentSurfPhase::setInterpolativeDependency(const
                                                            InterpolativeDependency&
                                                            int_deps)
{
    double hcov_last = 0.0;
    for (auto iter=int_deps.enthalpy_map.begin();iter!=int_deps.enthalpy_map.end();++iter){
        if (iter->first < hcov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Coverages are not in ascending order.");
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
    for (auto iter=int_deps.entropy_map.begin();iter!=int_deps.entropy_map.end();++iter){
        if (iter->first < scov_last) {
            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
            "Coverages are not in ascending order.");
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

void CoverageDependentSurfPhase::setHeatCapacityDependency(const HeatCapacityDependency&
                                                           cpcov_deps)
{
    m_HeatCapacityDependency.push_back(cpcov_deps);
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
                   throw CanteraError("CoverageDependentSurfPhase::initThermo",
                        "Unknown species '{}'.", item.first);
                }
                if (j == npos) {
                    throw CanteraError("CoverageDependentSurfPhase::initThermo",
                        "Unknown species '{}'.", item2.first);
                }
                auto& cov_map2 = item2.second.as<AnyMap>();
                // For linear model
                if (cov_map2["model"] == "linear") {
                    vector_fp h_coeffs (5, 0.0);
                    vector_fp s_coeffs (5, 0.0);
                    if (cov_map2.hasKey("enthalpy")) {
                        double h_slope = cov_map2.convert("enthalpy", "J/kmol");
                        h_coeffs[1] = h_slope;
                    }
                    if (cov_map2.hasKey("entropy")) {
                        double s_slope = cov_map2.convert("entropy", "J/kmol/K");
                        s_coeffs[1] = s_slope;
                    }

                    PolynomialDependency poly_deps(k, j, h_coeffs, s_coeffs);
                    setPolynomialDependency(poly_deps);
                // For polynomial(4th) model
                } else if (cov_map2["model"] == "polynomial") {
                    vector_fp h_coeffs (5, 0.0);
                    vector_fp s_coeffs (5, 0.0);
                    if (cov_map2.hasKey("enthalpy-coefficients")) {
                        h_coeffs = cov_map2.convertVector("enthalpy-coefficients", "J/kmol");
                        h_coeffs.insert(h_coeffs.begin(), 0.0);
                    }
                    if (cov_map2.hasKey("entropy-coefficients")) {
                        s_coeffs = cov_map2.convertVector("entropy-coefficients", "J/kmol/K");
                        s_coeffs.insert(s_coeffs.begin(), 0.0);
                    }

                    PolynomialDependency poly_deps(k, j, h_coeffs, s_coeffs);
                    setPolynomialDependency(poly_deps);
                // For piecewise linear model
                } else if (cov_map2["model"] == "piecewise-linear") {
                    vector_fp h_piecewise = {0.0, 0.0, 0.5};
                    vector_fp s_piecewise = {0.0, 0.0, 0.5};
                    if (cov_map2.hasKey("enthalpy-low")) {
                        h_piecewise[0] = cov_map2.convert("enthalpy-low", "J/kmol");
                        h_piecewise[1] = cov_map2.convert("enthalpy-high", "J/kmol");
                        h_piecewise[2] = cov_map2["enthalpy-change"].as<double>();
                    }
                    if (cov_map2.hasKey("entropy-low")) {
                        s_piecewise[0] = cov_map2.convert("entropy-low", "J/kmol/K");
                        s_piecewise[1] = cov_map2.convert("entropy-high", "J/kmol/K");
                        s_piecewise[2] = cov_map2["entropy-change"].as<double>();
                    }

                    PiecewiseDependency plin_deps(k, j, h_piecewise, s_piecewise);
                    setPiecewiseDependency(plin_deps);
                // For interpolative model
                } else if (cov_map2["model"] == "interpolative") {
                    std::map<double, double> hmap, smap;
                    if (cov_map2.hasKey("enthalpy-coverages") && cov_map2.hasKey("enthalpies")) {
                        auto hcovs = cov_map2["enthalpy-coverages"].as<vector_fp>();
                        vector_fp enthalpies = cov_map2.convertVector("enthalpies", "J/kmol");
                        if (hcovs.size() != enthalpies.size()) {
                            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
                            "Sizes of coverages array and enthalpies array are not equal.");
                        }
                        for (size_t i = 0; i < hcovs.size(); i++) {
                            hmap.insert({hcovs[i], enthalpies[i]});
                        }
                    } else {
                        hmap.insert({0.0, 0.0});
                        hmap.insert({1.0, 0.0});
                    }
                    if (cov_map2.hasKey("entropy-coverages") && cov_map2.hasKey("entropies")) {
                        auto scovs = cov_map2["entropy-coverages"].as<vector_fp>();
                        vector_fp entropies = cov_map2.convertVector("entropies", "J/kmol/K");
                        if (scovs.size() != entropies.size()) {
                            throw CanteraError("CoverageDependentSurfPhase::setInterpolativeDependency",
                            "Sizes of coverages array and entropies array are not equal.");
                        }
                        for (size_t i = 0; i < scovs.size(); i++) {
                            smap.insert({scovs[i], entropies[i]});
                        }
                    } else {
                        smap.insert({0.0, 0.0});
                        smap.insert({1.0, 0.0});
                    }

                    InterpolativeDependency int_deps(k, j, hmap, smap);
                    setInterpolativeDependency(int_deps);
                } else {
                    throw CanteraError("CoverageDependentSurfPhase::initThermo",
                        "Unrecognized coverage-dependencies model between '{}' and '{}'.",
                        item.first, item2.first);
                }
                // For coverage-dependent heat capacity parameters, if present
                if (cov_map2.hasKey("heat-capacity-a")) {
                    double cpcov_a = cov_map2.convert("heat-capacity-a", "J/kmol/K");
                    double cpcov_b = cov_map2.convert("heat-capacity-b", "J/kmol/K");

                    HeatCapacityDependency cpcov_deps(k, j, cpcov_a, cpcov_b);
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

// Functions calculating reference state thermodynamic properties--------------

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

// Functions calculating standard state thermodynamic properties---------------

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

// Functions calling partial molar thermodynamic properties----------------

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

// Functions calculating mixture thermodynamic properties--------------------------

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

void CoverageDependentSurfPhase::_updateCovDepThermo(bool force) const
{
    int stateNumnow = stateMFNumber();
    double tnow = temperature();
    if (m_stateNumlast != stateNumnow || m_tlast != tnow || force) {
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

        // For piecewise linear model
        for (auto& item : m_PiecewiseDependency) {
            double h_slope_low = item.enthalpy_params[0];
            double h_slope_high = item.enthalpy_params[1];
            double h_cov_change = item.enthalpy_params[2];
            if (m_cov[item.j] <= h_cov_change) {
                m_h_cov[item.k] += h_slope_low * m_cov[item.j];
            } else {
                m_h_cov[item.k] += h_slope_low * h_cov_change;
                m_h_cov[item.k] += h_slope_high * (m_cov[item.j] - h_cov_change);
            }

            double s_slope_low = item.entropy_params[0];
            double s_slope_high = item.entropy_params[1];
            double s_cov_change = item.entropy_params[2];
            if (m_cov[item.j] <= s_cov_change) {
                m_s_cov[item.k] += s_slope_low * m_cov[item.j];
            } else {
                m_s_cov[item.k] += s_slope_low * s_cov_change;
                m_s_cov[item.k] += s_slope_high * (m_cov[item.j] - s_cov_change);
            }
        }

        // For interpolative model
        for (auto& item : m_InterpolativeDependency) {
            auto h_iter = item.enthalpy_map.upper_bound(m_cov[item.j]);
            auto s_iter = item.entropy_map.upper_bound(m_cov[item.j]);
            AssertThrowMsg(h_iter != m_cov.end(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[iterm.j]);
            AssertThrowMsg(h_iter != m_cov.begin(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[iterm.j]);
            AssertThrowMsg(s_iter != m_cov.end(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[iterm.j]);
            AssertThrowMsg(s_iter != m_cov.begin(),
                           "CoverageDependentSurfPhase::_updateCovDepThermo",
                           "Coverage out of range: {}", m_cov[iterm.j]);

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
            m_h_cov[item.k] += (int_cp_tnow - int_cp_298) * m_cov[item.j] * m_cov[item.j];
            double int_cp_T_tnow = log(tnow) * (a * log(tnow) + 2 * b);
            double int_cp_T_298 = log(298.15) * (a * log(298.15) + 2 * b);
            m_s_cov[item.k] += 0.5 * (int_cp_T_tnow - int_cp_T_298) * m_cov[item.j] * m_cov[item.j];
        }

        for (size_t k = 0; k < m_kk; k++) {
            m_mu_cov[k] = m_h_cov[k] - tnow * m_s_cov[k];
        }
        m_stateNumlast = stateNumnow;
    }
}

void CoverageDependentSurfPhase::_updateTotalThermo(bool force) const
{
    _updateCovDepThermo(force);
    SurfPhase::_updateThermo(force);

    for (size_t k = 0; k < m_kk; k++) {
        m_enthalpy[k] = m_h0[k] + m_h_cov[k];
        m_entropy[k] = m_s0[k] + m_s_cov[k];
        m_heatcapacity[k] = m_cp0[k] + m_cp_cov[k];
        m_chempot[k] = m_mu0[k] + m_mu_cov[k];
    }
}

}
