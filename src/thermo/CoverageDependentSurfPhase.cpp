/**
 *  @file CoverageDependentSurfPhase.cpp
 *  Definitions for a thermodynamics model of a coverage-dependent surface
 *  phase derived from SurfPhase, applying adsorbate lateral interaction
 *  correction factors to the SurfPhase thermodynamic properties.
 *  (see @ref thermoprops and class
 *  @link Cantera::CoverageDependentSurfPhase CoverageDependentSurfPhase@endlink).
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

CoverageDependentSurfPhase::PolynomialDependency::PolynomialDependency(
    size_t k, size_t j, const AnyMap& dep_map
) :
    k(k),
    j(j),
    enthalpy_coeffs({0.0, 0.0, 0.0, 0.0, 0.0}),
    entropy_coeffs({0.0, 0.0, 0.0, 0.0, 0.0}),
    isLinear(false)
{
    // For linear model
    if (dep_map["model"] == "linear") {
        if (dep_map.hasKey("enthalpy")) {
            enthalpy_coeffs[1] = dep_map.convert("enthalpy", "J/kmol");
        }
        if (dep_map.hasKey("entropy")) {
            entropy_coeffs[1] = dep_map.convert("entropy", "J/kmol/K");
        }
        isLinear = true;
    // For polynomial(4th) model
    } else if (dep_map["model"] == "polynomial") {
        if (dep_map.hasKey("enthalpy-coefficients")) {
            enthalpy_coeffs = dep_map.convertVector(
                "enthalpy-coefficients", "J/kmol");
            enthalpy_coeffs.insert(enthalpy_coeffs.begin(), 0.0);
        }
        if (dep_map.hasKey("entropy-coefficients")) {
            entropy_coeffs = dep_map.convertVector(
                "entropy-coefficients", "J/kmol/K");
            entropy_coeffs.insert(entropy_coeffs.begin(), 0.0);
        }
    }
}

CoverageDependentSurfPhase::InterpolativeDependency::InterpolativeDependency(
    size_t k, size_t j, const AnyMap& dep_map, const AnyBase& node
) :
    k(k),
    j(j),
    enthalpy_map({{0.0, 0.0}, {1.0, 0.0}}),
    entropy_map({{0.0, 0.0}, {1.0, 0.0}}),
    isPiecewise(false)
{
    // For piecewise-linear model
    // Piecewise-linear model coefficients are converted into
    // a map <coverages: values>
    if (dep_map["model"] == "piecewise-linear") {
        if (dep_map.hasKey("enthalpy-low") ||
            dep_map.hasKey("enthalpy-change") ||
            dep_map.hasKey("enthalpy-high"))
        {
            auto cov_change = dep_map["enthalpy-change"].as<double>();
            enthalpy_map[cov_change] =
                dep_map.convert("enthalpy-low", "J/kmol") * cov_change;
            enthalpy_map[1.0] = (1.0 - cov_change)
                * dep_map.convert("enthalpy-high", "J/kmol")
                + enthalpy_map[cov_change];
        }
        if (dep_map.hasKey("entropy-low") ||
            dep_map.hasKey("entropy-change") ||
            dep_map.hasKey("entropy-high"))
        {
            auto cov_change = dep_map["entropy-change"].as<double>();
            entropy_map[cov_change] =
                dep_map.convert("entropy-low", "J/kmol/K") * cov_change;
            entropy_map[1.0] = (1.0 - cov_change)
                * dep_map.convert("entropy-high", "J/kmol/K")
                + entropy_map[cov_change];
        }
        isPiecewise = true;
    // For interpolative model
    } else if (dep_map["model"] == "interpolative") {
        if (dep_map.hasKey("enthalpy-coverages") || dep_map.hasKey("enthalpies")) {
            auto hcovs = dep_map["enthalpy-coverages"].as<vector<double>>();
            vector<double> enthalpies = dep_map.convertVector("enthalpies", "J/kmol");
            if (hcovs.size() != enthalpies.size()) {
                throw InputFileError("CoverageDependentSurfPhase::"
                    "addInterpolativeDependency", node,
                    "Sizes of coverages array and enthalpies array are not equal.");
            }
            for (size_t i = 0; i < hcovs.size(); i++) {
                enthalpy_map[hcovs[i]] = enthalpies[i];
            }
        }
        if (dep_map.hasKey("entropy-coverages") || dep_map.hasKey("entropies")) {
            auto scovs = dep_map["entropy-coverages"].as<vector<double>>();
            vector<double> entropies = dep_map.convertVector("entropies",
                                                        "J/kmol/K");
            if (scovs.size() != entropies.size()) {
                throw InputFileError("CoverageDependentSurfPhase::"
                    "addInterpolativeDependency", node,
                    "Sizes of coverages array and entropies array are not equal.");
            }
            for (size_t i = 0; i < scovs.size(); i++) {
                entropy_map[scovs[i]] = entropies[i];
            }
        }
    }
}

CoverageDependentSurfPhase::CoverageDependentSurfPhase(const string& infile,
                                                       const string& id_):
    m_theta_ref(1.0),
    m_stateNumlast(-2)
{
    setNDim(2);
    initThermoFile(infile, id_);
}

void CoverageDependentSurfPhase::addInterpolativeDependency(
    const InterpolativeDependency& int_deps)
{
    if (int_deps.enthalpy_map.begin()->first != 0.0) {
        throw CanteraError("CoverageDependentSurfPhase::addInterpolativeDependency",
            "The first element of enthalpy-coverages array must be 0.0.");
    }
    if (int_deps.enthalpy_map.rbegin()->first != 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::addInterpolativeDependency",
            "The last element of enthalpy-coverages array must be 1.0.");
    }

    if (int_deps.entropy_map.begin()->first != 0.0) {
        throw CanteraError("CoverageDependentSurfPhase::addInterpolativeDependency",
            "The first element of entropy-coverages array must be 0.0.");
    }
    if (int_deps.entropy_map.rbegin()->first != 1.0) {
        throw CanteraError("CoverageDependentSurfPhase::addInterpolativeDependency",
            "The last element of entropy-coverages array must be 1.0.");
    }

    m_InterpolativeDependency.push_back(int_deps);
}

void CoverageDependentSurfPhase::initThermo()
{
    SurfPhase::initThermo();
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
                size_t k = speciesIndex(item.first, false);
                size_t j = speciesIndex(item2.first, false);
                if (k == npos) {
                   throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unknown species '{}'.", item.first);
                }
                if (j == npos) {
                    throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unknown species '{}'.", item2.first);
                }
                auto& dep_map = item2.second.as<AnyMap>();
                // For linear model and polynomial model
                if (dep_map["model"] == "linear" || dep_map["model"] == "polynomial") {
                    PolynomialDependency poly_deps(k, j, dep_map);
                    m_PolynomialDependency.push_back(poly_deps);
                // For piecewise-linear model and interpolative model
                } else if (dep_map["model"] == "piecewise-linear" ||
                    dep_map["model"] == "interpolative") {
                    InterpolativeDependency int_deps(k, j, dep_map,
                                                     item.second->input);
                    addInterpolativeDependency(int_deps);
                } else {
                    throw InputFileError("CoverageDependentSurfPhase::initThermo",
                        item.second->input, "Unrecognized coverage dependency model. \
                        Model must be 'linear', 'piecewise-linear', 'polynomial', \
                        or 'interpolative'.");
                }
                // For coverage-dependent heat capacity parameters, if present
                if (dep_map.hasKey("heat-capacity-a")) {
                    HeatCapacityDependency cpcov_deps(k, j);
                    cpcov_deps.coeff_a = dep_map.convert("heat-capacity-a", "J/kmol/K");
                    cpcov_deps.coeff_b = dep_map.convert("heat-capacity-b", "J/kmol/K");

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

void CoverageDependentSurfPhase::getParameters(AnyMap& phaseNode) const
{
    SurfPhase::getParameters(phaseNode);
    phaseNode["reference-state-coverage"] = m_theta_ref;
}

void CoverageDependentSurfPhase::getSpeciesParameters(const string& name,
                                                      AnyMap& speciesNode) const
{
    SurfPhase::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name, true);
    // Get linear and polynomial model parameters from PolynomialDependency vector
    for (auto& item : m_PolynomialDependency) {
        if (item.k == k) {
            if (item.isLinear) {
                auto& covdepNode =
                    speciesNode["coverage-dependencies"][speciesName(item.j)]
                        .getMapWhere("model", "linear", true);
                covdepNode["enthalpy"].setQuantity(item.enthalpy_coeffs[1], "J/kmol");
                covdepNode["entropy"].setQuantity(item.entropy_coeffs[1], "J/kmol/K");
            } else {
                auto& covdepNode =
                    speciesNode["coverage-dependencies"][speciesName(item.j)]
                        .getMapWhere("model", "polynomial", true);
                vector<double> hvec (
                    item.enthalpy_coeffs.begin() + 1, item.enthalpy_coeffs.end());
                covdepNode["enthalpy-coefficients"].setQuantity(hvec, "J/kmol");
                vector<double> svec (
                    item.entropy_coeffs.begin() + 1, item.entropy_coeffs.end());
                covdepNode["entropy-coefficients"].setQuantity(svec, "J/kmol/K");
            }
        }
    }
    // Get piecewise-linear model parameters from InterpolativeDependency vector
    for (auto& item : m_InterpolativeDependency) {
        if (item.k == k) {
            if (item.isPiecewise) {
                auto& covdepNode =
                    speciesNode["coverage-dependencies"][speciesName(item.j)]
                        .getMapWhere("model", "piecewise-linear", true);
                vector<double> hcovs, enthalpies, scovs, entropies;
                for (const auto& hmap : item.enthalpy_map) {
                    hcovs.push_back(hmap.first);
                    enthalpies.push_back(hmap.second);
                }
                for (const auto& smap : item.entropy_map) {
                    scovs.push_back(smap.first);
                    entropies.push_back(smap.second);
                }
                covdepNode["enthalpy-change"] = hcovs[1];
                covdepNode["enthalpy-low"].setQuantity(
                    (enthalpies[1] - enthalpies[0]) / (hcovs[1] - hcovs[0]), "J/kmol");
                covdepNode["enthalpy-high"].setQuantity(
                    (enthalpies[2] - enthalpies[1]) / (hcovs[2] - hcovs[1]), "J/kmol");
                covdepNode["entropy-change"] = scovs[1];
                covdepNode["entropy-low"].setQuantity(
                    (entropies[1] - entropies[0]) / (scovs[1] - scovs[0]), "J/kmol/K");
                covdepNode["entropy-high"].setQuantity(
                    (entropies[2] - entropies[1]) / (scovs[2] - scovs[1]), "J/kmol/K");
            } else {
                auto& covdepNode =
                    speciesNode["coverage-dependencies"][speciesName(item.j)]
                        .getMapWhere("model", "interpolative", true);
                vector<double> hcovs, enthalpies, scovs, entropies;
                for (const auto& hmap : item.enthalpy_map) {
                    hcovs.push_back(hmap.first);
                    enthalpies.push_back(hmap.second);
                }
                for (const auto& smap : item.entropy_map) {
                    scovs.push_back(smap.first);
                    entropies.push_back(smap.second);
                }
                covdepNode["enthalpy-coverages"] = std::move(hcovs);
                covdepNode["enthalpies"].setQuantity(enthalpies, "J/kmol");
                covdepNode["entropy-coverages"] = std::move(scovs);
                covdepNode["entropies"].setQuantity(entropies, "J/kmol/K");
            }
        }
    }
    // Get heat capacity model parameters from HeatCapacityDependency vector
    for (auto& item : m_HeatCapacityDependency) {
        if (item.k == k) {
            auto& covdepNode =
                speciesNode["coverage-dependencies"][speciesName(item.j)]
                    .getMapWhere("heat-capacity-a", "", true);
            covdepNode["heat-capacity-a"].setQuantity(item.coeff_a, "J/kmol/K");
            covdepNode["heat-capacity-b"].setQuantity(item.coeff_b, "J/kmol/K");
        }
    }
}

void CoverageDependentSurfPhase::getGibbs_RT_ref(span<double> grt) const
{
    checkArraySize("CoverageDependentSurfPhase::getGibbs_RT_ref", grt.size(), m_kk);
    SurfPhase::_updateThermo();
    scale(m_mu0.begin(), m_mu0.end(), grt.begin(), 1.0/RT());
}

void CoverageDependentSurfPhase::getEnthalpy_RT_ref(span<double> hrt) const
{
    checkArraySize("CoverageDependentSurfPhase::getEnthalpy_RT_ref", hrt.size(), m_kk);
    SurfPhase::_updateThermo();
    scale(m_h0.begin(), m_h0.end(), hrt.begin(), 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R_ref(span<double> sr) const
{
    checkArraySize("CoverageDependentSurfPhase::getEntropy_R_ref", sr.size(), m_kk);
    SurfPhase::_updateThermo();
    scale(m_s0.begin(), m_s0.end(), sr.begin(), 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getCp_R_ref(span<double> cpr) const
{
    checkArraySize("CoverageDependentSurfPhase::getCp_R_ref", cpr.size(), m_kk);
    SurfPhase::_updateThermo();
    scale(m_cp0.begin(), m_cp0.end(), cpr.begin(), 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getEnthalpy_RT(span<double> hrt) const
{
    checkArraySize("CoverageDependentSurfPhase::getEnthalpy_RT", hrt.size(), m_kk);
    _updateTotalThermo();
    scale(m_enthalpy.begin(), m_enthalpy.end(), hrt.begin(), 1.0/RT());
}

void CoverageDependentSurfPhase::getEntropy_R(span<double> sr) const
{
    checkArraySize("CoverageDependentSurfPhase::getEntropy_R", sr.size(), m_kk);
    _updateTotalThermo();
    scale(m_entropy.begin(), m_entropy.end(), sr.begin(), 1.0/GasConstant);
    if (m_theta_ref != 1.0) {
        double tmp = -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            sr[k] -= tmp;
        }
    }
}

void CoverageDependentSurfPhase::getCp_R(span<double> cpr) const
{
    checkArraySize("CoverageDependentSurfPhase::getCp_R", cpr.size(), m_kk);
    _updateTotalThermo();
    scale(m_heatcapacity.begin(), m_heatcapacity.end(), cpr.begin(), 1.0/GasConstant);
}

void CoverageDependentSurfPhase::getGibbs_RT(span<double> grt) const
{
    checkArraySize("CoverageDependentSurfPhase::getGibbs_RT", grt.size(), m_kk);
    _updateTotalThermo();
    scale(m_chempot.begin(), m_chempot.end(), grt.begin(), 1.0/RT());
    if (m_theta_ref != 1.0) {
        double tmp = -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            grt[k] += tmp;
        }
    }
}

void CoverageDependentSurfPhase::getStandardChemPotentials(span<double> mu0) const
{
    checkArraySize("CoverageDependentSurfPhase::getStandardChemPotentials",
                   mu0.size(), m_kk);
    _updateTotalThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu0.begin());
    if (m_theta_ref != 1.0) {
        double tmp = RT() * -log(m_theta_ref);
        for (size_t k = 0; k < m_kk; k++) {
            mu0[k] += tmp;
        }
    }
}

void CoverageDependentSurfPhase::getPartialMolarEnthalpies(span<double> hbar) const
{
    checkArraySize("CoverageDependentSurfPhase::getPartialMolarEnthalpies",
                   hbar.size(), m_kk);
    _updateTotalThermo();
    copy(m_enthalpy.begin(), m_enthalpy.end(), hbar.begin());
}

void CoverageDependentSurfPhase::getPartialMolarEntropies(span<double> sbar) const
{
    checkArraySize("CoverageDependentSurfPhase::getPartialMolarEntropies",
                   sbar.size(), m_kk);
    _updateTotalThermo();
    copy(m_entropy.begin(), m_entropy.end(), sbar.begin());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= GasConstant * log(std::max(m_cov[k], SmallNumber) / m_theta_ref);
    }
}

void CoverageDependentSurfPhase::getPartialMolarCp(span<double> cpbar) const
{
    checkArraySize("CoverageDependentSurfPhase::getPartialMolarCp", cpbar.size(), m_kk);
    _updateTotalThermo();
    copy(m_heatcapacity.begin(), m_heatcapacity.end(), cpbar.begin());
}

void CoverageDependentSurfPhase::getChemPotentials(span<double> mu) const
{
    checkArraySize("CoverageDependentSurfPhase::getChemPotentials", mu.size(), m_kk);
    _updateTotalThermo();
    copy(m_chempot.begin(), m_chempot.end(), mu.begin());
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
        getCoverages(m_cov);

        // For linear and polynomial model
        for (auto& item : m_PolynomialDependency) {
            m_h_cov[item.k] += poly4(m_cov[item.j], item.enthalpy_coeffs);
            m_s_cov[item.k] += poly4(m_cov[item.j], item.entropy_coeffs);
        }

        // For piecewise-linear and interpolative model
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
            double a = item.coeff_a;
            double b = item.coeff_b;
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
