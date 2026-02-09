//! @file InterfaceRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/InterfaceRate.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void InterfaceData::update(double T)
{
    throw CanteraError("InterfaceData::update",
        "Missing state information: 'InterfaceData' requires species coverages.");
}

void InterfaceData::update(double T, span<const double> values)
{
    warn_user("InterfaceData::update",
        "This method does not update the site density.");
    ReactionData::update(T);
    sqrtT = sqrt(T);
    if (coverages.size() == 0) {
        coverages.assign(values.begin(), values.end());
        logCoverages.resize(values.size());
    } else if (values.size() == coverages.size()) {
        std::copy(values.begin(), values.end(), coverages.begin());
    } else {
        throw CanteraError("InterfaceData::update",
            "Incompatible lengths of coverage arrays: received {} elements while "
            "{} are required.", values.size(), coverages.size());
    }
    for (size_t n = 0; n < coverages.size(); n++) {
        logCoverages[n] = std::log(std::max(coverages[n], Tiny));
    }
}

bool InterfaceData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    int mf = 0;
    for (size_t n = 0; n < kin.nPhases(); n++) {
        mf += kin.thermo(n).stateMFNumber();
    }

    double T = phase.temperature();
    bool changed = false;
    const auto& surf = dynamic_cast<const SurfPhase&>(kin.thermo(0));
    double site_density = surf.siteDensity();
    if (density != site_density) {
        density = surf.siteDensity();
        changed = true;
    }
    if (T != temperature) {
        ReactionData::update(T);
        sqrtT = sqrt(T);
        changed = true;
    }
    if (changed || mf != m_state_mf_number) {
        surf.getCoverages(coverages);
        for (size_t n = 0; n < coverages.size(); n++) {
            logCoverages[n] = std::log(std::max(coverages[n], Tiny));
        }
        for (size_t n = 0; n < kin.nPhases(); n++) {
            size_t start = kin.kineticsSpeciesIndex(0, n);
            const auto& ph = kin.thermo(n);
            electricPotentials[n] = ph.electricPotential();
            ph.getPartialMolarEnthalpies(
                span<double>(partialMolarEnthalpies).subspan(start, ph.nSpecies()));
            ph.getStandardChemPotentials(
                span<double>(standardChemPotentials).subspan(start, ph.nSpecies()));
            size_t nsp = ph.nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                // only used for exchange current density formulation
                standardConcentrations[k + start] = ph.standardConcentration(k);
            }
        }
        m_state_mf_number = mf;
        changed = true;
    }
    return changed;
}

void InterfaceData::perturbTemperature(double deltaT)
{
    throw NotImplementedError("InterfaceData::perturbTemperature");
}

void InterfaceData::resize(Kinetics& kin) {
    coverages.resize(kin.thermo().nSpecies(), 0.);
    logCoverages.resize(kin.thermo().nSpecies(), 0.);
    partialMolarEnthalpies.resize(kin.nTotalSpecies(), 0.);
    electricPotentials.resize(kin.nPhases(), 0.);
    standardChemPotentials.resize(kin.nTotalSpecies(), 0.);
    standardConcentrations.resize(kin.nTotalSpecies(), 0.);
    ready = true;
}

InterfaceRateBase::InterfaceRateBase()
    : m_siteDensity(NAN)
    , m_acov(0.)
    , m_ecov(0.)
    , m_mcov(0.)
    , m_chargeTransfer(false)
    , m_exchangeCurrentDensityFormulation(false)
    , m_beta(0.5)
    , m_deltaPotential_RT(NAN)
    , m_deltaGibbs0_RT(NAN)
    , m_prodStandardConcentrations(NAN)
{
}

void InterfaceRateBase::setParameters(const AnyMap& node)
{
    if (node.hasKey("coverage-dependencies")) {
        setCoverageDependencies(
            node["coverage-dependencies"].as<AnyMap>(), node.units());
    }
    if (node.hasKey("beta")) {
        m_beta = node["beta"].asDouble();
    }
    m_exchangeCurrentDensityFormulation = node.getBool(
        "exchange-current-density-formulation", false);
}

void InterfaceRateBase::getParameters(AnyMap& node) const
{
    if (!m_cov.empty()) {
        AnyMap deps;
        getCoverageDependencies(deps);
        node["coverage-dependencies"] = std::move(deps);
    }
    if (m_chargeTransfer) {
        if (m_beta != 0.5) {
            node["beta"] = m_beta;
        }
        if (m_exchangeCurrentDensityFormulation) {
            node["exchange-current-density-formulation"] = true;
        }
    }
}

void InterfaceRateBase::setCoverageDependencies(const AnyMap& dependencies,
                                                const UnitSystem& units)
{
    m_cov.clear();
    m_ac.clear();
    m_ec.clear();
    m_mc.clear();
    m_lindep.clear();
    for (const auto& [species, coeffs] : dependencies) {
        double a, m;
        vector<double> E(5, 0.0);
        if (coeffs.is<AnyMap>()) {
            auto& cov_map = coeffs.as<AnyMap>();
            a = cov_map["a"].asDouble();
            m = cov_map["m"].asDouble();
            if (cov_map["E"].isScalar()) {
                m_lindep.push_back(true);
                E[1] = units.convertActivationEnergy(cov_map["E"], "K");
            } else {
                m_lindep.push_back(false);
                auto& E_temp = cov_map["E"].asVector<AnyValue>(1, 4);
                for (size_t i = 0; i < E_temp.size(); i++) {
                    E[i+1] = units.convertActivationEnergy(E_temp[i], "K");
                }
            }
        } else {
            auto& cov_vec = coeffs.asVector<AnyValue>(3);
            a = cov_vec[0].asDouble();
            m = cov_vec[1].asDouble();
            if (cov_vec[2].isScalar()) {
                m_lindep.push_back(true);
                E[1] = units.convertActivationEnergy(cov_vec[2], "K");
            } else {
                m_lindep.push_back(false);
                auto& E_temp = cov_vec[2].asVector<AnyValue>(1, 4);
                for (size_t i = 0; i < E_temp.size(); i++) {
                    E[i+1] = units.convertActivationEnergy(E_temp[i], "K");
                }
            }
        }
        addCoverageDependence(species, a, m, E);
    }
}

void InterfaceRateBase::getCoverageDependencies(AnyMap& dependencies) const
{
    for (size_t k = 0; k < m_cov.size(); k++) {
        AnyMap dep;
        dep["a"] = m_ac[k];
        dep["m"] = m_mc[k];
        if (m_lindep[k]) {
            dep["E"].setQuantity(m_ec[k][1], "K", true);
        } else {
            vector<AnyValue> E_temp(4);
            for (size_t i = 0; i < m_ec[k].size() - 1; i++) {
                E_temp[i].setQuantity(m_ec[k][i+1], "K", true);
            }
            dep["E"] = E_temp;
        }
        dependencies[m_cov[k]] = std::move(dep);
    }
}

void InterfaceRateBase::addCoverageDependence(const string& sp, double a, double m,
                                              span<const double> e)
{
    if (std::find(m_cov.begin(), m_cov.end(), sp) == m_cov.end()) {
        m_cov.push_back(sp);
        m_ac.push_back(a);
        m_ec.emplace_back(e.begin(), e.end());
        m_mc.push_back(m);
        m_indices.clear();
    } else {
        throw CanteraError("InterfaceRateBase::addCoverageDependence",
            "Coverage for species '{}' is already specified.", sp);
    }
}

void InterfaceRateBase::setSpecies(span<const string> species)
{
    m_indices.clear();
    for (size_t k = 0; k < m_cov.size(); k++) {
        auto it = find(species.begin(), species.end(), m_cov[k]);
        if (it != species.end()) {
            m_indices.emplace(k, it - species.begin());
        } else {
            throw CanteraError("InterfaceRateBase:setSpeciesIndices",
                "Species list does not contain '{}'.", m_cov[k]);
        }
    }
}

void InterfaceRateBase::updateFromStruct(const InterfaceData& shared_data) {
    if (shared_data.ready) {
        m_siteDensity = shared_data.density;
    }

    if (m_indices.size() != m_cov.size()) {
        // object is not set up correctly (setSpecies needs to be run)
        m_acov = NAN;
        m_ecov = NAN;
        m_mcov = NAN;
        return;
    }
    m_acov = 0.0;
    m_ecov = 0.0;
    m_mcov = 0.0;
    for (auto& [iCov, iKin] : m_indices) {
        m_acov += m_ac[iCov] * shared_data.coverages[iKin];
        if (m_lindep[iCov]) {
            m_ecov += m_ec[iCov][1] * shared_data.coverages[iKin];
        } else {
            m_ecov += poly4(shared_data.coverages[iKin], m_ec[iCov].data());
        }
        m_mcov += m_mc[iCov] * shared_data.logCoverages[iKin];
    }

    // Update change in electrical potential energy
    if (m_chargeTransfer) {
        m_deltaPotential_RT = 0.;
        for (const auto& [iPhase, netCharge] : m_netCharges) {
            m_deltaPotential_RT +=
                shared_data.electricPotentials[iPhase] * netCharge;
        }
        m_deltaPotential_RT /= GasConstant * shared_data.temperature;
    }

    // Update quantities used for exchange current density formulation
    if (m_exchangeCurrentDensityFormulation) {
        m_deltaGibbs0_RT = 0.;
        m_prodStandardConcentrations = 1.;
        for (const auto& [k, stoich] : m_stoichCoeffs) {
            m_deltaGibbs0_RT +=
                shared_data.standardChemPotentials[k] * stoich;
            if (stoich > 0.) {
                m_prodStandardConcentrations *=
                    shared_data.standardConcentrations[k];
            }
        }
        m_deltaGibbs0_RT /= GasConstant * shared_data.temperature;
    }
}

void InterfaceRateBase::setContext(const Reaction& rxn, const Kinetics& kin)
{
    setSpecies(kin.thermo().speciesNames());

    m_chargeTransfer = rxn.usesElectrochemistry(kin);
    if (!m_chargeTransfer) {
        return;
    }

    m_stoichCoeffs.clear();
    for (const auto& [name, stoich] : rxn.reactants) {
        m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(name), -stoich);
    }
    for (const auto& [name, stoich] : rxn.products) {
        m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(name), stoich);
    }

    m_netCharges.clear();
    for (const auto& [k, stoich] : m_stoichCoeffs) {
        size_t n = kin.speciesPhaseIndex(k);
        size_t start = kin.kineticsSpeciesIndex(0, n);
        double charge = kin.thermo(n).charge(k - start);
        m_netCharges.emplace_back(n, Faraday * charge * stoich);
    }
}

StickingCoverage::StickingCoverage()
    : m_motzWise(false)
    , m_explicitMotzWise(false)
    , m_stickingSpecies("")
    , m_explicitSpecies(false)
    , m_surfaceOrder(NAN)
    , m_multiplier(NAN)
    , m_factor(NAN)
{
}

void StickingCoverage::setStickingParameters(const AnyMap& node)
{
    m_motzWise = node.getBool("Motz-Wise", false);
    m_explicitMotzWise = node.hasKey("Motz-Wise");
    m_stickingSpecies = node.getString("sticking-species", "");
    m_explicitSpecies = node.hasKey("sticking-species");
}

void StickingCoverage::getStickingParameters(AnyMap& node) const
{
    if (m_explicitMotzWise) {
        node["Motz-Wise"] = m_motzWise;
    }
    if (m_explicitSpecies) {
        node["sticking-species"] = m_stickingSpecies;
    }
}

void StickingCoverage::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // Ensure that site density is initialized
    const ThermoPhase& phase = kin.thermo(0);
    const auto& surf = dynamic_cast<const SurfPhase&>(phase);
    m_siteDensity = surf.siteDensity();
    if (!m_explicitMotzWise) {
        m_motzWise = kin.thermo().input().getBool("Motz-Wise", false);
    }

    string sticking_species = m_stickingSpecies;
    if (sticking_species == "") {
        // Identify the sticking species if not explicitly given
        vector<string> gasSpecies;
        vector<string> anySpecies;
        for (const auto& [name, stoich] : rxn.reactants) {
            size_t iPhase = kin.speciesPhaseIndex(kin.kineticsSpeciesIndex(name));
            if (iPhase != 0) {
                // Non-interface species. There should be exactly one of these
                // (either in gas phase or other phase)
                if (kin.thermo(iPhase).phaseOfMatter() == "gas") {
                    gasSpecies.push_back(name);
                }
                anySpecies.push_back(name);
            }
        }
        if (gasSpecies.size() == 1) {
            // single sticking species in gas phase
            sticking_species = gasSpecies[0];
        } else if (anySpecies.size() == 1) {
            // single sticking species in any phase
            sticking_species = anySpecies[0];
        } else if (anySpecies.size() == 0) {
            throw InputFileError("StickingCoverage::setContext",
                rxn.input, "No non-interface species found "
                "in sticking reaction: '{}'", rxn.equation());
        } else {
            throw InputFileError("StickingCoverage::setContext",
                rxn.input, "Multiple non-interface species ({})\nfound in sticking "
                "reaction: '{}'.\nSticking species must be explicitly specified.",
                fmt::format("'{}'", fmt::join(anySpecies, "', '")), rxn.equation());
        }
    }
    m_stickingSpecies = sticking_species;

    double surface_order = 0.0;
    double multiplier = 1.0;
    // Adjust the A-factor
    for (const auto& [name, stoich] : rxn.reactants) {
        size_t iPhase = kin.speciesPhaseIndex(kin.kineticsSpeciesIndex(name));
        const ThermoPhase& p = kin.thermo(iPhase);
        size_t k = p.speciesIndex(name, true);
        if (name == sticking_species) {
            multiplier *= sqrt(GasConstant / (2 * Pi * p.molecularWeight(k)));
        } else {
            // Non-sticking species. Convert from coverages used in the
            // sticking probability expression to the concentration units
            // used in the mass action rate expression. For surface phases,
            // the dependence on the site density is incorporated when the
            // rate constant is evaluated, since we don't assume that the
            // site density is known at this time.
            double order = getValue(rxn.orders, name, stoich);
            if (&p == &surf) {
                multiplier *= pow(surf.size(k), order);
                surface_order += order;
            } else {
                multiplier *= pow(p.standardConcentration(k), -order);
            }
        }
    }
    m_surfaceOrder = surface_order;
    m_multiplier = multiplier;
}

}
