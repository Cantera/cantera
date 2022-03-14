// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

BulkKinetics::BulkKinetics(ThermoPhase* thermo) :
    m_ROP_ok(false),
    m_temp(0.0)
{
    if (thermo) {
        addPhase(*thermo);
    }
}

void BulkKinetics::resizeReactions()
{
    Kinetics::resizeReactions();

    m_multi_concm.resizeCoeffs(nTotalSpecies(), nReactions());
    for (auto& rates : m_bulk_rates) {
        rates->resize(nTotalSpecies(), nReactions(), nPhases());
        // @todo ensure that ReactionData are updated; calling rates->update
        //      blocks correct behavior in GasKinetics::update_rates_T
        //      and running updateROP() is premature
    }
}

bool BulkKinetics::isReversible(size_t i) {
    return std::find(m_revindex.begin(), m_revindex.end(), i) < m_revindex.end();
}

void BulkKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the ideal gas solution.
    thermo().getChemPotentials(m_grt.data());
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaG);
}

void BulkKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species in the ideal gas.
    thermo().getPartialMolarEnthalpies(m_grt.data());
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void BulkKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in the solid solution.
    thermo().getPartialMolarEntropies(m_grt.data());
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void BulkKinetics::getDeltaSSGibbs(doublereal* deltaG)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity. We define these here as
    // the chemical potentials of the pure species at the temperature and
    // pressure of the solution.
    thermo().getStandardChemPotentials(m_grt.data());
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaG);
}

void BulkKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    // Get the standard state enthalpies of the species.
    thermo().getEnthalpy_RT(m_grt.data());
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= thermo().RT();
    }
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void BulkKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    // Get the standard state entropy of the species. We define these here as
    // the entropies of the pure species at the temperature and pressure of the
    // solution.
    thermo().getEntropy_R(m_grt.data());
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= GasConstant;
    }
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void BulkKinetics::getRevRateConstants(double* krev, bool doIrreversible)
{
    // go get the forward rate constants. -> note, we don't really care about
    // speed or redundancy in these informational routines.
    getFwdRateConstants(krev);

    if (doIrreversible) {
        getEquilibriumConstants(m_ropnet.data());
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] /= m_ropnet[i];
        }
    } else {
        // m_rkcn[] is zero for irreversible reactions
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

bool BulkKinetics::addReaction(shared_ptr<Reaction> r, bool resize)
{
    bool added = Kinetics::addReaction(r, resize);
    if (!added) {
        // undeclared species, etc.
        return false;
    }
    double dn = 0.0;
    for (const auto& sp : r->products) {
        dn += sp.second;
    }
    for (const auto& sp : r->reactants) {
        dn -= sp.second;
    }

    m_dn.push_back(dn);

    if (r->reversible) {
        m_revindex.push_back(nReactions()-1);
    } else {
        m_irrev.push_back(nReactions()-1);
    }

    if (!(r->usesLegacy())) {
        shared_ptr<ReactionRate> rate = r->rate();
        // If necessary, add new MultiRate evaluator
        if (m_bulk_types.find(rate->type()) == m_bulk_types.end()) {
            m_bulk_types[rate->type()] = m_bulk_rates.size();
            m_bulk_rates.push_back(rate->newMultiRate());
            m_bulk_rates.back()->resize(m_kk, nReactions(), nPhases());
        }

        // Set index of rate to number of reaction within kinetics
        rate->setRateIndex(nReactions() - 1);
        rate->setContext(*r, *this);

        // Add reaction rate to evaluator
        size_t index = m_bulk_types[rate->type()];
        m_bulk_rates[index]->add(nReactions() - 1, *rate);

        // Add reaction to third-body evaluator
        if (r->thirdBody() != nullptr) {
            addThirdBody(r);
        }
    }

    m_concm.push_back(NAN);
        m_ready = resize;

    return true;
}

void BulkKinetics::addThirdBody(shared_ptr<Reaction> r)
{
    std::map<size_t, double> efficiencies;
    for (const auto& eff : r->thirdBody()->efficiencies) {
        size_t k = kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            efficiencies[k] = eff.second;
        } else if (!m_skipUndeclaredThirdBodies) {
            throw CanteraError("BulkKinetics::addThirdBody", "Found "
                "third-body efficiency for undefined species '" + eff.first +
                "' while adding reaction '" + r->equation() + "'");
        }
    }
    m_multi_concm.install(nReactions() - 1, efficiencies,
                          r->thirdBody()->default_efficiency,
                          r->thirdBody()->mass_action);
}

void BulkKinetics::addElementaryReaction(ElementaryReaction2& r)
{
    m_rates.install(nReactions()-1, r.rate);
}

void BulkKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    // operations common to all reaction types
    Kinetics::modifyReaction(i, rNew);

    if (!(rNew->usesLegacy())) {
        shared_ptr<ReactionRate> rate = rNew->rate();
        // Ensure that MultiRate evaluator is available
        if (m_bulk_types.find(rate->type()) == m_bulk_types.end()) {
            throw CanteraError("BulkKinetics::modifyReaction",
                 "Evaluator not available for type '{}'.", rate->type());
        }

        // Replace reaction rate to evaluator
        size_t index = m_bulk_types[rate->type()];
        rate->setRateIndex(i);
        rate->setContext(*rNew, *this);

        m_bulk_rates[index]->replace(i, *rate);
    }

    invalidateCache();
}

void BulkKinetics::modifyElementaryReaction(size_t i, ElementaryReaction2& rNew)
{
    m_rates.replace(i, rNew.rate);
}

void BulkKinetics::resizeSpecies()
{
    Kinetics::resizeSpecies();
    m_act_conc.resize(m_kk);
    m_phys_conc.resize(m_kk);
    m_grt.resize(m_kk);
    for (auto& rates : m_bulk_rates) {
        rates->resize(m_kk, nReactions(), nPhases());
    }
}

void BulkKinetics::setMultiplier(size_t i, double f) {
    Kinetics::setMultiplier(i, f);
    m_ROP_ok = false;
}

void BulkKinetics::invalidateCache()
{
    Kinetics::invalidateCache();
    m_ROP_ok = false;
    m_temp += 0.13579;
}

}
