// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

BulkKinetics::BulkKinetics(thermo_t* thermo) :
    m_ROP_ok(false),
    m_temp(0.0)
{
    if (thermo) {
        addPhase(*thermo);
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

void BulkKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
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

bool BulkKinetics::addReaction(shared_ptr<Reaction> r)
{
    bool added = Kinetics::addReaction(r);
    if (!added) {
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
    return true;
}

void BulkKinetics::addElementaryReaction(ElementaryReaction& r)
{
    m_rates.install(nReactions()-1, r.rate);
}

void BulkKinetics::modifyElementaryReaction(size_t i, ElementaryReaction& rNew)
{
    m_rates.replace(i, rNew.rate);
}

void BulkKinetics::resizeSpecies()
{
    Kinetics::resizeSpecies();
    m_conc.resize(m_kk);
    m_grt.resize(m_kk);
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
