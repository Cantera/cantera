/**
 *  @file Kinetics.cpp Declarations for the base class for kinetics managers
 *      (see \ref  kineticsmgr and class \link Cantera::Kinetics  Kinetics \endlink).
 *
 *  Kinetics managers calculate rates of progress of species due to
 *  homogeneous or heterogeneous kinetics.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/stringUtils.h"
#include <unordered_set>

using namespace std;

namespace Cantera
{
Kinetics::Kinetics() :
    m_kk(0),
    m_thermo(0),
    m_surfphase(npos),
    m_rxnphase(npos),
    m_mindim(4),
    m_skipUndeclaredSpecies(false),
    m_skipUndeclaredThirdBodies(false)
{
}

Kinetics::~Kinetics() {}

void Kinetics::checkReactionIndex(size_t i) const
{
    if (i >= nReactions()) {
        throw IndexError("Kinetics::checkReactionIndex", "reactions", i,
                         nReactions()-1);
    }
}

void Kinetics::checkReactionArraySize(size_t ii) const
{
    if (nReactions() > ii) {
        throw ArraySizeError("Kinetics::checkReactionArraySize", ii,
                             nReactions());
    }
}

void Kinetics::checkPhaseIndex(size_t m) const
{
    if (m >= nPhases()) {
        throw IndexError("Kinetics::checkPhaseIndex", "phase", m, nPhases()-1);
    }
}

void Kinetics::checkPhaseArraySize(size_t mm) const
{
    if (nPhases() > mm) {
        throw ArraySizeError("Kinetics::checkPhaseArraySize", mm, nPhases());
    }
}

void Kinetics::checkSpeciesIndex(size_t k) const
{
    if (k >= m_kk) {
        throw IndexError("Kinetics::checkSpeciesIndex", "species", k, m_kk-1);
    }
}

void Kinetics::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("Kinetics::checkSpeciesArraySize", kk, m_kk);
    }
}

std::pair<size_t, size_t> Kinetics::checkDuplicates(bool throw_err) const
{
    //! Map of (key indicating participating species) to reaction numbers
    std::map<size_t, std::vector<size_t> > participants;
    std::vector<std::map<int, double> > net_stoich;
    std::unordered_set<size_t> unmatched_duplicates;
    for (size_t i = 0; i < m_reactions.size(); i++) {
        if (m_reactions[i]->duplicate) {
            unmatched_duplicates.insert(i);
        }
    }

    for (size_t i = 0; i < m_reactions.size(); i++) {
        // Get data about this reaction
        unsigned long int key = 0;
        Reaction& R = *m_reactions[i];
        net_stoich.emplace_back();
        std::map<int, double>& net = net_stoich.back();
        for (const auto& sp : R.reactants) {
            int k = static_cast<int>(kineticsSpeciesIndex(sp.first));
            key += k*(k+1);
            net[-1 -k] -= sp.second;
        }
        for (const auto& sp : R.products) {
            int k = static_cast<int>(kineticsSpeciesIndex(sp.first));
            key += k*(k+1);
            net[1+k] += sp.second;
        }

        // Compare this reaction to others with similar participants
        vector<size_t>& related = participants[key];
        for (size_t m : related) {
            Reaction& other = *m_reactions[m];
            if (R.duplicate && other.duplicate) {
                // marked duplicates
                unmatched_duplicates.erase(i);
                unmatched_duplicates.erase(m);
                continue;
            } else if (R.reaction_type != other.reaction_type) {
                continue; // different reaction types
            }
            doublereal c = checkDuplicateStoich(net_stoich[i], net_stoich[m]);
            if (c == 0) {
                continue; // stoichiometries differ (not by a multiple)
            } else if (c < 0.0 && !R.reversible && !other.reversible) {
                continue; // irreversible reactions in opposite directions
            } else if (R.reaction_type == FALLOFF_RXN ||
                       R.reaction_type == CHEMACT_RXN) {
                ThirdBody& tb1 = dynamic_cast<FalloffReaction&>(R).third_body;
                ThirdBody& tb2 = dynamic_cast<FalloffReaction&>(other).third_body;
                bool thirdBodyOk = true;
                for (size_t k = 0; k < nTotalSpecies(); k++) {
                    string s = kineticsSpeciesName(k);
                    if (tb1.efficiency(s) * tb2.efficiency(s) != 0.0) {
                        // non-zero third body efficiencies for species `s` in
                        // both reactions
                        thirdBodyOk = false;
                        break;
                    }
                }
                if (thirdBodyOk) {
                    continue; // No overlap in third body efficiencies
                }
            } else if (R.reaction_type == THREE_BODY_RXN) {
                ThirdBody& tb1 = dynamic_cast<ThreeBodyReaction&>(R).third_body;
                ThirdBody& tb2 = dynamic_cast<ThreeBodyReaction&>(other).third_body;
                bool thirdBodyOk = true;
                for (size_t k = 0; k < nTotalSpecies(); k++) {
                    string s = kineticsSpeciesName(k);
                    if (tb1.efficiency(s) * tb2.efficiency(s) != 0.0) {
                        // non-zero third body efficiencies for species `s` in
                        // both reactions
                        thirdBodyOk = false;
                        break;
                    }
                }
                if (thirdBodyOk) {
                    continue; // No overlap in third body efficiencies
                }
            }
            if (throw_err) {
                throw InputFileError("Kinetics::checkDuplicates",
                        R.input, other.input,
                        "Undeclared duplicate reactions detected:\n"
                        "Reaction {}: {}\nReaction {}: {}\n",
                        i+1, other.equation(), m+1, R.equation());
            } else {
                return {i,m};
            }
        }
        participants[key].push_back(i);
    }
    if (unmatched_duplicates.size()) {
        size_t i = *unmatched_duplicates.begin();
        if (throw_err) {
            throw InputFileError("Kinetics::checkDuplicates",
                m_reactions[i]->input,
                "No duplicate found for declared duplicate reaction number {}"
                " ({})", i, m_reactions[i]->equation());
        } else {
            return {i, i};
        }
    }
    return {npos, npos};
}

double Kinetics::checkDuplicateStoich(std::map<int, double>& r1,
                                      std::map<int, double>& r2) const
{
    std::unordered_set<int> keys; // species keys (k+1 or -k-1)
    for (auto& r : r1) {
        keys.insert(r.first);
    }
    for (auto& r : r2) {
        keys.insert(r.first);
    }
    int k1 = r1.begin()->first;
    // check for duplicate written in the same direction
    doublereal ratio = 0.0;
    if (r1[k1] && r2[k1]) {
        ratio = r2[k1]/r1[k1];
        bool different = false;
        for (int k : keys) {
            if ((r1[k] && !r2[k]) ||
                (!r1[k] && r2[k]) ||
                (r1[k] && fabs(r2[k]/r1[k] - ratio) > 1.e-8)) {
                different = true;
                break;
            }
        }
        if (!different) {
            return ratio;
        }
    }

    // check for duplicate written in the reverse direction
    if (r1[k1] == 0.0 || r2[-k1] == 0.0) {
        return 0.0;
    }
    ratio = r2[-k1]/r1[k1];
    for (int k : keys) {
        if ((r1[k] && !r2[-k]) ||
            (!r1[k] && r2[-k]) ||
            (r1[k] && fabs(r2[-k]/r1[k] - ratio) > 1.e-8)) {
            return 0.0;
        }
    }
    return ratio;
}

void Kinetics::checkReactionBalance(const Reaction& R)
{
    Composition balr, balp;
    // iterate over the products
    for (const auto& sp : R.products) {
        const ThermoPhase& ph = speciesPhase(sp.first);
        size_t k = ph.speciesIndex(sp.first);
        double stoich = sp.second;
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] = 0.0; // so that balr contains all species
            balp[ph.elementName(m)] += stoich*ph.nAtoms(k,m);
        }
    }
    for (const auto& sp : R.reactants) {
        const ThermoPhase& ph = speciesPhase(sp.first);
        size_t k = ph.speciesIndex(sp.first);
        double stoich = sp.second;
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] += stoich*ph.nAtoms(k,m);
        }
    }

    string msg;
    bool ok = true;
    for (const auto& el : balr) {
        const string& elem = el.first;
        double elemsum = balr[elem] + balp[elem];
        double elemdiff = fabs(balp[elem] - balr[elem]);
        if (elemsum > 0.0 && elemdiff/elemsum > 1e-4) {
            ok = false;
            msg += fmt::format("  {}           {}           {}\n",
                               elem, balr[elem], balp[elem]);
        }
    }
    if (!ok) {
        throw InputFileError("Kinetics::checkReactionBalance", R.input,
            "The following reaction is unbalanced: {}\n"
            "  Element    Reactants    Products\n{}",
            R.equation(), msg);
    }
}

void Kinetics::selectPhase(const doublereal* data, const thermo_t* phase,
                           doublereal* phase_data)
{
    for (size_t n = 0; n < nPhases(); n++) {
        if (phase == m_thermo[n]) {
            size_t nsp = phase->nSpecies();
            copy(data + m_start[n],
                 data + m_start[n] + nsp, phase_data);
            return;
        }
    }
    throw CanteraError("Kinetics::selectPhase", "Phase not found.");
}

string Kinetics::kineticsSpeciesName(size_t k) const
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return thermo(n).speciesName(k - m_start[n]);
        }
    }
    return "<unknown>";
}

size_t Kinetics::kineticsSpeciesIndex(const std::string& nm) const
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        // Check the ThermoPhase object for a match
        size_t k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return k + m_start[n];
        }
    }
    return npos;
}

size_t Kinetics::kineticsSpeciesIndex(const std::string& nm,
                                      const std::string& ph) const
{
    if (ph == "<any>") {
        return kineticsSpeciesIndex(nm);
    }

    for (size_t n = 0; n < m_thermo.size(); n++) {
        string id = thermo(n).name();
        if (ph == id) {
            size_t k = thermo(n).speciesIndex(nm);
            if (k == npos) {
                return npos;
            }
            return k + m_start[n];
        }
    }
    return npos;
}

thermo_t& Kinetics::speciesPhase(const std::string& nm)
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        size_t k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return thermo(n);
        }
    }
    throw CanteraError("Kinetics::speciesPhase", "unknown species '{}'", nm);
}

const thermo_t& Kinetics::speciesPhase(const std::string& nm) const
{
    for (const auto thermo : m_thermo) {
        if (thermo->speciesIndex(nm) != npos) {
            return *thermo;
        }
    }
    throw CanteraError("Kinetics::speciesPhase", "unknown species '{}'", nm);
}

size_t Kinetics::speciesPhaseIndex(size_t k) const
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return n;
        }
    }
    throw CanteraError("Kinetics::speciesPhaseIndex",
                       "illegal species index: {}", k);
}

double Kinetics::reactantStoichCoeff(size_t kSpec, size_t irxn) const
{
    return getValue(m_reactions[irxn]->reactants, kineticsSpeciesName(kSpec),
                    0.0);
}

double Kinetics::productStoichCoeff(size_t kSpec, size_t irxn) const
{
    return getValue(m_reactions[irxn]->products, kineticsSpeciesName(kSpec),
                    0.0);
}

void Kinetics::getFwdRatesOfProgress(doublereal* fwdROP)
{
    updateROP();
    std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
}

void Kinetics::getRevRatesOfProgress(doublereal* revROP)
{
    updateROP();
    std::copy(m_ropr.begin(), m_ropr.end(), revROP);
}

void Kinetics::getNetRatesOfProgress(doublereal* netROP)
{
    updateROP();
    std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
}

void Kinetics::getReactionDelta(const double* prop, double* deltaProp)
{
    fill(deltaProp, deltaProp + nReactions(), 0.0);
    // products add
    m_revProductStoich.incrementReactions(prop, deltaProp);
    m_irrevProductStoich.incrementReactions(prop, deltaProp);
    // reactants subtract
    m_reactantStoich.decrementReactions(prop, deltaProp);
}

void Kinetics::getRevReactionDelta(const double* prop, double* deltaProp)
{
    fill(deltaProp, deltaProp + nReactions(), 0.0);
    // products add
    m_revProductStoich.incrementReactions(prop, deltaProp);
    // reactants subtract
    m_reactantStoich.decrementReactions(prop, deltaProp);
}

void Kinetics::getCreationRates(double* cdot)
{
    updateROP();

    // zero out the output array
    fill(cdot, cdot + m_kk, 0.0);

    // the forward direction creates product species
    m_revProductStoich.incrementSpecies(m_ropf.data(), cdot);
    m_irrevProductStoich.incrementSpecies(m_ropf.data(), cdot);

    // the reverse direction creates reactant species
    m_reactantStoich.incrementSpecies(m_ropr.data(), cdot);
}

void Kinetics::getDestructionRates(doublereal* ddot)
{
    updateROP();

    fill(ddot, ddot + m_kk, 0.0);
    // the reverse direction destroys products in reversible reactions
    m_revProductStoich.incrementSpecies(m_ropr.data(), ddot);
    // the forward direction destroys reactants
    m_reactantStoich.incrementSpecies(m_ropf.data(), ddot);
}

void Kinetics::getNetProductionRates(doublereal* net)
{
    updateROP();

    fill(net, net + m_kk, 0.0);
    // products are created for positive net rate of progress
    m_revProductStoich.incrementSpecies(m_ropnet.data(), net);
    m_irrevProductStoich.incrementSpecies(m_ropnet.data(), net);
    // reactants are destroyed for positive net rate of progress
    m_reactantStoich.decrementSpecies(m_ropnet.data(), net);
}

void Kinetics::addPhase(thermo_t& thermo)
{
    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo.nDim() <= m_mindim) {
        m_mindim = thermo.nDim();
        m_rxnphase = nPhases();
    }

    // there should only be one surface phase
    if (thermo.type() == kineticsType()) {
        m_surfphase = nPhases();
        m_rxnphase = nPhases();
    }
    m_thermo.push_back(&thermo);
    m_phaseindex[m_thermo.back()->name()] = nPhases();
    resizeSpecies();
}

void Kinetics::resizeSpecies()
{
    m_kk = 0;
    m_start.resize(nPhases());

    for (size_t i = 0; i < m_thermo.size(); i++) {
        m_start[i] = m_kk; // global index of first species of phase i
        m_kk += m_thermo[i]->nSpecies();
    }
    invalidateCache();
}

bool Kinetics::addReaction(shared_ptr<Reaction> r)
{
    r->validate();
    if (m_kk == 0) {
        init();
    }
    resizeSpecies();

    // If reaction orders are specified, then this reaction does not follow
    // mass-action kinetics, and is not an elementary reaction. So check that it
    // is not reversible, since computing the reverse rate from thermochemistry
    // only works for elementary reactions.
    if (r->reversible && !r->orders.empty()) {
        throw InputFileError("Kinetics::addReaction", r->input,
            "Reaction orders may only be given for irreversible reactions");
    }

    // Check for undeclared species
    for (const auto& sp : r->reactants) {
        if (kineticsSpeciesIndex(sp.first) == npos) {
            if (m_skipUndeclaredSpecies) {
                return false;
            } else {
                throw InputFileError("Kinetics::addReaction", r->input,
                    "Reaction '{}' contains the undeclared species '{}'",
                    r->equation(), sp.first);
            }
        }
    }
    for (const auto& sp : r->products) {
        if (kineticsSpeciesIndex(sp.first) == npos) {
            if (m_skipUndeclaredSpecies) {
                return false;
            } else {
                throw InputFileError("Kinetics::addReaction", r->input,
                    "Reaction '{}' contains the undeclared species '{}'",
                    r->equation(), sp.first);
            }
        }
    }
    for (const auto& sp : r->orders) {
        if (kineticsSpeciesIndex(sp.first) == npos) {
            if (m_skipUndeclaredSpecies) {
                return false;
            } else {
                throw InputFileError("Kinetics::addReaction", r->input,
                    "Reaction '{}' has a reaction order specified for the "
                    "undeclared species '{}'",
                    r->equation(), sp.first);
            }
        }
    }

    checkReactionBalance(*r);
    size_t irxn = nReactions(); // index of the new reaction

    // indices of reactant and product species within this Kinetics object
    std::vector<size_t> rk, pk;

    // Reactant and product stoichiometric coefficients, such that rstoich[i] is
    // the coefficient for species rk[i]
    vector_fp rstoich, pstoich;

    for (const auto& sp : r->reactants) {
        rk.push_back(kineticsSpeciesIndex(sp.first));
        rstoich.push_back(sp.second);
    }

    for (const auto& sp : r->products) {
        pk.push_back(kineticsSpeciesIndex(sp.first));
        pstoich.push_back(sp.second);
    }

    // The default order for each reactant is its stoichiometric coefficient,
    // which can be overridden by entries in the Reaction.orders map. rorder[i]
    // is the order for species rk[i].
    vector_fp rorder = rstoich;
    for (const auto& sp : r->orders) {
        size_t k = kineticsSpeciesIndex(sp.first);
        // Find the index of species k within rk
        auto rloc = std::find(rk.begin(), rk.end(), k);
        if (rloc != rk.end()) {
            rorder[rloc - rk.begin()] = sp.second;
        } else {
            // If the reaction order involves a non-reactant species, add an
            // extra term to the reactants with zero stoichiometry so that the
            // stoichiometry manager can be used to compute the global forward
            // reaction rate.
            rk.push_back(k);
            rstoich.push_back(0.0);
            rorder.push_back(sp.second);
        }
    }

    m_reactantStoich.add(irxn, rk, rorder, rstoich);
    // product orders = product stoichiometric coefficients
    if (r->reversible) {
        m_revProductStoich.add(irxn, pk, pstoich, pstoich);
    } else {
        m_irrevProductStoich.add(irxn, pk, pstoich, pstoich);
    }

    m_reactions.push_back(r);
    m_rfn.push_back(0.0);
    m_rkcn.push_back(0.0);
    m_ropf.push_back(0.0);
    m_ropr.push_back(0.0);
    m_ropnet.push_back(0.0);
    m_perturb.push_back(1.0);
    return true;
}

void Kinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    checkReactionIndex(i);
    shared_ptr<Reaction>& rOld = m_reactions[i];
    if (rNew->reaction_type != rOld->reaction_type) {
        throw CanteraError("Kinetics::modifyReaction",
            "Reaction types are different: {} != {}.",
            rOld->reaction_type, rNew->reaction_type);
    }

    if (rNew->reactants != rOld->reactants) {
        throw CanteraError("Kinetics::modifyReaction",
            "Reactants are different: '{}' != '{}'.",
            rOld->reactantString(), rNew->reactantString());
    }

    if (rNew->products != rOld->products) {
        throw CanteraError("Kinetics::modifyReaction",
            "Products are different: '{}' != '{}'.",
            rOld->productString(), rNew->productString());
    }
    m_reactions[i] = rNew;
    invalidateCache();
}

shared_ptr<Reaction> Kinetics::reaction(size_t i)
{
    checkReactionIndex(i);
    return m_reactions[i];
}

shared_ptr<const Reaction> Kinetics::reaction(size_t i) const
{
    checkReactionIndex(i);
    return m_reactions[i];
}

}
