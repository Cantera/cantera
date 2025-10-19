/**
 *  @file Kinetics.cpp Declarations for the base class for kinetics managers
 *      (see @ref  kineticsmgr and class @link Cantera::Kinetics  Kinetics @endlink).
 *
 *  Kinetics managers calculate rates of progress of species due to
 *  homogeneous or heterogeneous kinetics.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"
#include <unordered_set>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace Cantera
{

shared_ptr<Kinetics> Kinetics::clone(
    const vector<shared_ptr<ThermoPhase>>& phases) const
{
    vector<AnyMap> reactionDefs;
    for (size_t i = 0; i < nReactions(); i++) {
        reactionDefs.push_back(reaction(i)->parameters());
    }
    AnyMap phaseNode = parameters();
    phaseNode["__fix-duplicate-reactions__"] = true;
    AnyMap rootNode;
    rootNode["reactions"] = std::move(reactionDefs);
    rootNode.applyUnits();
    return newKinetics(phases, phaseNode, rootNode, phases[0]->root());
}

size_t Kinetics::checkReactionIndex(size_t i) const
{
    if (i < nReactions()) {
        return i;
    }
    throw IndexError("Kinetics::checkReactionIndex", "reactions", i, nReactions());
}

void Kinetics::resizeReactions()
{
    size_t nRxn = nReactions();

    // Stoichiometry managers
    m_reactantStoich.resizeCoeffs(m_kk, nRxn);
    m_productStoich.resizeCoeffs(m_kk, nRxn);
    m_revProductStoich.resizeCoeffs(m_kk, nRxn);

    m_rbuf.resize(nRxn);

    // products are created for positive net rate of progress
    m_stoichMatrix = m_productStoich.stoichCoeffs();
    // reactants are destroyed for positive net rate of progress
    m_stoichMatrix -= m_reactantStoich.stoichCoeffs();
}

void Kinetics::checkReactionArraySize(size_t ii) const
{
    warn_deprecated("Kinetics::checkReactionArraySize",
        "To be removed after Cantera 3.2. Only used by legacy CLib.");
    if (nReactions() > ii) {
        throw ArraySizeError("Kinetics::checkReactionArraySize", ii,
                             nReactions());
    }
}

size_t Kinetics::checkPhaseIndex(size_t m) const
{
    if (m < nPhases()) {
        return m;
    }
    throw IndexError("Kinetics::checkPhaseIndex", "phase", m, nPhases());
}

void Kinetics::checkPhaseArraySize(size_t mm) const
{
    warn_deprecated("Kinetics::checkPhaseArraySize",
        "To be removed after Cantera 3.2. Unused.");
    if (nPhases() > mm) {
        throw ArraySizeError("Kinetics::checkPhaseArraySize", mm, nPhases());
    }
}

size_t Kinetics::phaseIndex(const string& ph) const
{
    size_t ix = phaseIndex(ph, false);
    if (ix == npos) {
        warn_deprecated("Kinetics::phaseIndex", "'raise' argument not specified; "
            "Default behavior will change from returning -1 to throwing an "
            "exception after Cantera 3.2.");
    }
    return ix;
}

size_t Kinetics::phaseIndex(const string& ph, bool raise) const
{
    if (m_phaseindex.find(ph) == m_phaseindex.end()) {
        if (raise) {
            throw CanteraError("Kinetics::phaseIndex", "Phase '{}' not found", ph);
        }
        return npos;
    } else {
        return m_phaseindex.at(ph) - 1;
    }
}

shared_ptr<ThermoPhase> Kinetics::reactionPhase() const
{
    return m_thermo[0];
}

size_t Kinetics::checkSpeciesIndex(size_t k) const
{
    if (k < m_kk) {
        return k;
    }
    throw IndexError("Kinetics::checkSpeciesIndex", "species", k, m_kk);
}

void Kinetics::checkSpeciesArraySize(size_t kk) const
{
    warn_deprecated("Kinetics::checkSpeciesArraySize",
        "To be removed after Cantera 3.2. Only used by legacy CLib.");
    if (m_kk > kk) {
        throw ArraySizeError("Kinetics::checkSpeciesArraySize", kk, m_kk);
    }
}

void Kinetics::setExplicitThirdBodyDuplicateHandling(const string& flag)
{
    if (flag == "warn" || flag == "error" || flag == "mark-duplicate"
        || flag == "modify-efficiency")
    {
        m_explicit_third_body_duplicates = flag;
    } else {
        throw CanteraError("Kinetics::setExplicitThirdBodyDuplicateHandling",
            "Invalid flag '{}'", flag);
    }
}

pair<size_t, size_t> Kinetics::checkDuplicates(bool throw_err, bool fix)
{
    //! Map of (key indicating participating species) to reaction numbers
    map<size_t, vector<size_t>> participants;
    vector<map<int, double>> net_stoich;
    std::unordered_set<size_t> unmatched_duplicates;
    for (size_t i = 0; i < m_reactions.size(); i++) {
        if (m_reactions[i]->duplicate) {
            unmatched_duplicates.insert(i);
        }
    }

    vector<InputFileError> errs;
    for (size_t i = 0; i < m_reactions.size(); i++) {
        // Get data about this reaction
        unsigned long int key = 0;
        Reaction& R = *m_reactions[i];
        net_stoich.emplace_back();
        map<int, double>& net = net_stoich.back();
        for (const auto& [name, stoich] : R.reactants) {
            int k = static_cast<int>(kineticsSpeciesIndex(name));
            key += k*(k+1);
            net[-1 -k] -= stoich;
        }
        for (const auto& [name, stoich] : R.products) {
            int k = static_cast<int>(kineticsSpeciesIndex(name));
            key += k*(k+1);
            net[1+k] += stoich;
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
            } else if (R.type() != other.type()) {
                continue; // different reaction types
            } else if (R.rate() && other.rate()
                       && R.rate()->type() != other.rate()->type())
            {
                continue; // different rate parameterizations
            }
            double c = checkDuplicateStoich(net_stoich[i], net_stoich[m]);
            if (c == 0) {
                continue; // stoichiometries differ (not by a multiple)
            } else if (c < 0.0 && !R.reversible && !other.reversible) {
                continue; // irreversible reactions in opposite directions
            } else if (R.usesThirdBody() && other.usesThirdBody()) {
                ThirdBody& tb1 = *(R.thirdBody());
                ThirdBody& tb2 = *(other.thirdBody());
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
                } else if ((tb1.name() == "M") != (tb2.name() == "M")) {
                    // Exactly one of the reactions uses M as the third body
                    if (m_explicit_third_body_duplicates == "mark-duplicate") {
                        R.duplicate = true;
                        other.duplicate = true;
                        continue;
                    } else if (m_explicit_third_body_duplicates == "modify-efficiency") {
                        if (tb1.name() == "M") {
                            tb1.efficiencies[tb2.name()] = 0.0;
                        } else {
                            tb2.efficiencies[tb1.name()] = 0.0;
                        }
                        continue;
                    } else if (m_explicit_third_body_duplicates == "warn") {
                        InputFileError msg("Kinetics::checkDuplicates",
                            R.input, other.input,
                            "Undeclared duplicate third body reactions with a common "
                            "third body detected.\nAdd the field "
                            "'explicit-third-body-duplicates: mark-duplicate' or\n"
                            "'explicit-third-body-duplicates: modify-efficiency' to "
                            "the YAML phase entry\nto choose how these reactions "
                            "should be handled and suppress this warning.\n"
                            "Reaction {}: {}\nReaction {}: {}\n",
                            m+1, R.equation(), i+1, other.equation());
                        if (!fix) {
                            warn_user("Kinetics::checkDuplicates", msg.what());
                        }
                        continue;
                    } // else m_explicit_third_body_duplicates == "error"
                }
            }
            if (throw_err) {
                errs.emplace_back("Kinetics::checkDuplicates",
                        R.input, other.input,
                        "Undeclared duplicate reactions detected:\n"
                        "Reaction {}: {}\nReaction {}: {}\n",
                        m+1, R.equation(), i+1, other.equation());
            } else if (fix) {
                R.duplicate = true;
                other.duplicate = true;
                unmatched_duplicates.erase(i);
                unmatched_duplicates.erase(m);
            } else {
                return {i,m};
            }
        }
        participants[key].push_back(i);
    }
    if (unmatched_duplicates.size()) {
        for (auto i : unmatched_duplicates) {
            if (throw_err) {
                errs.emplace_back("Kinetics::checkDuplicates",
                    m_reactions[i]->input,
                    "No duplicate found for declared duplicate reaction number {}"
                    " ({})", i, m_reactions[i]->equation());
            } else if (fix) {
                m_reactions[i]->duplicate = false;
            } else {
                return {i, i};
            }
        }
    }
    if (errs.empty()) {
        return {npos, npos};
    } else if (errs.size() == 1) {
        throw errs[0];
    } else {
        fmt::memory_buffer msg;
        for (const auto& err : errs) {
            fmt_append(msg, "\n{}\n", err.getMessage());
        }
        throw CanteraError("Kinetics::checkDuplicates", to_string(msg));
    }
}

double Kinetics::checkDuplicateStoich(map<int, double>& r1, map<int, double>& r2) const
{
    std::unordered_set<int> keys; // species keys (k+1 or -k-1)
    for (auto& [speciesKey, stoich] : r1) {
        keys.insert(speciesKey);
    }
    for (auto& [speciesKey, stoich] : r2) {
        keys.insert(speciesKey);
    }
    int k1 = r1.begin()->first;
    // check for duplicate written in the same direction
    double ratio = 0.0;
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

string Kinetics::kineticsSpeciesName(size_t k) const
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return thermo(n).speciesName(k - m_start[n]);
        }
    }
    return "<unknown>";
}

size_t Kinetics::kineticsSpeciesIndex(const string& nm) const
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        // Check the ThermoPhase object for a match
        size_t k = thermo(n).speciesIndex(nm, false);
        if (k != npos) {
            return k + m_start[n];
        }
    }
    return npos;
}

ThermoPhase& Kinetics::speciesPhase(const string& nm)
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        size_t k = thermo(n).speciesIndex(nm, false);
        if (k != npos) {
            return thermo(n);
        }
    }
    throw CanteraError("Kinetics::speciesPhase", "unknown species '{}'", nm);
}

const ThermoPhase& Kinetics::speciesPhase(const string& nm) const
{
    for (const auto& thermo : m_thermo) {
        if (thermo->speciesIndex(nm, false) != npos) {
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
    return m_reactantStoich.stoichCoeffs().coeff(kSpec, irxn);
}

double Kinetics::productStoichCoeff(size_t kSpec, size_t irxn) const
{
    return m_productStoich.stoichCoeffs().coeff(kSpec, irxn);
}

void Kinetics::getFwdRatesOfProgress(double* fwdROP)
{
    updateROP();
    std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
}

void Kinetics::getRevRatesOfProgress(double* revROP)
{
    updateROP();
    std::copy(m_ropr.begin(), m_ropr.end(), revROP);
}

void Kinetics::getNetRatesOfProgress(double* netROP)
{
    updateROP();
    std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
}

void Kinetics::getReactionDelta(const double* prop, double* deltaProp) const
{
    fill(deltaProp, deltaProp + nReactions(), 0.0);
    // products add
    m_productStoich.incrementReactions(prop, deltaProp);
    // reactants subtract
    m_reactantStoich.decrementReactions(prop, deltaProp);
}

void Kinetics::getRevReactionDelta(const double* prop, double* deltaProp) const
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
    m_productStoich.incrementSpecies(m_ropf.data(), cdot);

    // the reverse direction creates reactant species
    m_reactantStoich.incrementSpecies(m_ropr.data(), cdot);
}

void Kinetics::getDestructionRates(double* ddot)
{
    updateROP();

    fill(ddot, ddot + m_kk, 0.0);
    // the reverse direction destroys products in reversible reactions
    m_revProductStoich.incrementSpecies(m_ropr.data(), ddot);
    // the forward direction destroys reactants
    m_reactantStoich.incrementSpecies(m_ropf.data(), ddot);
}

void Kinetics::getNetProductionRates(double* net)
{
    updateROP();

    fill(net, net + m_kk, 0.0);
    // products are created for positive net rate of progress
    m_productStoich.incrementSpecies(m_ropnet.data(), net);
    // reactants are destroyed for positive net rate of progress
    m_reactantStoich.decrementSpecies(m_ropnet.data(), net);
}

void Kinetics::getCreationRates_ddT(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the forward direction creates product species
    getFwdRatesOfProgress_ddT(buf.data());
    out = m_productStoich.stoichCoeffs() * buf;
    // the reverse direction creates reactant species
    getRevRatesOfProgress_ddT(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

void Kinetics::getCreationRates_ddP(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the forward direction creates product species
    getFwdRatesOfProgress_ddP(buf.data());
    out = m_productStoich.stoichCoeffs() * buf;
    // the reverse direction creates reactant species
    getRevRatesOfProgress_ddP(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

void Kinetics::getCreationRates_ddC(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the forward direction creates product species
    getFwdRatesOfProgress_ddC(buf.data());
    out = m_productStoich.stoichCoeffs() * buf;
    // the reverse direction creates reactant species
    getRevRatesOfProgress_ddC(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

Eigen::SparseMatrix<double> Kinetics::creationRates_ddX()
{
    Eigen::SparseMatrix<double> jac;
    // the forward direction creates product species
    jac = m_productStoich.stoichCoeffs() * fwdRatesOfProgress_ddX();
    // the reverse direction creates reactant species
    jac += m_reactantStoich.stoichCoeffs() * revRatesOfProgress_ddX();
    return jac;
}

Eigen::SparseMatrix<double> Kinetics::creationRates_ddCi()
{
    Eigen::SparseMatrix<double> jac;
    // the forward direction creates product species
    jac = m_productStoich.stoichCoeffs() * fwdRatesOfProgress_ddCi();
    // the reverse direction creates reactant species
    jac += m_reactantStoich.stoichCoeffs() * revRatesOfProgress_ddCi();
    return jac;
}

void Kinetics::getDestructionRates_ddT(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the reverse direction destroys products in reversible reactions
    getRevRatesOfProgress_ddT(buf.data());
    out = m_revProductStoich.stoichCoeffs() * buf;
    // the forward direction destroys reactants
    getFwdRatesOfProgress_ddT(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

void Kinetics::getDestructionRates_ddP(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the reverse direction destroys products in reversible reactions
    getRevRatesOfProgress_ddP(buf.data());
    out = m_revProductStoich.stoichCoeffs() * buf;
    // the forward direction destroys reactants
    getFwdRatesOfProgress_ddP(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

void Kinetics::getDestructionRates_ddC(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    // the reverse direction destroys products in reversible reactions
    getRevRatesOfProgress_ddC(buf.data());
    out = m_revProductStoich.stoichCoeffs() * buf;
    // the forward direction destroys reactants
    getFwdRatesOfProgress_ddC(buf.data());
    out += m_reactantStoich.stoichCoeffs() * buf;
}

Eigen::SparseMatrix<double> Kinetics::destructionRates_ddX()
{
    Eigen::SparseMatrix<double> jac;
    // the reverse direction destroys products in reversible reactions
    jac = m_revProductStoich.stoichCoeffs() * revRatesOfProgress_ddX();
    // the forward direction destroys reactants
    jac += m_reactantStoich.stoichCoeffs() * fwdRatesOfProgress_ddX();
    return jac;
}

Eigen::SparseMatrix<double> Kinetics::destructionRates_ddCi()
{
    Eigen::SparseMatrix<double> jac;
    // the reverse direction destroys products in reversible reactions
    jac = m_revProductStoich.stoichCoeffs() * revRatesOfProgress_ddCi();
    // the forward direction destroys reactants
    jac += m_reactantStoich.stoichCoeffs() * fwdRatesOfProgress_ddCi();
    return jac;
}

void Kinetics::getNetProductionRates_ddT(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    getNetRatesOfProgress_ddT(buf.data());
    out = m_stoichMatrix * buf;
}

void Kinetics::getNetProductionRates_ddP(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    getNetRatesOfProgress_ddP(buf.data());
    out = m_stoichMatrix * buf;
}

void Kinetics::getNetProductionRates_ddC(double* dwdot)
{
    Eigen::Map<Eigen::VectorXd> out(dwdot, m_kk);
    Eigen::Map<Eigen::VectorXd> buf(m_rbuf.data(), nReactions());
    getNetRatesOfProgress_ddC(buf.data());
    out = m_stoichMatrix * buf;
}

Eigen::SparseMatrix<double> Kinetics::netProductionRates_ddX()
{
    return m_stoichMatrix * netRatesOfProgress_ddX();
}

Eigen::SparseMatrix<double> Kinetics::netProductionRates_ddCi()
{
    return m_stoichMatrix * netRatesOfProgress_ddCi();
}

void Kinetics::addThermo(shared_ptr<ThermoPhase> thermo)
{
    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo->nDim() <= m_mindim) {
        if (!m_thermo.empty()) {
            throw CanteraError("Kinetics::addThermo",
                "The reacting (lowest dimensional) phase must be added first.");
        }
        m_mindim = thermo->nDim();
    }

    m_thermo.push_back(thermo);
    m_phaseindex[m_thermo.back()->name()] = nPhases();
    resizeSpecies();
}

void Kinetics::setParameters(const AnyMap& phaseNode) {
    skipUndeclaredThirdBodies(phaseNode.getBool("skip-undeclared-third-bodies", false));
    setExplicitThirdBodyDuplicateHandling(
        phaseNode.getString("explicit-third-body-duplicates", "warn"));

    if (phaseNode.hasKey("rate-multipliers")) {
        const auto& defaultMultipliers = phaseNode["rate-multipliers"];
        for (auto& [key, val] : defaultMultipliers) {
            if (key == "default") {
                m_defaultPerturb[-1] = val.asDouble();
            } else {
                m_defaultPerturb[stoi(key)] = val.asDouble();
            }
        }
    }
}

AnyMap Kinetics::parameters() const
{
    AnyMap out;
    string name = KineticsFactory::factory()->canonicalize(kineticsType());
    if (name != "none") {
        out["kinetics"] = name;
        if (nReactions() == 0) {
            out["reactions"] = "none";
        }
        if (m_hasUndeclaredThirdBodies) {
            out["skip-undeclared-third-bodies"] = true;
        }
        if (m_explicit_third_body_duplicates == "error") {
            // "warn" is the default, and does not need to be added. "mark-duplicate"
            // and "modify-efficiency" do not need to be propagated here as their
            // effects are already applied to the corresponding reactions.
            out["explicit-third-body-duplicates"] = "error";
        }
        map<double, int> multipliers;
        for (auto m : m_perturb) {
            multipliers[m] += 1;
        }
        if (multipliers[1.0] != nReactions()) {
            int defaultCount = 0;
            double defaultMultiplier = 1.0;
            for (auto& [m, count] : multipliers) {
                if (count > defaultCount) {
                    defaultCount = count;
                    defaultMultiplier = m;
                }
            }
            AnyMap multiplierMap;
            multiplierMap["default"] = defaultMultiplier;
            for (size_t i = 0; i < nReactions(); i++) {
                if (m_perturb[i] != defaultMultiplier) {
                    multiplierMap[to_string(i)] = m_perturb[i];
                }
            }
            out["rate-multipliers"] = multiplierMap;
        }
    }
    return out;
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

bool Kinetics::addReaction(shared_ptr<Reaction> r, bool resize)
{
    r->check();
    r->validate(*this);

    if (m_kk == 0) {
        init();
    }
    resizeSpecies();

    // Check validity of reaction within the context of the Kinetics object
    if (!r->checkSpecies(*this)) {
        // Do not add reaction
        return false;
    }

    // For reactions created outside the context of a Kinetics object, the units
    // of the rate coefficient can't be determined in advance. Do that here.
    if (r->rate_units.factor() == 0) {
        r->rate()->setRateUnits(r->calculateRateCoeffUnits(*this));
    }

    size_t irxn = nReactions(); // index of the new reaction

    // indices of reactant and product species within this Kinetics object
    vector<size_t> rk, pk;

    // Reactant and product stoichiometric coefficients, such that rstoich[i] is
    // the coefficient for species rk[i]
    vector<double> rstoich, pstoich;

    for (const auto& [name, stoich] : r->reactants) {
        rk.push_back(kineticsSpeciesIndex(name));
        rstoich.push_back(stoich);
    }

    for (const auto& [name, stoich] : r->products) {
        pk.push_back(kineticsSpeciesIndex(name));
        pstoich.push_back(stoich);
    }

    // The default order for each reactant is its stoichiometric coefficient,
    // which can be overridden by entries in the Reaction.orders map. rorder[i]
    // is the order for species rk[i].
    vector<double> rorder = rstoich;
    for (const auto& [name, order] : r->orders) {
        size_t k = kineticsSpeciesIndex(name);
        // Find the index of species k within rk
        auto rloc = std::find(rk.begin(), rk.end(), k);
        if (rloc != rk.end()) {
            rorder[rloc - rk.begin()] = order;
        } else {
            // If the reaction order involves a non-reactant species, add an
            // extra term to the reactants with zero stoichiometry so that the
            // stoichiometry manager can be used to compute the global forward
            // reaction rate.
            rk.push_back(k);
            rstoich.push_back(0.0);
            rorder.push_back(order);
        }
    }

    m_reactantStoich.add(irxn, rk, rorder, rstoich);
    // product orders = product stoichiometric coefficients
    m_productStoich.add(irxn, pk, pstoich, pstoich);
    if (r->reversible) {
        m_revindex.push_back(irxn);
        m_revProductStoich.add(irxn, pk, pstoich, pstoich);
    } else {
        m_irrev.push_back(irxn);
    }

    m_reactions.push_back(r);
    m_rfn.push_back(0.0);
    m_delta_gibbs0.push_back(0.0);
    m_rkcn.push_back(0.0);
    m_ropf.push_back(0.0);
    m_ropr.push_back(0.0);
    m_ropnet.push_back(0.0);
    m_perturb.push_back(getValue(m_defaultPerturb, irxn, m_defaultPerturb[-1]));
    m_dH.push_back(0.0);

    if (resize) {
        resizeReactions();
    }

    for (const auto& [id, callback] : m_reactionAddedCallbacks) {
        callback();
    }

    return true;
}

void Kinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    checkReactionIndex(i);
    shared_ptr<Reaction>& rOld = m_reactions[i];

    if (rNew->rate()->type() == "electron-collision-plasma") {
        throw CanteraError("Kinetics::modifyReaction",
            "Type electron-collision-plasma is not supported. "
            "Use the rate object of the reaction to modify the data.");
    }

    if (rNew->type() != rOld->type()) {
        throw CanteraError("Kinetics::modifyReaction",
            "Reaction types are different: {} != {}.",
            rOld->type(), rNew->type());
    }

    if (rNew->rate()->type() != rOld->rate()->type()) {
        throw CanteraError("Kinetics::modifyReaction",
            "ReactionRate types are different: {} != {}.",
            rOld->rate()->type(), rNew->rate()->type());
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
