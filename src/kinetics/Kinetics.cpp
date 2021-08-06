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
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/StoichManager.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"
#include <unordered_set>
#include <boost/algorithm/string/join.hpp>

using namespace std;

namespace Cantera
{
//! Extract components describing sparse matrix
size_t sparseComponents(const Eigen::SparseMatrix<double>& mat,
    std::vector<std::pair<int, int>>& indices, vector_fp& values)
{
    size_t count = 0;
    for (int i = 0; i < mat.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
            if (count < indices.size()) {
                indices[count] = std::make_pair(it.row(), it.col());
            }
            if (count < values.size()) {
                values[count] = it.value();
            }
            count++;
        }
    }
    if (count > indices.size() || count > values.size()) {
        throw CanteraError("sparseComponents",
            "Output vectors have insufficient length. Required size is {}, "
            "while provided lengths are:\nindices.size()={} and "
            "values.size()={}.", count, indices.size(), values.size());
    }
    return count;
}

Kinetics::Kinetics() :
    m_finalized(false),
    m_kk(0),
    m_thermo(0),
    m_surfphase(npos),
    m_rxnphase(npos),
    m_mindim(4),
    m_skipUndeclaredSpecies(false),
    m_skipUndeclaredThirdBodies(false)
{
    m_reactantStoich = std::unique_ptr<StoichManagerN>(new StoichManagerN());
    m_productStoich = std::unique_ptr<StoichManagerN>(new StoichManagerN());
    m_revProductStoich = std::unique_ptr<StoichManagerN>(new StoichManagerN());
}

Kinetics::~Kinetics() {}

void Kinetics::checkReactionIndex(size_t i) const
{
    if (i >= nReactions()) {
        throw IndexError("Kinetics::checkReactionIndex", "reactions", i,
                         nReactions()-1);
    }
}

void Kinetics::finalizeSetup()
{
    size_t nRxn = nReactions();

    // Stoichiometry managers
    m_reactantStoich->finalizeSetup(m_kk, nRxn);
    m_productStoich->finalizeSetup(m_kk, nRxn);
    m_revProductStoich->finalizeSetup(m_kk, nRxn);

    // products are created for positive net rate of progress
    m_stoichMatrix = m_productStoich->stoichCoeffs();
    // reactants are destroyed for positive net rate of progress
    m_stoichMatrix -= m_reactantStoich->stoichCoeffs();

    m_finalized = true;
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
            } else if (R.type() != other.type()) {
                continue; // different reaction types
            }
            doublereal c = checkDuplicateStoich(net_stoich[i], net_stoich[m]);
            if (c == 0) {
                continue; // stoichiometries differ (not by a multiple)
            } else if (c < 0.0 && !R.reversible && !other.reversible) {
                continue; // irreversible reactions in opposite directions
            } else if (R.type() == "falloff" ||
                       R.type() == "chemically-activated") {
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
            } else if (R.type() == "three-body") {
                ThirdBody& tb1 = *(dynamic_cast<ThreeBodyReaction3&>(R).thirdBody());
                ThirdBody& tb2 = *(dynamic_cast<ThreeBodyReaction3&>(other).thirdBody());
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
            } else if (R.type() == "three-body-legacy") {
                ThirdBody& tb1 = dynamic_cast<ThreeBodyReaction2&>(R).third_body;
                ThirdBody& tb2 = dynamic_cast<ThreeBodyReaction2&>(other).third_body;
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
    R.checkBalance(*this);
    warn_deprecated("Kinetics::checkReactionBalance",
        "To be removed after Cantera 2.6. Replacable by Reaction::checkBalance.");
}

void Kinetics::selectPhase(const double* data, const ThermoPhase* phase,
                           double* phase_data)
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

ThermoPhase& Kinetics::speciesPhase(const std::string& nm)
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        size_t k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return thermo(n);
        }
    }
    throw CanteraError("Kinetics::speciesPhase", "unknown species '{}'", nm);
}

const ThermoPhase& Kinetics::speciesPhase(const std::string& nm) const
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

double Kinetics::reactantStoichCoeff(size_t kSpec, size_t irxn)
{
    if (kSpec >= m_kk || irxn >= nReactions()) {
        throw CanteraError("Kinetics::reactantStoichCoeff",
            "Matrix index ({},{}) is out of range; size is ({},{}).",
            kSpec, irxn, m_kk, nReactions());
    }
    if (!m_finalized) {
        finalizeSetup();
    }
    return m_reactantStoich->stoichCoeffs().coeff(kSpec, irxn);
}

size_t Kinetics::reactantStoichCoeffs(
    std::vector<std::pair<int, int>>& indices, vector_fp& coeffs)
{
    if (!m_finalized) {
        finalizeSetup();
    }
    return m_reactantStoich->sparseStoichCoeffs(indices, coeffs);
}

double Kinetics::productStoichCoeff(size_t kSpec, size_t irxn, bool irreversible)
{
    if (kSpec >= m_kk || irxn >= nReactions()) {
        throw CanteraError("Kinetics::productStoichCoeff",
            "Matrix index ({},{}) is out of range; size is ({},{}).",
            kSpec, irxn, m_kk, nReactions());
    }
    if (!m_finalized) {
        finalizeSetup();
    }
    if (irreversible) {
        return m_productStoich->stoichCoeffs().coeff(kSpec, irxn);
    }
    return m_revProductStoich->stoichCoeffs().coeff(kSpec, irxn);
}

size_t Kinetics::productStoichCoeffs(
    std::vector<std::pair<int, int>>& indices, vector_fp& coeffs, bool irreversible)
{
    if (!m_finalized) {
        finalizeSetup();
    }
    if (irreversible) {
        return m_productStoich->sparseStoichCoeffs(indices, coeffs);
    }
    return m_revProductStoich->sparseStoichCoeffs(indices, coeffs);
}

int Kinetics::reactionType(size_t i) const {
    warn_deprecated("Kinetics::reactionType",
        "To be changed after Cantera 2.6. "
        "Return string instead of magic number; use "
        "Kinetics::reactionTypeStr during transition.");
    return m_reactions[i]->reaction_type;
}

std::string Kinetics::reactionTypeStr(size_t i) const {
    return m_reactions[i]->type();
}

std::string Kinetics::reactionString(size_t i) const
{
    return m_reactions[i]->equation();
}

//! Returns a string containing the reactants side of the reaction equation.
std::string Kinetics::reactantString(size_t i) const
{
    return m_reactions[i]->reactantString();
}

//! Returns a string containing the products side of the reaction equation.
std::string Kinetics::productString(size_t i) const
{
    return m_reactions[i]->productString();
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
    m_productStoich->incrementReactions(prop, deltaProp);
    // reactants subtract
    m_reactantStoich->decrementReactions(prop, deltaProp);
}

void Kinetics::getRevReactionDelta(const double* prop, double* deltaProp)
{
    fill(deltaProp, deltaProp + nReactions(), 0.0);
    // products add
    m_revProductStoich->incrementReactions(prop, deltaProp);
    // reactants subtract
    m_reactantStoich->decrementReactions(prop, deltaProp);
}

void Kinetics::getCreationRates(doublereal* cdot)
{
    updateROP();

    // zero out the output array
    fill(cdot, cdot + m_kk, 0.0);

    // the forward direction creates product species
    m_productStoich->incrementSpecies(m_ropf.data(), cdot);

    // the reverse direction creates reactant species
    m_reactantStoich->incrementSpecies(m_ropr.data(), cdot);
}

Eigen::SparseMatrix<double> Kinetics::getCreationRateSpeciesDerivatives(
    bool thirdbodies)
{
    Eigen::SparseMatrix<double> jac;
    // the forward direction creates product species
    jac = m_productStoich->stoichCoeffs()
        * getFwdRopSpeciesDerivatives(thirdbodies);
    // the reverse direction creates reactant species
    jac += m_reactantStoich->stoichCoeffs()
        * getRevRopSpeciesDerivatives(thirdbodies);
    return jac;
}

Eigen::VectorXd Kinetics::getCreationRateTemperatureDerivatives()
{
    Eigen::VectorXd jac;
    // the forward direction creates product species
    jac = m_productStoich->stoichCoeffs() * getFwdRopTemperatureDerivatives();
    // the reverse direction creates reactant species
    jac += m_reactantStoich->stoichCoeffs() * getRevRopTemperatureDerivatives();
    return jac;
}

void Kinetics::getDestructionRates(doublereal* ddot)
{
    updateROP();

    fill(ddot, ddot + m_kk, 0.0);
    // the reverse direction destroys products in reversible reactions
    m_revProductStoich->incrementSpecies(m_ropr.data(), ddot);
    // the forward direction destroys reactants
    m_reactantStoich->incrementSpecies(m_ropf.data(), ddot);
}

Eigen::SparseMatrix<double> Kinetics::getDestructionRateSpeciesDerivatives(
    bool thirdbodies)
{
    Eigen::SparseMatrix<double> jac;
    // the reverse direction destroys products in reversible reactions
    jac = m_revProductStoich->stoichCoeffs()
        * getRevRopSpeciesDerivatives(thirdbodies);
    // the forward direction destroys reactants
    jac += m_reactantStoich->stoichCoeffs()
        * getFwdRopSpeciesDerivatives(thirdbodies);
    return jac;
}

Eigen::VectorXd Kinetics::getDestructionRateTemperatureDerivatives()
{
    Eigen::VectorXd jac;
    // the reverse direction destroys products in reversible reactions
    jac = m_revProductStoich->stoichCoeffs() * getRevRopTemperatureDerivatives();
    // the forward direction destroys reactants
    jac += m_reactantStoich->stoichCoeffs() * getFwdRopTemperatureDerivatives();
    return jac;
}

void Kinetics::getNetProductionRates(doublereal* net)
{
    updateROP();

    fill(net, net + m_kk, 0.0);
    // products are created for positive net rate of progress
    m_productStoich->incrementSpecies(m_ropnet.data(), net);
    // reactants are destroyed for positive net rate of progress
    m_reactantStoich->decrementSpecies(m_ropnet.data(), net);
}

Eigen::SparseMatrix<double> Kinetics::getNetProductionRateSpeciesDerivatives(
    bool thirdbodies)
{
    return m_stoichMatrix * (
        getFwdRopSpeciesDerivatives(thirdbodies)
        - getRevRopSpeciesDerivatives(thirdbodies));
}

Eigen::VectorXd Kinetics::getNetProductionRateTemperatureDerivatives()
{
    return m_stoichMatrix * (
        getFwdRopTemperatureDerivatives() - getRevRopTemperatureDerivatives());
}

size_t Kinetics::getRopSpeciesDerivatives(
    std::vector<std::pair<int, int>>& indices, vector_fp& values,
    bool forward, bool reverse, bool thirdbodies)
{
    Eigen::SparseMatrix<double> ret;
    if (forward && reverse) {
        ret = getFwdRopSpeciesDerivatives(thirdbodies)
            - getRevRopSpeciesDerivatives(thirdbodies);
    } else if (forward) {
        ret = getFwdRopSpeciesDerivatives(thirdbodies);
    } else if (reverse) {
        ret = getRevRopSpeciesDerivatives(thirdbodies);
    }
    return sparseComponents(ret, indices, values);
}

size_t Kinetics::getProductionRateSpeciesDerivatives(
    std::vector<std::pair<int, int>>& indices, vector_fp& values,
    bool creation, bool destruction, bool thirdbodies)
{
    Eigen::SparseMatrix<double> ret;
    if (creation && destruction) {
        ret = getNetProductionRateSpeciesDerivatives(thirdbodies);
    } else if (creation) {
        ret = getCreationRateSpeciesDerivatives(thirdbodies);
    } else if (destruction) {
        ret = getDestructionRateSpeciesDerivatives(thirdbodies);
    }
    return sparseComponents(ret, indices, values);
}

void Kinetics::getRopTemperatureDerivatives(
    vector_fp& values, bool forward, bool reverse)
{
    MappedVector mapped(values.data(), nReactions());
    if (forward && reverse) {
        mapped = getFwdRopTemperatureDerivatives() - getRevRopTemperatureDerivatives();
    } else if (forward) {
        mapped = getFwdRopTemperatureDerivatives();
    } else if (reverse) {
        mapped = getRevRopTemperatureDerivatives();
    } else {
        mapped.setZero();
    }
}

void Kinetics::getProductionRateTemperatureDerivatives(
    vector_fp& values, bool creation, bool destruction)
{
    MappedVector mapped(values.data(), nReactions());
    if (creation && destruction) {
        mapped = getNetProductionRateTemperatureDerivatives();
    } else if (creation) {
        mapped = getCreationRateTemperatureDerivatives();
    } else if (destruction) {
        mapped = getDestructionRateTemperatureDerivatives();
    } else {
        mapped.setZero();
    }
}

void Kinetics::addPhase(ThermoPhase& thermo)
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

AnyMap Kinetics::parameters()
{
    AnyMap out;
    string name = KineticsFactory::factory()->canonicalize(kineticsType());
    if (name != "none") {
        out["kinetics"] = name;
        if (nReactions() == 0) {
            out["reactions"] = "none";
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

bool Kinetics::addReaction(shared_ptr<Reaction> r)
{
    r->validate();
    if (m_kk == 0) {
        init();
    }
    resizeSpecies();
    m_finalized = false;

    // Check validity of reaction within the context of the Kinetics object
    if (!r->checkSpecies(*this)) {
        // Do not add reaction
        return false;
    }

    // For reactions created outside the context of a Kinetics object, the units
    // of the rate coefficient can't be determined in advance. Do that here.
    if (r->rate_units.factor() == 0) {
        r->calculateRateCoeffUnits(*this);
    }

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

    m_reactantStoich->add(irxn, rk, rorder, rstoich);
    // product orders = product stoichiometric coefficients
    m_productStoich->add(irxn, pk, pstoich, pstoich);
    if (r->reversible) {
        m_revProductStoich->add(irxn, pk, pstoich, pstoich);
    }

    m_reactions.push_back(r);
    m_rfn.push_back(0.0);
    m_rkcn.push_back(0.0);
    m_ropf.push_back(0.0);
    m_ropr.push_back(0.0);
    m_ropnet.push_back(0.0);
    m_perturb.push_back(1.0);
    m_dH.push_back(0.0);
    return true;
}

void Kinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    checkReactionIndex(i);
    shared_ptr<Reaction>& rOld = m_reactions[i];
    if (rNew->type() != rOld->type()) {
        throw CanteraError("Kinetics::modifyReaction",
            "Reaction types are different: {} != {}.",
            rOld->type(), rNew->type());
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
