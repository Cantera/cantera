/**
 *  @file Kinetics.cpp Declarations for the base class for kinetics managers
 *      (see \ref  kineticsmgr and class \link Cantera::Kinetics  Kinetics \endlink).
 *
 *  Kinetics managers calculate rates of progress of species due to
 *  homogeneous or heterogeneous kinetics.
 */
// Copyright 2001-2004  California Institute of Technology

#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/stringUtils.h"

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

Kinetics::Kinetics(const Kinetics& right)
{
    // Call the assignment operator
    *this = right;
}

Kinetics& Kinetics::operator=(const Kinetics& right)
{
    // Check for self assignment.
    if (this == &right) {
        return *this;
    }

    m_reactantStoich = right.m_reactantStoich;
    m_revProductStoich = right.m_revProductStoich;
    m_irrevProductStoich = right.m_irrevProductStoich;
    m_kk = right.m_kk;
    m_perturb = right.m_perturb;
    m_reactions = right.m_reactions;
    m_thermo = right.m_thermo; // DANGER -> shallow pointer copy
    m_start = right.m_start;
    m_phaseindex = right.m_phaseindex;
    m_surfphase = right.m_surfphase;
    m_rxnphase = right.m_rxnphase;
    m_mindim = right.m_mindim;
    m_rgroups = right.m_rgroups;
    m_pgroups = right.m_pgroups;
    m_rfn = right.m_rfn;
    m_rkcn = right.m_rkcn;
    m_ropf = right.m_ropf;
    m_ropr = right.m_ropr;
    m_ropnet = right.m_ropnet;
    m_skipUndeclaredSpecies = right.m_skipUndeclaredSpecies;

    return *this;
}

Kinetics* Kinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    Kinetics* ko = new Kinetics(*this);
    ko->assignShallowPointers(tpVector);
    return ko;
}

int Kinetics::type() const
{
    return 0;
}

void Kinetics::checkReactionIndex(size_t i) const
{
    if (i >= nReactions()) {
        throw IndexError("checkReactionIndex", "reactions", i, nReactions()-1);
    }
}

void Kinetics::checkReactionArraySize(size_t ii) const
{
    if (nReactions() > ii) {
        throw ArraySizeError("checkReactionArraySize", ii, nReactions());
    }
}

void Kinetics::checkPhaseIndex(size_t m) const
{
    if (m >= nPhases()) {
        throw IndexError("checkPhaseIndex", "phase", m, nPhases()-1);
    }
}

void Kinetics::checkPhaseArraySize(size_t mm) const
{
    if (nPhases() > mm) {
        throw ArraySizeError("checkPhaseArraySize", mm, nPhases());
    }
}

void Kinetics::checkSpeciesIndex(size_t k) const
{
    if (k >= m_kk) {
        throw IndexError("checkSpeciesIndex", "species", k, m_kk-1);
    }
}

void Kinetics::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_kk);
    }
}

void Kinetics::assignShallowPointers(const std::vector<thermo_t*> & tpVector)
{
    size_t ns = tpVector.size();
    if (ns != m_thermo.size()) {
        throw CanteraError(" Kinetics::assignShallowPointers",
                           " Number of ThermoPhase objects arent't the same");
    }
    for (size_t i = 0; i < ns; i++) {
        ThermoPhase* ntp = tpVector[i];
        ThermoPhase* otp = m_thermo[i];
        if (ntp->id() != otp->id()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " id() of the ThermoPhase objects isn't the same");
        }
        if (ntp->eosType() != otp->eosType()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " eosType() of the ThermoPhase objects isn't the same");
        }
        if (ntp->nSpecies() != otp->nSpecies()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " Number of ThermoPhase objects isn't the same");
        }
        m_thermo[i] = tpVector[i];
    }
}

std::pair<size_t, size_t> Kinetics::checkDuplicates(bool throw_err) const
{
    //! Map of (key indicating participating species) to reaction numbers
    std::map<size_t, std::vector<size_t> > participants;
    std::vector<std::map<int, double> > net_stoich;

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
        for (size_t m = 0; m < related.size(); m++) {
            Reaction& other = *m_reactions[related[m]];
            if (R.reaction_type != other.reaction_type) {
                continue; // different reaction types
            } else if (R.duplicate && other.duplicate) {
                continue; // marked duplicates
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
                throw CanteraError("installReaction",
                        "Undeclared duplicate reactions detected:\n"
                        "Reaction {}: {}\nReaction {}: {}\n",
                        i+1, other.equation(), m+1, R.equation());
            } else {
                return {i,m};
            }
        }
        participants[key].push_back(i);
    }
    return {npos, npos};
}

double Kinetics::checkDuplicateStoich(std::map<int, double>& r1,
                                      std::map<int, double>& r2) const
{
    auto b = r1.begin(), e = r1.end();
    int k1 = b->first;
    // check for duplicate written in the same direction
    doublereal ratio = 0.0;
    if (r1[k1] && r2[k1]) {
        ratio = r2[k1]/r1[k1];
        ++b;
        bool different = false;
        for (; b != e; ++b) {
            k1 = b->first;
            if (!r1[k1] || !r2[k1] || fabs(r2[k1]/r1[k1] - ratio) > 1.e-8) {
                different = true;
                break;
            }
        }
        if (!different) {
            return ratio;
        }
    }

    // check for duplicate written in the reverse direction
    b = r1.begin();
    k1 = b->first;
    if (r1[k1] == 0.0 || r2[-k1] == 0.0) {
        return 0.0;
    }
    ratio = r2[-k1]/r1[k1];
    ++b;
    for (; b != e; ++b) {
        k1 = b->first;
        if (!r1[k1] || !r2[-k1] || fabs(r2[-k1]/r1[k1] - ratio) > 1.e-8) {
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
        msg = "The following reaction is unbalanced: " + R.equation() + "\n" +
              "  Element    Reactants    Products\n" + msg;
        throw CanteraError("checkReactionBalance", msg);
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
        string id = thermo(n).id();
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
        string id = thermo(n).id();
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
    size_t np = m_thermo.size();
    size_t k;
    string id;
    for (size_t n = 0; n < np; n++) {
        k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return thermo(n);
        }
    }
    throw CanteraError("speciesPhase", "unknown species "+nm);
    return thermo(0);
}

size_t Kinetics::speciesPhaseIndex(size_t k)
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return n;
        }
    }
    throw CanteraError("speciesPhaseIndex", "illegal species index: {}", k);
    return npos;
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
    // if not the first thermo object, set the start position
    // to that of the last object added + the number of its species
    if (m_thermo.size() > 0) {
        m_start.push_back(m_start.back()
                          + m_thermo.back()->nSpecies());
    } else {
        // otherwise start at 0
        m_start.push_back(0);
    }

    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo.nDim() <= m_mindim) {
        m_mindim = thermo.nDim();
        m_rxnphase = nPhases();
    }

    // there should only be one surface phase
    int ptype = -100;
    if (type() == cEdgeKinetics) {
        ptype = cEdge;
    } else if (type() == cInterfaceKinetics) {
        ptype = cSurf;
    }
    if (thermo.eosType() == ptype) {
        m_surfphase = nPhases();
        m_rxnphase = nPhases();
    }
    m_thermo.push_back(&thermo);
    m_phaseindex[m_thermo.back()->id()] = nPhases();
}

void Kinetics::finalize()
{
    m_kk = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        size_t nsp = m_thermo[n]->nSpecies();
        m_kk += nsp;
    }
}

bool Kinetics::addReaction(shared_ptr<Reaction> r)
{
    r->validate();

    // If reaction orders are specified, then this reaction does not follow
    // mass-action kinetics, and is not an elementary reaction. So check that it
    // is not reversible, since computing the reverse rate from thermochemistry
    // only works for elementary reactions.
    if (r->reversible && !r->orders.empty()) {
        throw CanteraError("Kinetics::addReaction", "Reaction orders may only "
            "be given for irreversible reactions");
    }

    // Check for undeclared species
    for (const auto& sp : r->reactants) {
        if (kineticsSpeciesIndex(sp.first) == npos) {
            if (m_skipUndeclaredSpecies) {
                return false;
            } else {
                throw CanteraError("Kinetics::addReaction", "Reaction '" +
                    r->equation() + "' contains the undeclared species '" +
                    sp.first + "'");
            }
        }
    }
    for (const auto& sp : r->products) {
        if (kineticsSpeciesIndex(sp.first) == npos) {
            if (m_skipUndeclaredSpecies) {
                return false;
            } else {
                throw CanteraError("Kinetics::addReaction", "Reaction '" +
                    r->equation() + "' contains the undeclared species '" +
                    sp.first + "'");
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
        size_t k = kineticsSpeciesIndex(sp.first);
        rk.push_back(k);
        rstoich.push_back(sp.second);
    }

    for (const auto& sp : r->products) {
        size_t k = kineticsSpeciesIndex(sp.first);
        pk.push_back(k);
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
}

shared_ptr<Reaction> Kinetics::reaction(size_t i)
{
    checkReactionIndex(i);
    return m_reactions[i];
}

}
