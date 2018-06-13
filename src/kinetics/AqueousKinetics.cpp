/**
 *  @file AqueousKinetics.cpp
 *
 * Homogeneous kinetics in an aqueous phase, either condensed or dilute in salts
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/AqueousKinetics.h"
#include "cantera/kinetics/Reaction.h"

using namespace std;

namespace Cantera
{

AqueousKinetics::AqueousKinetics(thermo_t* thermo) :
    BulkKinetics(thermo)
{
    warn_deprecated("Class AqueousKinetics", "To be removed after Cantera 2.4");
}

void AqueousKinetics::_update_rates_T()
{
    doublereal T = thermo().temperature();
    m_rates.update(T, log(T), m_rfn.data());

    m_temp = T;
    updateKc();
    m_ROP_ok = false;
}

void AqueousKinetics::_update_rates_C()
{
    thermo().getActivityConcentrations(m_conc.data());
    m_ROP_ok = false;
}

void AqueousKinetics::updateKc()
{
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);
    for (size_t k = 0; k < thermo().nSpecies(); k++) {
        doublereal logStandConc_k = thermo().logStandardConc(k);
        m_grt[k] -= GasConstant * m_temp * logStandConc_k;
    }

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
    }

    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void AqueousKinetics::getEquilibriumConstants(doublereal* kc)
{
    _update_rates_T();

    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);
    for (size_t k = 0; k < thermo().nSpecies(); k++) {
        doublereal logStandConc_k = thermo().logStandardConc(k);
        m_grt[k] -= GasConstant * m_temp * logStandConc_k;
    }

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-m_rkcn[i]*rrt);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void AqueousKinetics::updateROP()
{
    _update_rates_T();
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    m_ropr = m_ropf;

    // for reverse rates computed from thermochemistry, multiply the forward
    // rates copied into m_ropr by the reciprocals of the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich.multiply(m_conc.data(), m_ropr.data());

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    m_ROP_ok = true;
}

void AqueousKinetics::getFwdRateConstants(doublereal* kfwd)
{
    _update_rates_T();
    _update_rates_C();

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    for (size_t i = 0; i < nReactions(); i++) {
        kfwd[i] = m_ropf[i];
    }
}

bool AqueousKinetics::addReaction(shared_ptr<Reaction> r)
{
    bool added = BulkKinetics::addReaction(r);
    if (!added) {
        return false;
    }
    if (r->reaction_type == ELEMENTARY_RXN) {
        addElementaryReaction(dynamic_cast<ElementaryReaction&>(*r));
    } else {
        throw CanteraError("AqueousKinetics::addReaction",
            "Invalid reaction type: {}", r->reaction_type);
    }
    return true;
}

void AqueousKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    BulkKinetics::modifyReaction(i, rNew);
    modifyElementaryReaction(i, dynamic_cast<ElementaryReaction&>(*rNew));

    // invalidate all cached data
    m_ROP_ok = false;
    m_temp += 0.1234;
}

}
