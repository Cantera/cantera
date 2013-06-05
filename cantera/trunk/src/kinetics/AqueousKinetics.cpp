/**
 *  @file AqueousKinetics.cpp
 *
 * Homogeneous kinetics in an aqueous phase, either condensed
 * or dilute in salts
 *
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/kinetics/AqueousKinetics.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

AqueousKinetics::AqueousKinetics(thermo_t* thermo) :
    Kinetics(),
    m_nfall(0),
    m_nirrev(0),
    m_nrev(0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_finalized(false)
{
    warn_deprecated("AqueousKinetics",
                    "Unfinished implementation of this class will be removed.");
    if (thermo != 0) {
        addPhase(*thermo);
    }
}

AqueousKinetics::AqueousKinetics(const AqueousKinetics& right) :
    Kinetics(),
    m_nfall(0),
    m_nirrev(0),
    m_nrev(0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_finalized(false)
{
    *this = right;
}

AqueousKinetics& AqueousKinetics::operator=(const AqueousKinetics& right)
{
    if (this == &right) {
        return *this;
    }

    Kinetics::operator=(right);

    m_nfall = right.m_nfall;
    m_rates = right.m_rates;
    m_index = right.m_index;
    m_irrev = right.m_irrev;

    m_rxnstoich = right.m_rxnstoich;

    m_fwdOrder = right.m_fwdOrder;
    m_nirrev = right.m_nirrev;
    m_nrev = right.m_nrev;
    m_rgroups = right.m_rgroups;
    m_pgroups = right.m_pgroups;
    m_rxntype = right.m_rxntype;
    m_rrxn = right.m_rrxn;
    m_prxn = right.m_prxn;
    m_dn = right.m_dn;
    m_revindex = right.m_revindex;
    m_rxneqn = right.m_rxneqn;

    m_ropf = right.m_ropf;
    m_ropr = right.m_ropr;
    m_ropnet = right.m_ropnet;
    m_ROP_ok  = right.m_ROP_ok;
    m_temp = right.m_temp;
    m_rfn = right.m_rfn;
    m_rkcn = right.m_rkcn;

    m_conc = right.m_conc;
    m_grt = right.m_grt;
    m_finalized = right.m_finalized;

    throw CanteraError("GasKinetics::operator=()",
                       "Unfinished implementation");

    return *this;

}

Kinetics* AqueousKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    AqueousKinetics* gK = new AqueousKinetics(*this);
    gK->assignShallowPointers(tpVector);
    return gK;
}

void AqueousKinetics::
update_T() {}

void AqueousKinetics::
update_C() {}

void AqueousKinetics::_update_rates_T()
{
    doublereal T = thermo().temperature();
    doublereal logT = log(T);
    m_rates.update(T, logT, &m_rfn[0]);

    m_temp = T;
    updateKc();
    m_ROP_ok = false;
}

void AqueousKinetics::
_update_rates_C()
{
    thermo().getActivityConcentrations(&m_conc[0]);

    m_ROP_ok = false;
}

void AqueousKinetics::updateKc()
{
    doublereal rt = GasConstant * m_temp;

    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);
    for (size_t k = 0; k < thermo().nSpecies(); k++) {
        doublereal logStandConc_k = thermo().logStandardConc(k);
        m_grt[k] -= rt * logStandConc_k;
    }

    // compute Delta G^0 for all reversible reactions
    m_rxnstoich.getRevReactionDelta(m_ii, &m_grt[0], &m_rkcn[0]);

    //doublereal logStandConc = m_kdata->m_logStandConc;
    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_nrev; i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
    }

    for (size_t i = 0; i != m_nirrev; ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void AqueousKinetics::getEquilibriumConstants(doublereal* kc)
{
    _update_rates_T();

    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);
    doublereal rt = GasConstant * m_temp;
    for (size_t k = 0; k < thermo().nSpecies(); k++) {
        doublereal logStandConc_k = thermo().logStandardConc(k);
        m_grt[k] -= rt * logStandConc_k;
    }

    // compute Delta G^0 for all reactions
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], &m_rkcn[0]);

    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-m_rkcn[i]*rrt);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void AqueousKinetics::getDeltaGibbs(doublereal* deltaG)
{
    /*
     * Get the chemical potentials of the species in the
     * ideal gas solution.
     */
    thermo().getChemPotentials(&m_grt[0]);
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaG);
}

void AqueousKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    /*
     * Get the partial molar enthalpy of all species in the
     * ideal gas.
     */
    thermo().getPartialMolarEnthalpies(&m_grt[0]);
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaH);
}

void AqueousKinetics::getDeltaEntropy(doublereal* deltaS)
{
    /*
     * Get the partial molar entropy of all species in the
     * solid solution.
     */
    thermo().getPartialMolarEntropies(&m_grt[0]);
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaS);
}

void AqueousKinetics::getDeltaSSGibbs(doublereal* deltaG)
{
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    thermo().getStandardChemPotentials(&m_grt[0]);
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaG);
}

void AqueousKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    /*
     *  Get the standard state enthalpies of the species.
     *  This is the array of chemical potentials at unit activity
     *  We define these here as the enthalpies of the pure
     *  species at the temperature and pressure of the solution.
     */
    thermo().getEnthalpy_RT(&m_grt[0]);
    doublereal RT = thermo().temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= RT;
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaH);
}

void AqueousKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    /*
     *  Get the standard state entropy of the species.
     *  We define these here as the entropies of the pure
     *  species at the temperature and pressure of the solution.
     */
    thermo().getEntropy_R(&m_grt[0]);
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= R;
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaS);
}

void AqueousKinetics::updateROP()
{
    _update_rates_T();
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // for reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by concentration products
    m_rxnstoich.multiplyReactants(&m_conc[0], &m_ropf[0]);
    //m_reactantStoich.multiply(m_conc.begin(), ropf.begin());

    // for reversible reactions, multiply ropr by concentration
    // products
    m_rxnstoich.multiplyRevProducts(&m_conc[0], &m_ropr[0]);
    //m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

    for (size_t j = 0; j != m_ii; ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    m_ROP_ok = true;
}

void AqueousKinetics::
getFwdRateConstants(doublereal* kfwd)
{
    _update_rates_T();
    _update_rates_C();

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    for (size_t i = 0; i < m_ii; i++) {
        kfwd[i] = m_ropf[i];
    }
}

void AqueousKinetics::
getRevRateConstants(doublereal* krev, bool doIrreversible)
{
    /*
     * go get the forward rate constants. -> note, we don't
     * really care about speed or redundancy in these
     * informational routines.
     */
    getFwdRateConstants(krev);

    if (doIrreversible) {
        getEquilibriumConstants(&m_ropnet[0]);
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] /=  m_ropnet[i];
        }
    } else {
        /*
         * m_rkcn[] is zero for irreversible reactions
         */
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

void AqueousKinetics::addReaction(ReactionData& r)
{
    if (r.reactionType == ELEMENTARY_RXN) {
        addElementaryReaction(r);
    }

    // operations common to all reaction types
    installReagents(r);
    installGroups(reactionNumber(), r.rgroups, r.pgroups);
    incrementRxnCount();
    m_rxneqn.push_back(r.equation);
}

void AqueousKinetics::addElementaryReaction(ReactionData& r)
{
    size_t iloc;

    // install rate coeff calculator
    iloc = m_rates.install(reactionNumber(), r);

    // add constant term to rate coeff value vector
    m_rfn.push_back(r.rateCoeffParameters[0]);

    // forward rxn order equals number of reactants
    m_fwdOrder.push_back(r.reactants.size());
    registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);
}

void AqueousKinetics::installReagents(const ReactionData& r)
{

    m_ropf.push_back(0.0);     // extend by one for new rxn
    m_ropr.push_back(0.0);
    m_ropnet.push_back(0.0);
    size_t n, ns, m;
    doublereal nsFlt;
    doublereal reactantGlobalOrder = 0.0;
    doublereal productGlobalOrder  = 0.0;
    size_t rnum = reactionNumber();

    std::vector<size_t> rk;
    size_t nr = r.reactants.size();
    for (n = 0; n < nr; n++) {
        nsFlt = r.rstoich[n];
        reactantGlobalOrder += nsFlt;
        ns = (size_t) nsFlt;
        if ((doublereal) ns != nsFlt) {
            if (ns < 1) {
                ns = 1;
            }
        }
        if (r.rstoich[n] != 0.0) {
            m_rrxn[r.reactants[n]][rnum] += r.rstoich[n];
        }
        for (m = 0; m < ns; m++) {
            rk.push_back(r.reactants[n]);
        }
    }
    m_reactants.push_back(rk);

    std::vector<size_t> pk;
    size_t np = r.products.size();
    for (n = 0; n < np; n++) {
        nsFlt = r.pstoich[n];
        productGlobalOrder += nsFlt;
        ns = (size_t) nsFlt;
        if ((double) ns != nsFlt) {
            if (ns < 1) {
                ns = 1;
            }
        }
        if (r.pstoich[n] != 0.0) {
            m_prxn[r.products[n]][rnum] += r.pstoich[n];
        }
        for (m = 0; m < ns; m++) {
            pk.push_back(r.products[n]);
        }
    }
    m_products.push_back(pk);

    m_rkcn.push_back(0.0);

    m_rxnstoich.add(reactionNumber(), r);

    if (r.reversible) {
        m_dn.push_back(productGlobalOrder - reactantGlobalOrder);
        m_revindex.push_back(reactionNumber());
        m_nrev++;
    } else {
        m_dn.push_back(productGlobalOrder - reactantGlobalOrder);
        m_irrev.push_back(reactionNumber());
        m_nirrev++;
    }
}

void AqueousKinetics::installGroups(size_t irxn,
                                    const vector<grouplist_t>& r,
                                    const vector<grouplist_t>& p)
{
    if (!r.empty()) {
        writelog("installing groups for reaction "+int2str(reactionNumber()));
        m_rgroups[reactionNumber()] = r;
        m_pgroups[reactionNumber()] = p;
    }
}

void AqueousKinetics::init()
{
    m_kk = thermo().nSpecies();
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_conc.resize(m_kk);
    m_grt.resize(m_kk);
}

void AqueousKinetics::finalize()
{
    if (!m_finalized) {
        m_finalized = true;

        // Guarantee that these arrays can be converted to double* even in the
        // special case where there are no reactions defined.
        if (!m_ii) {
            m_perturb.resize(1, 1.0);
            m_ropf.resize(1, 0.0);
            m_ropr.resize(1, 0.0);
            m_ropnet.resize(1, 0.0);
            m_rkcn.resize(1, 0.0);
        }
    }
}

bool AqueousKinetics::ready() const
{
    return m_finalized;
}

}
