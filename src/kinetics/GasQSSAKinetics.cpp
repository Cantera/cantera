/**
 *  @file GasQSSAKinetics.cpp Homogeneous kinetics in ideal gases w/ QSSA
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/GasQSSAKinetics.h"

using namespace std;

namespace Cantera
{
GasQSSAKinetics::GasQSSAKinetics(thermo_t *thermo) :
    GasKinetics(thermo), m_QSS_init(false), m_QSS_ok(false)
{

}

GasQSSAKinetics::~GasQSSAKinetics()
{

}

void GasQSSAKinetics::getEquilibriumConstants(doublereal* kc)
{
    update_rates_T();
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_grt.data() + m_start[n]);
    }
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo(0).RT();
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-m_rkcn[i]*rrt + m_dn[i]*m_logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void GasQSSAKinetics::getFwdRateConstants(doublereal* kfwd)
{
    GasKinetics::getFwdRateConstants(kfwd);
}

void GasQSSAKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the all of the phases used
    // in the kinetics mechanism
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaG);
}

void GasQSSAKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void GasQSSAKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in all of the phases
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void GasQSSAKinetics::getDeltaSSGibbs(doublereal* deltaG)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity We define these here as the
    // chemical potentials of the pure species at the temperature and pressure
    // of the solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaG);
}

void GasQSSAKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    // Get the standard state enthalpies of the species. This is the array of
    // chemical potentials at unit activity We define these here as the
    // enthalpies of the pure species at the temperature and pressure of the
    // solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEnthalpy_RT(m_grt.data() + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= thermo(0).RT();
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void GasQSSAKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    // Get the standard state entropy of the species. We define these here as
    // the entropies of the pure species at the temperature and pressure of the
    // solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEntropy_R(m_grt.data() + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= GasConstant;
    }

    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void GasQSSAKinetics::init()
{
    GasKinetics::init();
    // GasQSSAKinetics have two phases (ideal_gas, qssa_gas)
    assert(nPhases() == 2);
    m_nSpeciesQSS = thermo(1).nSpecies();
    m_rodf_qss.resize(m_nSpeciesQSS);
    m_rodr_qss.resize(m_nSpeciesQSS);
    m_ropf_noqss.resize(m_nSpeciesQSS);
    m_ropr_noqss.resize(m_nSpeciesQSS);
    m_ropf_qss = vector<vector<vector<size_t>>>(
            m_nSpeciesQSS, vector<vector<size_t>>(m_nSpeciesQSS));
    m_ropr_qss = vector<vector<vector<size_t>>>(
            m_nSpeciesQSS, vector<vector<size_t>>(m_nSpeciesQSS));
}

bool GasQSSAKinetics::addReaction(shared_ptr<Reaction> r)
{
    // operations common to all reaction types
    bool added = GasKinetics::addReaction(r);
    if (!added) {
        return false;
    }
    return addReactionQSS(r);
}

void GasQSSAKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    GasKinetics::modifyReaction(i, rNew);
}

void GasQSSAKinetics::updateROP()
{
    update_rates_C();
    update_rates_T();
    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(m_ropf.data(), concm_3b_values.data());
    }

    if (m_falloff_high_rates.nReactions()) {
        processFalloffReactions();
    }

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

    //

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    for (size_t i = 0; i < m_rfn.size(); i++) {
        AssertFinite(m_rfn[i], "GasKinetics::updateROP",
                     "m_rfn[{}] is not finite.", i);
        AssertFinite(m_ropf[i], "GasKinetics::updateROP",
                     "m_ropf[{}] is not finite.", i);
        AssertFinite(m_ropr[i], "GasKinetics::updateROP",
                     "m_ropr[{}] is not finite.", i);
    }
    m_ROP_ok = true;
}

void GasQSSAKinetics::update_rates_T()
{
    doublereal T = thermo(0).temperature();
    doublereal P = thermo(0).pressure();
    m_logStandConc = log(thermo(0).standardConcentration());
    doublereal logT = log(T);

    if (T != m_temp) {
        if (!m_rfn.empty()) {
            m_rates.update(T, logT, m_rfn.data());
        }

        if (!m_rfn_low.empty()) {
            m_falloff_low_rates.update(T, logT, m_rfn_low.data());
            m_falloff_high_rates.update(T, logT, m_rfn_high.data());
        }
        if (!falloff_work.empty()) {
            m_falloffn.updateTemp(T, falloff_work.data());
        }
        updateKc();
        m_ROP_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void GasQSSAKinetics::update_rates_C()
{
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getActivityConcentrations(m_conc.data() + m_start[n]);
    }
    doublereal ctot = thermo(0).molarDensity();

    // 3-body reactions
    if (!concm_3b_values.empty()) {
        m_3b_concm.update(m_conc, ctot, concm_3b_values.data());
    }

    // Falloff reactions
    if (!concm_falloff_values.empty()) {
        m_falloff_concm.update(m_conc, ctot, concm_falloff_values.data());
    }

    // P-log reactions
    if (m_plog_rates.nReactions()) {
        double logP = log(thermo(0).pressure());
        m_plog_rates.update_C(&logP);
    }

    // Chebyshev reactions
    if (m_cheb_rates.nReactions()) {
        double log10P = log10(thermo(0).pressure());
        m_cheb_rates.update_C(&log10P);
    }

    m_ROP_ok = false;
}

void GasQSSAKinetics::updateKc()
{
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_grt.data() + m_start[n]);
    }
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = std::min(exp(m_rkcn[irxn]*rrt - m_dn[irxn]*m_logStandConc),
                                BigNumber);
    }

    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

bool GasQSSAKinetics::addReactionQSS(shared_ptr<Reaction> r)
{
    vector<size_t> qss_rts, qss_pds;
    // check if reactants are QSS
    for (const auto& sp : r->reactants) {
        size_t k = thermo(1).speciesIndex(sp.first);
        if (k != npos) qss_rts.push_back(k);
    }
    // check if products are QSS
    for (const auto& sp : r->products) {
        size_t k = thermo(1).speciesIndex(sp.first);
        if (k != npos) qss_pds.push_back(k);
    }

    // return if no QSS on either side
    if (!qss_rts.size() && !qss_pds.size()) {
        return true;
    }
    // warn if more than one QSS on one side
    if (qss_rts.size() > 1 || qss_pds.size() > 1)  {
        writelog(
            "WARNING: more than one QSS species on"
            " one side of reaction: {}.\n", r->equation());
    }
    // fill m_rodf_qss
    for (const auto rt : qss_rts) {
        m_rodf_qss[rt].push_back(nReactions()- 1);
    }
    // fill m_rodr_qss
    if (r->reversible) {
        for (const auto pd : qss_pds) {
            m_rodr_qss[pd].push_back(nReactions()- 1);
        }
    }
    // fill m_ropf_noqss
    if (!qss_rts.size()) { // if no qss reactants
        for (const auto pd : qss_pds) {
            m_ropf_noqss[pd].push_back(nReactions()- 1);
        }
    }
    // fill m_ropr_noqss
    if (!qss_pds.size() && r->reversible) { // if no qss products
        for (const auto rt : qss_rts) {
            m_ropr_noqss[rt].push_back(nReactions()- 1);
        }
    }
    // fill m_ropf_qss and m_ropr_qss
    if (qss_rts.size() && qss_pds.size()) {
        for (const auto rt : qss_rts) {
            for (const auto pd : qss_pds) {
                m_ropf_qss[rt][pd].push_back(nReactions()- 1);
                if (r->reversible)
                    m_ropr_qss[pd][rt].push_back(nReactions()- 1);
            }
        }
    }
    return true;
}

void GasQSSAKinetics::init_QSS()
{
    
}

}
