/**
 *  @file GasQSSKinetics.cpp Homogeneous kinetics in ideal gases w/ QSSA
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.
// Created by Hao Wu (wuhao@stanford.edu), Feb. 2017

#include "cantera/kinetics/GasQSSKinetics.h"

using namespace std;
using namespace Eigen;

namespace Cantera
{
GasQSSKinetics::GasQSSKinetics(doublereal *_buff_full, thermo_t *thermo) :
    GasKinetics(thermo), m_buff_full(_buff_full),
    m_rel_density_qss(1.e-12), m_QSS_init(false), m_QSS_ok(false)
{

}

GasQSSKinetics::~GasQSSKinetics()
{

}

void GasQSSKinetics::getEquilibriumConstants(doublereal* kc)
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

void GasQSSKinetics::getFwdRateConstants(doublereal* kfwd)
{
    GasKinetics::getFwdRateConstants(kfwd);
}

void GasQSSKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the all of the phases used
    // in the kinetics mechanism
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaG);
}

void GasQSSKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void GasQSSKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in all of the phases
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(m_grt.data() + m_start[n]);
    }
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void GasQSSKinetics::getDeltaSSGibbs(doublereal* deltaG)
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

void GasQSSKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
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

void GasQSSKinetics::getDeltaSSEntropy(doublereal* deltaS)
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

void GasQSSKinetics::getCreationRates(doublereal* cdot)
{
    GasKinetics::getCreationRates(m_buff_full);
    copy(m_buff_full, m_buff_full + m_start[1], cdot);
}

void GasQSSKinetics::getDestructionRates(doublereal* ddot)
{
    GasKinetics::getDestructionRates(m_buff_full);
    copy(m_buff_full, m_buff_full + m_start[1], ddot);
}

void GasQSSKinetics::getNetProductionRates(doublereal* net)
{
    GasKinetics::getNetProductionRates(m_buff_full);
    copy(m_buff_full, m_buff_full + m_start[1], net);
}

void GasQSSKinetics::init()
{
    GasKinetics::init();
    // GasQSSKinetics have two phases (ideal_gas, qssa_gas)
    assert(nPhases() == 2);
    // Resize QSS objects according to m_nSpeciesQSS
    m_nSpeciesQSS = thermo(1).nSpecies();
    m_rodf_qss.resize(m_nSpeciesQSS);
    m_rodr_qss.resize(m_nSpeciesQSS);
    m_ropf_noqss.resize(m_nSpeciesQSS);
    m_ropr_noqss.resize(m_nSpeciesQSS);
    m_ropf_qss_tmp = vector<vector<vector<size_t>>>(
            m_nSpeciesQSS, vector<vector<size_t>>(m_nSpeciesQSS));
    m_ropr_qss_tmp = vector<vector<vector<size_t>>>(
            m_nSpeciesQSS, vector<vector<size_t>>(m_nSpeciesQSS));
    m_rod_qss = VectorXd::Zero(m_nSpeciesQSS);
    m_rop_noqss = VectorXd::Zero(m_nSpeciesQSS);
    m_rop_qss.resize(m_nSpeciesQSS, m_nSpeciesQSS);
    // Make sure m_buff_full is not null
    assert(m_buff_full != nullptr);
}

bool GasQSSKinetics::addReaction(shared_ptr<Reaction> r)
{
    // operations common to all reaction types
    bool added = GasKinetics::addReaction(r);
    if (!added) {
        return false;
    }
    return addReactionQSS(r);
}

void GasQSSKinetics::updateROP()
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

    // make QSS species unit concentration
    Map<VectorXd>(m_conc.data() + m_start[1], m_nSpeciesQSS).setOnes();
    // fill(m_conc.begin() + m_start[1], m_conc.end(), 0.0);

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich.multiply(m_conc.data(), m_ropr.data());

    // calculate concentration of QSS species
    calc_conc_QSS(m_conc.data() + m_start[1]);

    // calculate concentration of QSS species
    update_ROP_QSS(m_conc.data() + m_start[1]);

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

void GasQSSKinetics::update_rates_T()
{
    doublereal T = thermo(0).temperature();
    doublereal P = thermo(0).pressure();

    thermo(1).setTemperature(T);
    thermo(1).setPressure(P);

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
        m_QSS_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
            m_QSS_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
            m_QSS_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void GasQSSKinetics::update_rates_C()
{
    thermo(0).getActivityConcentrations(m_conc.data());
    fill(m_conc.begin() + m_start[1], m_conc.end(), SmallNumber);
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
    m_QSS_ok = false;
}

void GasQSSKinetics::updateKc()
{
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_grt.data() + m_start[n]);
    }
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);
    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_rkcn.data());
    doublereal rrt = 1.0 / thermo(0).RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = std::min(exp(m_rkcn[irxn]*rrt - m_dn[irxn]*m_logStandConc),
                                BigNumber);
    }
    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

bool GasQSSKinetics::addReactionQSS(shared_ptr<Reaction> r)
{
    vector<size_t> qss_rts, qss_pds;
    // check if reactants are QSS
    for (const auto& sp : r->reactants) {
        size_t k = thermo(1).speciesIndex(sp.first);
        if (k == npos) continue;
        for (size_t i = 0; i < sp.second; ++i) qss_rts.push_back(k);
    }
    // check if products are QSS
    for (const auto& sp : r->products) {
        size_t k = thermo(1).speciesIndex(sp.first);
        if (k == npos) continue;
        for (size_t i = 0; i < sp.second; ++i) qss_pds.push_back(k);
    }

    // return if no QSS on either side
    if (!qss_rts.size() && !qss_pds.size()) return true;

    // warn if more than one QSS on one side
    if (qss_rts.size() > 1 || (qss_pds.size() > 1 && r->reversible))  {
        writelog(
            "WARNING: more than one QSS species on"
            " the destruction side of reaction {}: {}.\n",
            nReactions() - 1, r->equation());
    }
    // fill m_rodf_qss
    if (qss_rts.size() > 1) m_rodf_mqss.push_back(nReactions()- 1);
    else if (qss_rts.size() == 1)
        m_rodf_qss[qss_rts[0]].push_back(nReactions()- 1);
    // fill m_rodr_qss
    if (r->reversible) {
        if (qss_pds.size() > 1) m_rodr_mqss.push_back(nReactions()- 1);
        else if (qss_pds.size() == 1)
            m_rodr_qss[qss_pds[0]].push_back(nReactions()- 1);
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
    // fill m_ropf_qss_tmp
    if (qss_rts.size() == 1 && qss_pds.size()) {
        for (const auto pd : qss_pds) {
            m_ropf_qss_tmp[qss_rts[0]][pd].push_back(nReactions()- 1);
        }
    }
    // fill m_ropr_qss_tmp
    if (r->reversible && qss_rts.size() && qss_pds.size() == 1) {
        for (const auto rt : qss_rts) {
            m_ropr_qss_tmp[qss_pds[0]][rt].push_back(nReactions()- 1);
        }
    }
    return true;
}

void GasQSSKinetics::init_QSS()
{
    if (m_QSS_init) return;
    // make sparsity pattern for m_rop_qss
    vector< Triplet<double> > tripletList;
    for (size_t k = 0; k < m_nSpeciesQSS; k++) {
        tripletList.push_back(Triplet<double>(k, k, 1.));
        for (size_t i = 0; i < m_nSpeciesQSS; i++) {
            if (m_ropf_qss_tmp[k][i].size() > 0 ||
                m_ropr_qss_tmp[k][i].size() > 0) {
                tripletList.push_back(Triplet<double>(i, k, 1.));
            }
        }
    }
    m_rop_qss.setFromTriplets(tripletList.begin(), tripletList.end());
    m_rop_qss.makeCompressed();
    // make m_ropf_qss and m_ropr_qss from temp
    m_ropf_qss.clear();
    m_ropr_qss.clear();
    for (int k = 0; k < m_rop_qss.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(m_rop_qss,k); it; ++it) {
            if (it.col() == it.row()) continue;
            char _ifr_qss = 0;
            assert(m_ropf_qss_tmp[it.col()][it.row()].size() > 0 ||
                   m_ropr_qss_tmp[it.col()][it.row()].size() > 0);
            if (m_ropf_qss_tmp[it.col()][it.row()].size() > 0) {
                m_ropf_qss.push_back(m_ropf_qss_tmp[it.col()][it.row()]);
                _ifr_qss |= IROPF;
            }
            if (m_ropr_qss_tmp[it.col()][it.row()].size() > 0) {
                m_ropr_qss.push_back(m_ropr_qss_tmp[it.col()][it.row()]);
                _ifr_qss |= IROPR;
            }
            m_ifr_qss.push_back(_ifr_qss);
        }
    // clear m_ropf_qss_tmp and m_ropr_qss_tmp to save memory
    m_ropf_qss_tmp.clear();
    m_ropr_qss_tmp.clear();
    //
    m_solver_qss.analyzePattern(m_rop_qss);
    // done
    m_QSS_init = true;
}

void GasQSSKinetics::calc_conc_QSS(doublereal* conc_qss)
{
    init_QSS();
    if (m_QSS_ok) return;
    // zero-lize m_ropf and m_ropr for m_rodf_mqss and m_rodr_mqss
    for (auto const r : m_rodf_mqss) m_ropf[r] = 0.;
    for (auto const r : m_rodr_mqss) m_ropr[r] = 0.;
    // compute rate of destruction
    m_rod_qss.setZero();
    for (size_t i = 0; i < m_nSpeciesQSS; i++) {
        for (const auto r : m_rodf_qss[i])
            m_rod_qss[i] += m_ropf[r];
        for (const auto r : m_rodr_qss[i])
            m_rod_qss[i] += m_ropr[r];
    }
    // compute rate of production from non-qss
    m_rop_noqss.setZero();
    for (size_t i = 0; i < m_nSpeciesQSS; i++) {
        for (const auto r : m_ropf_noqss[i])
            m_rop_noqss[i] += m_ropf[r];
        for (const auto r : m_ropr_noqss[i])
            m_rop_noqss[i] += m_ropr[r];
    }
    // compute rate of production from qss
    vector<char>::iterator it_ift = m_ifr_qss.begin();
    vector<vector<size_t>>::iterator it_ropf_qss = m_ropf_qss.begin();
    vector<vector<size_t>>::iterator it_ropr_qss = m_ropr_qss.begin();
    for (int k = 0; k < m_rop_qss.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(m_rop_qss,k); it; ++it) {
            it.valueRef() = 0.;
            if (it.col() == it.row()) {
                it.valueRef() +=  m_rod_qss[it.row()];
                continue;
            }
            if (*it_ift & IROPF) { // forward
                for (const auto r : *it_ropf_qss) it.valueRef() -= m_ropf[r];
                ++it_ropf_qss;
            }
            if (*it_ift & IROPR) { // reverse
                for (const auto r : *it_ropr_qss) it.valueRef() -= m_ropr[r];
                ++it_ropr_qss;
            }
            ++it_ift;
        }
    // solve linear system
    m_solver_qss.factorize(m_rop_qss);
    Map<VectorXd>(conc_qss, m_nSpeciesQSS) = m_solver_qss.solve(m_rop_noqss);

    m_QSS_ok = true;
}

void GasQSSKinetics::update_ROP_QSS(const doublereal* conc_qss)
{
    for (size_t i = 0; i < m_nSpeciesQSS; i++) {
        for (const auto r : m_rodf_qss[i])
            m_ropf[r] *= conc_qss[i];
        for (const auto r : m_rodr_qss[i])
            m_ropr[r] *= conc_qss[i];
    }
}

void GasQSSKinetics::printQSS()
{
    // QSS species
    cout << "QSSA Species: ";
    for (const auto& i : thermo(1).speciesNames())
        cout << i << " ";
    cout << endl;

    // m_rodf_qss
    cout << "m_rodf_qss: " << endl;
    for (size_t i = 0; i < m_rodf_qss.size(); ++i) {
        cout << i << ": ";
        for (const auto r : m_rodf_qss[i])
            cout << r << " ";
        cout << endl;
    }

    // m_rodr_qss
    cout << "m_rodr_qss: " << endl;
    for (size_t i = 0; i < m_rodr_qss.size(); ++i) {
        cout << i << ": ";
        for (const auto r : m_rodr_qss[i])
            cout << r << " ";
        cout << endl;
    }

    // m_ropf_noqss
    cout << "m_ropf_noqss: " << endl;
    for (size_t i = 0; i < m_ropf_noqss.size(); ++i) {
        cout << i << ": ";
        for (const auto r : m_ropf_noqss[i])
            cout << r << " ";
        cout << endl;
    }

    // m_ropr_noqss
    cout << "m_ropr_noqss: " << endl;
    for (size_t i = 0; i < m_ropr_noqss.size(); ++i) {
        cout << i << ": ";
        for (const auto r : m_ropr_noqss[i])
            cout << r << " ";
        cout << endl;
    }

    // m_rop_qss
    cout << "m_rop_qss: " << endl;
    vector<char>::iterator it_ift = m_ifr_qss.begin();
    vector<vector<size_t>>::iterator it_ropf_qss = m_ropf_qss.begin();
    vector<vector<size_t>>::iterator it_ropr_qss = m_ropr_qss.begin();
    for (int k = 0; k < m_rop_qss.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(m_rop_qss,k); it; ++it) {
            if (it.col() == it.row()) continue;
            if (*it_ift & IROPF) { // forward
                cout << "forward " << it.col() << ", " << it.row() << ": ";
                for (const auto r : *it_ropf_qss) {
                    cout << r << " ";
                }
                ++it_ropf_qss;
                cout << endl;
            }
            if (*it_ift & IROPR) { // reverse
                cout << "reverse " << it.col() << ", " << it.row() << ": ";
                for (const auto r : *it_ropr_qss) {
                    cout << r << " ";
                }
                ++it_ropr_qss;
                cout << endl;
            }
            ++it_ift;
        }
}

}
