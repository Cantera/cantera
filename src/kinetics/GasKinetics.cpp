/**
 *  @file GasKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/GasKinetics.h"
#include "cantera/thermo/ThermoPhase.h"

using namespace std;

namespace Cantera
{
GasKinetics::GasKinetics(ThermoPhase* thermo) :
    BulkKinetics(thermo),
    m_logStandConc(0.0),
    m_pres(0.0)
{
    setDerivativeSettings(AnyMap()); // use default settings
}

void GasKinetics::resizeReactions()
{
    m_rbuf0.resize(nReactions());
    m_rbuf1.resize(nReactions());
    m_rbuf2.resize(nReactions());
    m_sbuf0.resize(nTotalSpecies());
    m_state.resize(thermo().stateSize());

    BulkKinetics::resizeReactions();
}

void GasKinetics::getThirdBodyConcentrations(double* concm)
{
    updateROP();
    std::copy(m_concm.begin(), m_concm.end(), concm);
}

void GasKinetics::update_rates_T()
{
    double T = thermo().temperature();
    double P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());

    if (T != m_temp) {
        updateKc();
        m_ROP_ok = false;
    }

    // loop over MultiRate evaluators for each reaction type
    for (auto& rates : m_bulk_rates) {
        bool changed = rates->update(thermo(), *this);
        if (changed) {
            rates->getRateConstants(m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void GasKinetics::update_rates_C()
{
    thermo().getActivityConcentrations(m_act_conc.data());
    thermo().getConcentrations(m_phys_conc.data());
    doublereal ctot = thermo().molarDensity();

    // Third-body objects interacting with MultiRate evaluator
    m_multi_concm.update(m_phys_conc, ctot, m_concm.data());
    m_ROP_ok = false;
}

void GasKinetics::updateKc()
{
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_delta_gibbs0.begin(), m_delta_gibbs0.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_delta_gibbs0.data());

    double rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = std::min(
            exp(m_delta_gibbs0[irxn] * rrt - m_dn[irxn] * m_logStandConc),
            BigNumber);
    }

    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void GasKinetics::processFwdRateCoefficients(double* ropf)
{
    update_rates_C();
    update_rates_T();

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), ropf);

    // Scale the forward rate coefficient by the perturbation factor
    for (size_t i = 0; i < nReactions(); ++i) {
        ropf[i] *= m_perturb[i];
    }
}

void GasKinetics::processThirdBodies(double* rop)
{
    // reactions involving third body
    if (!m_concm.empty()) {
        m_multi_concm.multiply(rop, m_concm.data());
    }
}

void GasKinetics::processEquilibriumConstants(double* rop)
{
    // For reverse rates computed from thermochemistry, multiply the forward
    // rate coefficients by the reciprocals of the equilibrium constants
    for (size_t i = 0; i < nReactions(); ++i) {
        rop[i] *= m_rkcn[i];
    }
}

void GasKinetics::getEquilibriumConstants(doublereal* kc)
{
    update_rates_T(); // this step ensures that m_grt is updated

    vector_fp& delta_gibbs0 = m_rbuf0;
    fill(delta_gibbs0.begin(), delta_gibbs0.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), delta_gibbs0.data());

    double rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-delta_gibbs0[i] * rrt + m_dn[i] * m_logStandConc);
    }
}

void GasKinetics::updateROP()
{
    processFwdRateCoefficients(m_ropf.data());
    processThirdBodies(m_ropf.data());
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_act_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    processEquilibriumConstants(m_ropr.data());
    m_revProductStoich.multiply(m_act_conc.data(), m_ropr.data());
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

void GasKinetics::getFwdRateConstants(double* kfwd)
{
    processFwdRateCoefficients(m_ropf.data());

    if (legacy_rate_constants_used()) {
        processThirdBodies(m_ropf.data());
    }

    // copy result
    copy(m_ropf.begin(), m_ropf.end(), kfwd);
}

void GasKinetics::getDerivativeSettings(AnyMap& settings) const
{
    settings["skip-third-bodies"] = m_jac_skip_third_bodies;
    settings["skip-falloff"] = m_jac_skip_falloff;
    settings["rtol-delta"] = m_jac_rtol_delta;
}

void GasKinetics::setDerivativeSettings(const AnyMap& settings)
{
    bool force = settings.empty();
    if (force || settings.hasKey("skip-third-bodies")) {
        m_jac_skip_third_bodies = settings.getBool("skip-third-bodies", false);
    }
    if (force || settings.hasKey("skip-falloff")) {
        m_jac_skip_falloff = settings.getBool("skip-falloff", false);
    }
    if (force || settings.hasKey("rtol-delta")) {
        m_jac_rtol_delta = settings.getDouble("rtol-delta", 1e-8);
    }
}

void GasKinetics::assertDerivativesValid(const std::string& name)
{
    if (!thermo().isIdeal()) {
        throw NotImplementedError(name,
            "Not supported for non-ideal ThermoPhase models.");
    }
}

void GasKinetics::processEquilibriumConstants_ddT(double* drkcn)
{
    double T = thermo().temperature();
    double P = thermo().pressure();
    double rrt = 1. / thermo().RT();

    vector_fp& grt = m_sbuf0;
    vector_fp& delta_gibbs0 = m_rbuf1;
    fill(delta_gibbs0.begin(), delta_gibbs0.end(), 0.0);

    // compute perturbed Delta G^0 for all reversible reactions
    thermo().saveState(m_state);
    thermo().setState_TP(T * (1. + m_jac_rtol_delta), P);
    thermo().getStandardChemPotentials(grt.data());
    getRevReactionDelta(grt.data(), delta_gibbs0.data());

    // apply scaling for derivative of inverse equilibrium constant
    double Tinv = 1. / T;
    double rrt_dTinv = rrt * Tinv / m_jac_rtol_delta;
    double rrtt = rrt * Tinv;
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        double factor = delta_gibbs0[irxn] - m_delta_gibbs0[irxn];
        factor *= rrt_dTinv;
        factor += m_dn[irxn] * Tinv - m_delta_gibbs0[irxn] * rrtt;
        drkcn[irxn] *= factor;
    }

    for (size_t i = 0; i < m_irrev.size(); ++i) {
        drkcn[m_irrev[i]] = 0.0;
    }

    thermo().restoreState(m_state);
}

void GasKinetics::process_ddT(const vector_fp& in, double* drop)
{
    // apply temperature derivative
    copy(in.begin(), in.end(), drop);
    for (auto& rates : m_bulk_rates) {
        rates->processRateConstants_ddT(drop, m_rfn.data(), m_jac_rtol_delta);
    }
}

void GasKinetics::getFwdRateConstants_ddT(double* dkfwd)
{
    assertDerivativesValid("GasKinetics::getFwdRateConstants_ddT");
    updateROP();
    process_ddT(m_rfn, dkfwd);
}

void GasKinetics::getFwdRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("GasKinetics::getFwdRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropf, drop);
}

void GasKinetics::getRevRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("GasKinetics::getRevRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropr, drop);
    Eigen::Map<Eigen::VectorXd> dRevRop(drop, nReactions());

    // reverse rop times scaled inverse equilibrium constant derivatives
    Eigen::Map<Eigen::VectorXd> dRevRop2(m_rbuf2.data(), nReactions());
    copy(m_ropr.begin(), m_ropr.end(), m_rbuf2.begin());
    processEquilibriumConstants_ddT(dRevRop2.data());
    dRevRop += dRevRop2;
}

void GasKinetics::getNetRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("GasKinetics::getNetRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropnet, drop);
    Eigen::Map<Eigen::VectorXd> dNetRop(drop, nReactions());

    // reverse rop times scaled inverse equilibrium constant derivatives
    Eigen::Map<Eigen::VectorXd> dNetRop2(m_rbuf2.data(), nReactions());
    copy(m_ropr.begin(), m_ropr.end(), m_rbuf2.begin());
    processEquilibriumConstants_ddT(dNetRop2.data());
    dNetRop -= dNetRop2;
}

void GasKinetics::process_ddP(const vector_fp& in, double* drop)
{
    // apply pressure derivative
    copy(in.begin(), in.end(), drop);
    for (auto& rates : m_bulk_rates) {
        rates->processRateConstants_ddP(drop, m_rfn.data(), m_jac_rtol_delta);
    }
}

void GasKinetics::getFwdRateConstants_ddP(double* dkfwd)
{
    assertDerivativesValid("GasKinetics::getFwdRateConstants_ddP");
    updateROP();
    process_ddP(m_rfn, dkfwd);
}

void GasKinetics::getFwdRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("GasKinetics::getFwdRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropf, drop);
}

void GasKinetics::getRevRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("GasKinetics::getRevRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropr, drop);
}

void GasKinetics::getNetRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("GasKinetics::getNetRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropnet, drop);
}

void GasKinetics::process_ddC(
    StoichManagerN& stoich, const vector_fp& in,
    double* drop, bool mass_action)
{
    Eigen::Map<Eigen::VectorXd> out(drop, nReactions());
    out.setZero();
    double ctot_inv = 1. / thermo().molarDensity();

    // derivatives due to concentrations in law of mass action
    if (mass_action) {
        stoich.scale(in.data(), out.data(), ctot_inv);
    }
    if (m_jac_skip_third_bodies || m_multi_concm.empty()) {
        return;
    }

    // derivatives due to third-body colliders in law of mass action
    Eigen::Map<Eigen::VectorXd> outM(m_rbuf1.data(), nReactions());
    if (mass_action) {
        outM.fill(0.);
        m_multi_concm.scale(in.data(), outM.data(), ctot_inv);
        out += outM;
    }

    // derivatives due to reaction rates depending on third-body colliders
    if (!m_jac_skip_falloff) {
        m_multi_concm.scaleM(in.data(), outM.data(), m_concm.data(), ctot_inv);
        for (auto& rates : m_bulk_rates) {
            // processing step assigns zeros to entries not dependent on M
            rates->processRateConstants_ddM(
                outM.data(), m_rfn.data(), m_jac_rtol_delta);
        }
        out += outM;
    }
}

void GasKinetics::getFwdRateConstants_ddC(double* dkfwd)
{
    assertDerivativesValid("GasKinetics::getFwdRateConstants_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_rfn, dkfwd, false);
}

void GasKinetics::getFwdRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("GasKinetics::getFwdRatesOfProgress_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_ropf, drop);
}

void GasKinetics::getRevRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("GasKinetics::getRevRatesOfProgress_ddC");
    updateROP();
    return process_ddC(m_revProductStoich, m_ropr, drop);
}

void GasKinetics::getNetRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("GasKinetics::getNetRatesOfProgress_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_ropf, drop);
    Eigen::Map<Eigen::VectorXd> dNetRop(drop, nReactions());

    process_ddC(m_revProductStoich, m_ropr, m_rbuf2.data());
    Eigen::Map<Eigen::VectorXd> dNetRop2(m_rbuf2.data(), nReactions());
    dNetRop -= dNetRop2;
}

Eigen::SparseMatrix<double> GasKinetics::process_ddX(
    StoichManagerN& stoich, const vector_fp& in)
{
    Eigen::SparseMatrix<double> out;
    vector_fp& scaled = m_rbuf1;
    vector_fp& outV = m_rbuf2;

    // convert from concentration to mole fraction output
    double ctot = thermo().molarDensity();
    for (size_t i = 0; i < nReactions(); ++i) {
        scaled[i] = ctot * in[i];
    }

    // derivatives handled by StoichManagerN
    copy(scaled.begin(), scaled.end(), outV.begin());
    processThirdBodies(outV.data());
    out = stoich.derivatives(m_act_conc.data(), outV.data());
    if (m_jac_skip_third_bodies || m_multi_concm.empty()) {
        return out;
    }

    // derivatives due to law of mass action
    copy(scaled.begin(), scaled.end(), outV.begin());
    stoich.multiply(m_act_conc.data(), outV.data());

    // derivatives due to reaction rates depending on third-body colliders
    if (!m_jac_skip_falloff) {
        for (auto& rates : m_bulk_rates) {
            // processing step does not modify entries not dependent on M
            rates->processRateConstants_ddM(
                outV.data(), m_rfn.data(), m_jac_rtol_delta, false);
        }
    }

    // derivatives handled by ThirdBodyCalc
    out += m_multi_concm.derivatives(outV.data());

    return out;
}

Eigen::SparseMatrix<double> GasKinetics::fwdRatesOfProgress_ddX()
{
    assertDerivativesValid("GasKinetics::fwdRatesOfProgress_ddX");

    // forward reaction rate coefficients
    vector_fp& rop_rates = m_rbuf0;
    processFwdRateCoefficients(rop_rates.data());
    return process_ddX(m_reactantStoich, rop_rates);
}

Eigen::SparseMatrix<double> GasKinetics::revRatesOfProgress_ddX()
{
    assertDerivativesValid("GasKinetics::revRatesOfProgress_ddX");

    // reverse reaction rate coefficients
    vector_fp& rop_rates = m_rbuf0;
    processFwdRateCoefficients(rop_rates.data());
    processEquilibriumConstants(rop_rates.data());
    return process_ddX(m_revProductStoich, rop_rates);
}

Eigen::SparseMatrix<double> GasKinetics::netRatesOfProgress_ddX()
{
    assertDerivativesValid("GasKinetics::netRatesOfProgress_ddX");

    // forward reaction rate coefficients
    vector_fp& rop_rates = m_rbuf0;
    processFwdRateCoefficients(rop_rates.data());
    Eigen::SparseMatrix<double> jac = process_ddX(m_reactantStoich, rop_rates);

    // reverse reaction rate coefficients
    processEquilibriumConstants(rop_rates.data());
    return jac - process_ddX(m_revProductStoich, rop_rates);
}

bool GasKinetics::addReaction(shared_ptr<Reaction> r, bool resize)
{
    // operations common to all reaction types
    bool added = BulkKinetics::addReaction(r, resize);
    if (!added) {
        return false;
    }
    return true;
}

void GasKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    // operations common to all bulk reaction types
    BulkKinetics::modifyReaction(i, rNew);

    // invalidate all cached data
    invalidateCache();
}

void GasKinetics::invalidateCache()
{
    BulkKinetics::invalidateCache();
    m_pres += 0.13579;
}

}
