// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

BulkKinetics::BulkKinetics() {
    setDerivativeSettings(AnyMap()); // use default settings
}

bool BulkKinetics::isReversible(size_t i) {
    return std::find(m_revindex.begin(), m_revindex.end(), i) < m_revindex.end();
}

bool BulkKinetics::addReaction(shared_ptr<Reaction> r, bool resize)
{
    bool added = Kinetics::addReaction(r, resize);
    if (!added) {
        // undeclared species, etc.
        return false;
    }
    double dn = 0.0;
    for (const auto& [name, stoich] : r->products) {
        dn += stoich;
    }
    for (const auto& [name, stoich] : r->reactants) {
        dn -= stoich;
    }

    m_dn.push_back(dn);

    if (r->reversible) {
        m_revindex.push_back(nReactions()-1);
    } else {
        m_irrev.push_back(nReactions()-1);
    }

    shared_ptr<ReactionRate> rate = r->rate();
    string rtype = rate->subType();
    if (rtype == "") {
        rtype = rate->type();
    }

    // If necessary, add new MultiRate evaluator
    if (m_bulk_types.find(rtype) == m_bulk_types.end()) {
        m_bulk_types[rtype] = m_bulk_rates.size();
        m_bulk_rates.push_back(rate->newMultiRate());
        m_bulk_rates.back()->resize(m_kk, nReactions(), nPhases());
    }

    // Set index of rate to number of reaction within kinetics
    rate->setRateIndex(nReactions() - 1);
    rate->setContext(*r, *this);

    // Add reaction rate to evaluator
    size_t index = m_bulk_types[rtype];
    m_bulk_rates[index]->add(nReactions() - 1, *rate);

    // Add reaction to third-body evaluator
    if (r->thirdBody() != nullptr) {
        addThirdBody(r);
    }

    m_concm.push_back(NAN);
    m_ready = resize;
    return true;
}

void BulkKinetics::addThirdBody(shared_ptr<Reaction> r)
{
    map<size_t, double> efficiencies;
    for (const auto& [name, efficiency] : r->thirdBody()->efficiencies) {
        size_t k = kineticsSpeciesIndex(name);
        if (k != npos) {
            efficiencies[k] = efficiency;
        } else if (!m_skipUndeclaredThirdBodies) {
            throw CanteraError("BulkKinetics::addThirdBody", "Found third-body"
                " efficiency for undefined species '{}' while adding reaction '{}'",
                name, r->equation());
        } else {
            m_hasUndeclaredThirdBodies = true;
        }
    }
    m_multi_concm.install(nReactions() - 1, efficiencies,
                          r->thirdBody()->default_efficiency,
                          r->thirdBody()->mass_action);
}

void BulkKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    // operations common to all reaction types
    Kinetics::modifyReaction(i, rNew);

    shared_ptr<ReactionRate> rate = rNew->rate();
    string rtype = rate->subType();
    if (rtype == "") {
        rtype = rate->type();
    }

    // Ensure that MultiRate evaluator is available
    if (m_bulk_types.find(rtype) == m_bulk_types.end()) {
        throw CanteraError("BulkKinetics::modifyReaction",
                "Evaluator not available for type '{}'.", rtype);
    }

    // Replace reaction rate to evaluator
    size_t index = m_bulk_types[rtype];
    rate->setRateIndex(i);
    rate->setContext(*rNew, *this);
    m_bulk_rates[index]->replace(i, *rate);
    invalidateCache();
}

void BulkKinetics::resizeSpecies()
{
    Kinetics::resizeSpecies();
    m_act_conc.resize(m_kk);
    m_phys_conc.resize(m_kk);
    m_grt.resize(m_kk);
    for (auto& rates : m_bulk_rates) {
        rates->resize(m_kk, nReactions(), nPhases());
    }
}

void BulkKinetics::resizeReactions()
{
    Kinetics::resizeReactions();
    m_rbuf0.resize(nReactions());
    m_rbuf1.resize(nReactions());
    m_rbuf2.resize(nReactions());
    m_kf0.resize(nReactions());
    m_sbuf0.resize(nTotalSpecies());
    m_state.resize(thermo().stateSize());
    m_multi_concm.resizeCoeffs(nTotalSpecies(), nReactions());
    for (auto& rates : m_bulk_rates) {
        rates->resize(nTotalSpecies(), nReactions(), nPhases());
        // @todo ensure that ReactionData are updated; calling rates->update
        //      blocks correct behavior in update_rates_T
        //      and running updateROP() is premature
    }
}

void BulkKinetics::setMultiplier(size_t i, double f)
{
    Kinetics::setMultiplier(i, f);
    m_ROP_ok = false;
}

void BulkKinetics::invalidateCache()
{
    Kinetics::invalidateCache();
    m_ROP_ok = false;
}

void BulkKinetics::getFwdRateConstants(double* kfwd)
{
    updateROP();
    copy(m_rfn.begin(), m_rfn.end(), kfwd);
    if (legacy_rate_constants_used()) {
        processThirdBodies(kfwd);
    }
}

void BulkKinetics::getEquilibriumConstants(double* kc)
{
    updateROP();

    vector<double>& delta_gibbs0 = m_rbuf0;
    fill(delta_gibbs0.begin(), delta_gibbs0.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), delta_gibbs0.data());

    double rrt = 1.0 / thermo().RT();
    double logStandConc = log(thermo().standardConcentration());
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-delta_gibbs0[i] * rrt + m_dn[i] * logStandConc);
    }
}

void BulkKinetics::getRevRateConstants(double* krev, bool doIrreversible)
{
    // go get the forward rate constants. -> note, we don't really care about
    // speed or redundancy in these informational routines.
    getFwdRateConstants(krev);

    if (doIrreversible) {
        getEquilibriumConstants(m_rbuf0.data());
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] /= m_rbuf0[i];
        }
    } else {
        // m_rkcn[] is zero for irreversible reactions
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

void BulkKinetics::getDeltaGibbs(double* deltaG)
{
    // Get the chemical potentials for each species
    thermo().getChemPotentials(m_sbuf0.data());
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaG);
}

void BulkKinetics::getDeltaEnthalpy(double* deltaH)
{
    // Get the partial molar enthalpy for each species
    thermo().getPartialMolarEnthalpies(m_sbuf0.data());
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaH);
}

void BulkKinetics::getDeltaEntropy(double* deltaS)
{
    // Get the partial molar entropy for each species
    thermo().getPartialMolarEntropies(m_sbuf0.data());
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaS);
}

void BulkKinetics::getDeltaSSGibbs(double* deltaG)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity. We define these here as
    // the chemical potentials of the pure species at the temperature and
    // pressure of the solution.
    thermo().getStandardChemPotentials(m_sbuf0.data());
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaG);
}

void BulkKinetics::getDeltaSSEnthalpy(double* deltaH)
{
    // Get the standard state enthalpies of the species.
    thermo().getEnthalpy_RT(m_sbuf0.data());
    for (size_t k = 0; k < m_kk; k++) {
        m_sbuf0[k] *= thermo().RT();
    }
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaH);
}

void BulkKinetics::getDeltaSSEntropy(double* deltaS)
{
    // Get the standard state entropy of the species. We define these here as
    // the entropies of the pure species at the temperature and pressure of the
    // solution.
    thermo().getEntropy_R(m_sbuf0.data());
    for (size_t k = 0; k < m_kk; k++) {
        m_sbuf0[k] *= GasConstant;
    }
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_sbuf0.data(), deltaS);
}

void BulkKinetics::getDerivativeSettings(AnyMap& settings) const
{
    settings["skip-third-bodies"] = m_jac_skip_third_bodies;
    settings["skip-falloff"] = m_jac_skip_falloff;
    settings["rtol-delta"] = m_jac_rtol_delta;
}

void BulkKinetics::setDerivativeSettings(const AnyMap& settings)
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

void BulkKinetics::getFwdRateConstants_ddT(double* dkfwd)
{
    assertDerivativesValid("BulkKinetics::getFwdRateConstants_ddT");
    updateROP();
    process_ddT(m_rfn, dkfwd);
}

void BulkKinetics::getFwdRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("BulkKinetics::getFwdRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropf, drop);
}

void BulkKinetics::getRevRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("BulkKinetics::getRevRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropr, drop);
    Eigen::Map<Eigen::VectorXd> dRevRop(drop, nReactions());

    // reverse rop times scaled inverse equilibrium constant derivatives
    Eigen::Map<Eigen::VectorXd> dRevRop2(m_rbuf2.data(), nReactions());
    copy(m_ropr.begin(), m_ropr.end(), m_rbuf2.begin());
    applyEquilibriumConstants_ddT(dRevRop2.data());
    dRevRop += dRevRop2;
}

void BulkKinetics::getNetRatesOfProgress_ddT(double* drop)
{
    assertDerivativesValid("BulkKinetics::getNetRatesOfProgress_ddT");
    updateROP();
    process_ddT(m_ropnet, drop);
    Eigen::Map<Eigen::VectorXd> dNetRop(drop, nReactions());

    // reverse rop times scaled inverse equilibrium constant derivatives
    Eigen::Map<Eigen::VectorXd> dNetRop2(m_rbuf2.data(), nReactions());
    copy(m_ropr.begin(), m_ropr.end(), m_rbuf2.begin());
    applyEquilibriumConstants_ddT(dNetRop2.data());
    dNetRop -= dNetRop2;
}

void BulkKinetics::getFwdRateConstants_ddP(double* dkfwd)
{
    assertDerivativesValid("BulkKinetics::getFwdRateConstants_ddP");
    updateROP();
    process_ddP(m_rfn, dkfwd);
}

void BulkKinetics::getFwdRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("BulkKinetics::getFwdRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropf, drop);
}

void BulkKinetics::getRevRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("BulkKinetics::getRevRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropr, drop);
}

void BulkKinetics::getNetRatesOfProgress_ddP(double* drop)
{
    assertDerivativesValid("BulkKinetics::getNetRatesOfProgress_ddP");
    updateROP();
    process_ddP(m_ropnet, drop);
}

void BulkKinetics::getFwdRateConstants_ddC(double* dkfwd)
{
    assertDerivativesValid("BulkKinetics::getFwdRateConstants_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_rfn, dkfwd, false);
}

void BulkKinetics::getFwdRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("BulkKinetics::getFwdRatesOfProgress_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_ropf, drop);
}

void BulkKinetics::getRevRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("BulkKinetics::getRevRatesOfProgress_ddC");
    updateROP();
    return process_ddC(m_revProductStoich, m_ropr, drop);
}

void BulkKinetics::getNetRatesOfProgress_ddC(double* drop)
{
    assertDerivativesValid("BulkKinetics::getNetRatesOfProgress_ddC");
    updateROP();
    process_ddC(m_reactantStoich, m_ropf, drop);
    Eigen::Map<Eigen::VectorXd> dNetRop(drop, nReactions());

    process_ddC(m_revProductStoich, m_ropr, m_rbuf2.data());
    Eigen::Map<Eigen::VectorXd> dNetRop2(m_rbuf2.data(), nReactions());
    dNetRop -= dNetRop2;
}

Eigen::SparseMatrix<double> BulkKinetics::fwdRatesOfProgress_ddX()
{
    assertDerivativesValid("BulkKinetics::fwdRatesOfProgress_ddX");

    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_reactantStoich, rop_rates);
}

Eigen::SparseMatrix<double> BulkKinetics::revRatesOfProgress_ddX()
{
    assertDerivativesValid("BulkKinetics::revRatesOfProgress_ddX");

    // reverse reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    applyEquilibriumConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_revProductStoich, rop_rates);
}

Eigen::SparseMatrix<double> BulkKinetics::netRatesOfProgress_ddX()
{
    assertDerivativesValid("BulkKinetics::netRatesOfProgress_ddX");

    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    auto jac = calculateCompositionDerivatives(m_reactantStoich, rop_rates);

    // reverse reaction rate coefficients
    applyEquilibriumConstants(rop_rates.data());
    return jac - calculateCompositionDerivatives(m_revProductStoich, rop_rates);
}

Eigen::SparseMatrix<double> BulkKinetics::fwdRatesOfProgress_ddCi()
{
    assertDerivativesValid("BulkKinetics::fwdRatesOfProgress_ddCi");

    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_reactantStoich, rop_rates, false);
}

Eigen::SparseMatrix<double> BulkKinetics::revRatesOfProgress_ddCi()
{
    assertDerivativesValid("BulkKinetics::revRatesOfProgress_ddCi");

    // reverse reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    applyEquilibriumConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_revProductStoich, rop_rates, false);
}

Eigen::SparseMatrix<double> BulkKinetics::netRatesOfProgress_ddCi()
{
    assertDerivativesValid("BulkKinetics::netRatesOfProgress_ddCi");

    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    auto jac = calculateCompositionDerivatives(m_reactantStoich, rop_rates, false);

    // reverse reaction rate coefficients
    applyEquilibriumConstants(rop_rates.data());
    return jac - calculateCompositionDerivatives(m_revProductStoich, rop_rates, false);
}

void BulkKinetics::updateROP()
{
    static const int cacheId = m_cache.getId();
    CachedScalar last = m_cache.getScalar(cacheId);
    double T = thermo().temperature();
    double rho = thermo().density();
    int statenum = thermo().stateMFNumber();

    if (last.state1 != T || last.state2 != rho) {
        // Update properties that are independent of the composition
        thermo().getStandardChemPotentials(m_grt.data());
        fill(m_delta_gibbs0.begin(), m_delta_gibbs0.end(), 0.0);
        double logStandConc = log(thermo().standardConcentration());

        // compute Delta G^0 for all reversible reactions
        getRevReactionDelta(m_grt.data(), m_delta_gibbs0.data());

        double rrt = 1.0 / thermo().RT();
        for (size_t i = 0; i < m_revindex.size(); i++) {
            size_t irxn = m_revindex[i];
            m_rkcn[irxn] = std::min(
                exp(m_delta_gibbs0[irxn] * rrt - m_dn[irxn] * logStandConc), BigNumber);
        }

        for (size_t i = 0; i != m_irrev.size(); ++i) {
            m_rkcn[ m_irrev[i] ] = 0.0;
        }
    }

    if (!last.validate(T, rho, statenum)) {
        // Update terms dependent on species concentrations and temperature
        thermo().getActivityConcentrations(m_act_conc.data());
        thermo().getConcentrations(m_phys_conc.data());
        double ctot = thermo().molarDensity();

        // Third-body objects interacting with MultiRate evaluator
        m_multi_concm.update(m_phys_conc, ctot, m_concm.data());
        m_ROP_ok = false;
    }

    // loop over MultiRate evaluators for each reaction type
    for (auto& rates : m_bulk_rates) {
        bool changed = rates->update(thermo(), *this);
        if (changed) {
            rates->getRateConstants(m_kf0.data());
            m_ROP_ok = false;
        }
    }

    if (m_ROP_ok) {
        // rates of progress are up-to-date only if both the thermodynamic state
        // and m_perturb are unchanged
        return;
    }

    // Scale the forward rate coefficient by the perturbation factor
    for (size_t i = 0; i < nReactions(); ++i) {
        m_rfn[i] = m_kf0[i] * m_perturb[i];
    }

    copy(m_rfn.begin(), m_rfn.end(), m_ropf.data());
    processThirdBodies(m_ropf.data());
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_act_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    applyEquilibriumConstants(m_ropr.data());
    m_revProductStoich.multiply(m_act_conc.data(), m_ropr.data());
    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    for (size_t i = 0; i < m_rfn.size(); i++) {
        AssertFinite(m_rfn[i], "BulkKinetics::updateROP",
                     "m_rfn[{}] is not finite.", i);
        AssertFinite(m_ropf[i], "BulkKinetics::updateROP",
                     "m_ropf[{}] is not finite.", i);
        AssertFinite(m_ropr[i], "BulkKinetics::updateROP",
                     "m_ropr[{}] is not finite.", i);
    }
    m_ROP_ok = true;
}

void BulkKinetics::getThirdBodyConcentrations(double* concm)
{
    updateROP();
    std::copy(m_concm.begin(), m_concm.end(), concm);
}

void BulkKinetics::processThirdBodies(double* rop)
{
    // reactions involving third body
    if (!m_concm.empty()) {
        m_multi_concm.multiply(rop, m_concm.data());
    }
}

void BulkKinetics::applyEquilibriumConstants(double* rop)
{
    // For reverse rates computed from thermochemistry, multiply the forward
    // rate coefficients by the reciprocals of the equilibrium constants
    for (size_t i = 0; i < nReactions(); ++i) {
        rop[i] *= m_rkcn[i];
    }
}

void BulkKinetics::applyEquilibriumConstants_ddT(double* drkcn)
{
    double T = thermo().temperature();
    double P = thermo().pressure();
    double rrt = 1. / thermo().RT();

    vector<double>& grt = m_sbuf0;
    vector<double>& delta_gibbs0 = m_rbuf1;
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

void BulkKinetics::process_ddT(const vector<double>& in, double* drop)
{
    // apply temperature derivative
    copy(in.begin(), in.end(), drop);
    for (auto& rates : m_bulk_rates) {
        rates->processRateConstants_ddT(drop, m_rfn.data(), m_jac_rtol_delta);
    }
}

void BulkKinetics::process_ddP(const vector<double>& in, double* drop)
{
    // apply pressure derivative
    copy(in.begin(), in.end(), drop);
    for (auto& rates : m_bulk_rates) {
        rates->processRateConstants_ddP(drop, m_rfn.data(), m_jac_rtol_delta);
    }
}

void BulkKinetics::process_ddC(StoichManagerN& stoich, const vector<double>& in,
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

Eigen::SparseMatrix<double> BulkKinetics::calculateCompositionDerivatives(
    StoichManagerN& stoich, const vector<double>& in, bool ddX)
{
    Eigen::SparseMatrix<double> out;
    vector<double>& scaled = m_rbuf1;
    vector<double>& outV = m_rbuf2;

    // convert from concentration to mole fraction output
    copy(in.begin(), in.end(), scaled.begin());
    if (ddX) {
        double ctot = thermo().molarDensity();
        for (size_t i = 0; i < nReactions(); ++i) {
            scaled[i] *= ctot;
        }
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

void BulkKinetics::assertDerivativesValid(const string& name)
{
    if (!thermo().isIdeal()) {
        throw NotImplementedError(name,
            "Not supported for non-ideal ThermoPhase models.");
    }
}

}
