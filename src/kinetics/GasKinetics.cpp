/**
 *  @file GasKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <numeric>
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/StoichManager.h"
#include "cantera/kinetics/ThirdBodyCalc.h"
#include "cantera/thermo/ThermoPhase.h"

using namespace std;

namespace Cantera
{
GasKinetics::GasKinetics(ThermoPhase* thermo) :
    BulkKinetics(thermo),
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_pres(0.0)
{
    m_3b_concm = std::unique_ptr<ThirdBodyCalc>(new ThirdBodyCalc());
    m_falloff_concm = std::unique_ptr<ThirdBodyCalc>(new ThirdBodyCalc());
}

void GasKinetics::finalizeSetup()
{
    // Third-body calculators
    m_3b_concm->finalizeSetup(m_kk, nReactions());
    m_falloff_concm->finalizeSetup(m_kk, nReactions());

    m_rop_stoich.resize(nReactions());
    m_rop_3b.resize(nReactions());
    m_rop_rates.resize(nReactions());

    BulkKinetics::finalizeSetup();
}

void GasKinetics::update_rates_T()
{
    double T = thermo().temperature();
    double P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());
    double logT = log(T);

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
        if (m_blowersmasel_rates.nReactions()) {
            thermo().getPartialMolarEnthalpies(m_grt.data());
            getReactionDelta(m_grt.data(), m_dH.data());
            m_blowersmasel_rates.updateBlowersMasel(T, logT, m_rfn.data(), m_dH.data());
        }
    }

    if (T != m_temp || P != m_pres) {

        // loop over MultiBulkRates evaluators
        for (auto& rates : m_bulk_rates) {
            rates->update(thermo(), m_concm.data());
            rates->getRateConstants(thermo(), m_rfn.data(), m_concm.data());
        }

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

void GasKinetics::update_rates_C()
{
    thermo().getActivityConcentrations(m_act_conc.data());
    thermo().getConcentrations(m_phys_conc.data());
    doublereal ctot = thermo().molarDensity();

    // 3-body reactions
    if (!concm_3b_values.empty()) {
        m_3b_concm->update(m_phys_conc, ctot, concm_3b_values.data());
    }

    // Falloff reactions
    if (!concm_falloff_values.empty()) {
        m_falloff_concm->update(m_phys_conc, ctot, concm_falloff_values.data());
    }

    // Third-body objects interacting with MultiRate evaluator
    if (!concm_multi_values.empty()) {
        // using pre-existing third-body handlers requires copying;
        m_multi_concm->update(m_phys_conc, ctot, concm_multi_values.data());
        for (size_t i = 0; i < m_multi_indices.size(); i++) {
            m_concm[m_multi_indices[i]] = concm_multi_values[i];
        }
    }

    // P-log reactions
    if (m_plog_rates.nReactions()) {
        double logP = log(thermo().pressure());
        m_plog_rates.update_C(&logP);
    }

    // Chebyshev reactions
    if (m_cheb_rates.nReactions()) {
        double log10P = log10(thermo().pressure());
        m_cheb_rates.update_C(&log10P);
    }

    m_ROP_ok = false;
}

void GasKinetics::updateKc()
{
    thermo().getStandardChemPotentials(m_grt.data());
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

void GasKinetics::getEquilibriumConstants(doublereal* kc)
{
    update_rates_T();
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-m_rkcn[i]*rrt + m_dn[i]*m_logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void GasKinetics::processFalloffReactions(vector_fp& ropf)
{
    // use m_ropr for temporary storage of reduced pressure
    vector_fp& pr = m_ropr;

    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        pr[i] = concm_falloff_values[i] * m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
        AssertFinite(pr[i], "GasKinetics::processFalloffReactions",
                     "pr[{}] is not finite.", i);
    }

    m_falloffn.pr_to_falloff(pr.data(), falloff_work.data());

    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        if (reactionTypeStr(m_fallindx[i]) == "falloff") {
            pr[i] *= m_rfn_high[i];
        } else { // CHEMACT_RXN
            pr[i] *= m_rfn_low[i];
        }
        ropf[m_fallindx[i]] = pr[i];
    }
}

Eigen::SparseMatrix<double> GasKinetics::getFwdRopSpeciesDerivatives()
{
    Eigen::SparseMatrix<double> jac;

    if (!m_finalized) {
        finalizeSetup();
    }

    // forward reaction rate coefficients
    calcFwdRateCoefficients(m_rop_rates);

    // derivatives handled by StoichManagerN
    std::copy(m_rop_rates.begin(), m_rop_rates.end(), m_rop_stoich.begin());
    applyThirdBodies(m_rop_stoich);
    jac = m_reactantStoich->speciesDerivatives(
        m_act_conc.data(), m_rop_stoich.data());

    if (m_jac_skip_third_bodies) {
        return jac;
    }

    // derivatives handled by ThirdBodyCalc
    std::copy(m_rop_rates.begin(), m_rop_rates.end(), m_rop_3b.begin());
    m_reactantStoich->multiply(m_act_conc.data(), m_rop_3b.data());
    if (!concm_3b_values.empty()) {
        // Do not support legacy CTI/XML-based reaction rate evaluators
        throw CanteraError("GasKinetics::getFwdRopSpeciesDerivatives",
            "Not supported for legacy input format.");
    }
    if (!concm_multi_values.empty()) {
        jac += m_multi_concm->speciesDerivatives(m_rop_3b.data());
    }

    return jac;
}

Eigen::SparseMatrix<double> GasKinetics::getRevRopSpeciesDerivatives()
{
    Eigen::SparseMatrix<double> jac;

    if (!m_finalized) {
        finalizeSetup();
    }

    // forward reaction rate coefficients
    calcFwdRateCoefficients(m_rop_rates);
    calcRevRateCoefficients(m_rop_rates, m_rop_rates);

    // derivatives handled by StoichManagerN
    std::copy(m_rop_rates.begin(), m_rop_rates.end(), m_rop_stoich.begin());
    applyThirdBodies(m_rop_stoich);
    jac = m_revProductStoich->speciesDerivatives(
        m_act_conc.data(), m_rop_stoich.data());

    if (m_jac_skip_third_bodies) {
        return jac;
    }

    // derivatives handled by ThirdBodyCalc
    std::copy(m_rop_rates.begin(), m_rop_rates.end(), m_rop_3b.begin());
    m_productStoich->multiply(m_act_conc.data(), m_rop_3b.data());
    if (!concm_3b_values.empty()) {
        // Do not support legacy CTI/XML-based reaction rate evaluators
        throw CanteraError("GasKinetics::getRevRopSpeciesDerivatives",
            "Not supported for legacy input format.");
    }
    if (!concm_multi_values.empty()) {
        jac += m_multi_concm->speciesDerivatives(m_rop_3b.data());
    }

    return jac;
}

Eigen::SparseMatrix<double> GasKinetics::getNetRopSpeciesDerivatives()
{
    return getFwdRopSpeciesDerivatives() - getRevRopSpeciesDerivatives();
}

Eigen::VectorXd GasKinetics::getFwdRopTemperatureDerivatives()
{
    if (!m_jac_exact_temperature_derivatives) {
        return Kinetics::getFwdRopTemperatureDerivatives();
    }

    Eigen::VectorXd out(nReactions());
    updateROP();
    out.fill(NAN);
    for (auto& rates : m_bulk_rates) {
        rates->getScaledRateConstantTemperatureDerivatives(
            thermo(), out.data(), m_concm.data());
    }
    out.array() *= MappedVector(m_ropf.data(), nReactions()).array();
    return out;
}

Eigen::VectorXd GasKinetics::getRevRopTemperatureDerivatives()
{
    if (!m_jac_exact_temperature_derivatives) {
        return Kinetics::getRevRopTemperatureDerivatives();
    }

    // @TODO look at temperature derivatives of equilibrium constants
    throw CanteraError("GasKinetics::getRevRopTemperatureDerivatives",
        "Exact reverse rates-of-progress temperature derivatives "
        "are not implemented.\nUse numerical approximation instead.");
}

Eigen::VectorXd GasKinetics::getNetRopTemperatureDerivatives()
{
    if (!m_jac_exact_temperature_derivatives) {
        return Kinetics::getNetRopTemperatureDerivatives();
    }

    // @TODO look at temperature derivatives of equilibrium constants
    throw CanteraError("GasKinetics::getRevNetTemperatureDerivatives",
        "Exact net rates-of-progress temperature derivatives "
        "are not implemented.\nUse numerical approximation instead.");
}

void GasKinetics::updateROP()
{
    if (!m_finalized) {
        finalizeSetup();
    }

    calcFwdRateCoefficients(m_ropf);
    applyThirdBodies(m_ropf);
    calcRevRateCoefficients(m_ropf, m_ropr);

    // multiply ropf by concentration products
    m_reactantStoich->multiply(m_act_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich->multiply(m_act_conc.data(), m_ropr.data());

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

 void GasKinetics::calcFwdRateCoefficients(vector_fp& ropf)
 {
    update_rates_C();
    update_rates_T();

    // copy rate coefficients into ropf
    ropf = m_rfn;

    if (m_falloff_high_rates.nReactions()) {
        processFalloffReactions(ropf);
    }

    // Scale the forward rate coefficient by the perturbation factor
    for (size_t i = 0; i < nReactions(); i++) {
        ropf[i] *= m_perturb[i];
    }
}

void GasKinetics::applyThirdBodies(vector_fp& ropf)
{
    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm->multiply(ropf.data(), concm_3b_values.data());
    }

    // reactions involving third body
    if (!concm_multi_values.empty()) {
        m_multi_concm->multiply(ropf.data(), concm_multi_values.data());
    }
}

void GasKinetics::calcRevRateCoefficients(const vector_fp& ropf, vector_fp& ropr)
{
    // For reverse rates computed from thermochemistry, multiply the forward
    // rate coefficients by the reciprocals of the equilibrium constants
    for (size_t i = 0; i < nReactions(); i++) {
        ropr[i] = ropf[i] * m_rkcn[i];
    }
}

void GasKinetics::getFwdRateConstants(doublereal* kfwd)
{
    calcFwdRateCoefficients(m_ropf);
    if (legacy_rate_constants_used()) {
        warn_deprecated("GasKinetics::getFwdRateConstants",
            "Behavior to change after Cantera 2.6;\nresults will no longer include "
            "third-body concentrations for three-body reactions.\nTo switch to new "
            "behavior, use 'cantera.use_legacy_rate_constants(False)' (Python),\n"
            "'useLegacyRateConstants(0)' (MATLAB), 'Cantera::use_legacy_rate_constants"
            "(false)' (C++),\nor 'ct_use_legacy_rate_constants(0)' (clib).");

        applyThirdBodies(m_ropf);
    }

    // copy result
    for (size_t i = 0; i < nReactions(); i++) {
        kfwd[i] = m_ropf[i];
    }
}

bool GasKinetics::addReaction(shared_ptr<Reaction> r)
{
    // operations common to all reaction types
    bool added = BulkKinetics::addReaction(r);
    if (!added) {
        return false;
    } else if (!(r->usesLegacy())) {
        // Rate object already added in BulkKinetics::addReaction
        return true;
    }

    if (r->type() == "elementary-legacy") {
        addElementaryReaction(dynamic_cast<ElementaryReaction2&>(*r));
    } else if (r->type() == "three-body-legacy") {
        addThreeBodyReaction(dynamic_cast<ThreeBodyReaction2&>(*r));
    } else if (r->type() == "falloff") {
        addFalloffReaction(dynamic_cast<FalloffReaction&>(*r));
    } else if (r->type() == "chemically-activated") {
        addFalloffReaction(dynamic_cast<FalloffReaction&>(*r));
    } else if (r->type() == "pressure-dependent-Arrhenius-legacy") {
        addPlogReaction(dynamic_cast<PlogReaction2&>(*r));
    } else if (r->type() == "Chebyshev-legacy") {
        addChebyshevReaction(dynamic_cast<ChebyshevReaction2&>(*r));
    } else if (r->type() == "Blowers-Masel") {
        addBlowersMaselReaction(dynamic_cast<BlowersMaselReaction&>(*r));
    } else {
        throw CanteraError("GasKinetics::addReaction",
            "Unknown reaction type specified: '{}'", r->type());
    }
    return true;
}

void GasKinetics::addFalloffReaction(FalloffReaction& r)
{
    // install high and low rate coeff calculators and extend the high and low
    // rate coeff value vectors
    size_t nfall = m_falloff_high_rates.nReactions();
    m_falloff_high_rates.install(nfall, r.high_rate);
    m_rfn_high.push_back(0.0);
    m_falloff_low_rates.install(nfall, r.low_rate);
    m_rfn_low.push_back(0.0);

    // add this reaction number to the list of falloff reactions
    m_fallindx.push_back(nReactions()-1);
    m_rfallindx[nReactions()-1] = nfall;

    // install the enhanced third-body concentration calculator
    map<size_t, double> efficiencies;
    for (const auto& eff : r.third_body.efficiencies) {
        size_t k = kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            efficiencies[k] = eff.second;
        }
    }
    m_falloff_concm->install(nfall, efficiencies,
                            r.third_body.default_efficiency);
    concm_falloff_values.resize(m_falloff_concm->workSize());

    // install the falloff function calculator for this reaction
    m_falloffn.install(nfall, r.type(), r.falloff);
    falloff_work.resize(m_falloffn.workSize());
}

void GasKinetics::addThreeBodyReaction(ThreeBodyReaction2& r)
{
    m_rates.install(nReactions()-1, r.rate);
    map<size_t, double> efficiencies;
    for (const auto& eff : r.third_body.efficiencies) {
        size_t k = kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            efficiencies[k] = eff.second;
        }
    }
    m_3b_concm->install(nReactions()-1, efficiencies,
                       r.third_body.default_efficiency);
    concm_3b_values.resize(m_3b_concm->workSize());
}

void GasKinetics::addPlogReaction(PlogReaction2& r)
{
    m_plog_rates.install(nReactions()-1, r.rate);
}

void GasKinetics::addChebyshevReaction(ChebyshevReaction2& r)
{
    m_cheb_rates.install(nReactions()-1, r.rate);
}

void GasKinetics::addBlowersMaselReaction(BlowersMaselReaction& r)
{
    m_blowersmasel_rates.install(nReactions()-1, r.rate);
}

void GasKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    // operations common to all bulk reaction types
    BulkKinetics::modifyReaction(i, rNew);

    if (!(rNew->usesLegacy())) {
        // Rate object already modified in BulkKinetics::modifyReaction
        return;
    }

    if (rNew->type() == "elementary-legacy") {
        modifyElementaryReaction(i, dynamic_cast<ElementaryReaction2&>(*rNew));
    } else if (rNew->type() == "three-body-legacy") {
        modifyThreeBodyReaction(i, dynamic_cast<ThreeBodyReaction2&>(*rNew));
    } else if (rNew->type() == "falloff") {
        modifyFalloffReaction(i, dynamic_cast<FalloffReaction&>(*rNew));
    } else if (rNew->type() == "chemically-activated") {
        modifyFalloffReaction(i, dynamic_cast<FalloffReaction&>(*rNew));
    } else if (rNew->type() == "pressure-dependent-Arrhenius-legacy") {
        modifyPlogReaction(i, dynamic_cast<PlogReaction2&>(*rNew));
    } else if (rNew->type() == "Chebyshev-legacy") {
        modifyChebyshevReaction(i, dynamic_cast<ChebyshevReaction2&>(*rNew));
    } else if (rNew->type() == "Blowers-Masel") {
        modifyBlowersMaselReaction(i, dynamic_cast<BlowersMaselReaction&>(*rNew));
    } else {
        throw CanteraError("GasKinetics::modifyReaction",
            "Unknown reaction type specified: '{}'", rNew->type());
    }

    // invalidate all cached data
    m_ROP_ok = false;
    m_temp += 0.1234;
    m_pres += 0.1234;
}

void GasKinetics::modifyThreeBodyReaction(size_t i, ThreeBodyReaction2& r)
{
    m_rates.replace(i, r.rate);
}

void GasKinetics::modifyFalloffReaction(size_t i, FalloffReaction& r)
{
    size_t iFall = m_rfallindx[i];
    m_falloff_high_rates.replace(iFall, r.high_rate);
    m_falloff_low_rates.replace(iFall, r.low_rate);
    m_falloffn.replace(iFall, r.falloff);
}

void GasKinetics::modifyPlogReaction(size_t i, PlogReaction2& r)
{
    m_plog_rates.replace(i, r.rate);
}

void GasKinetics::modifyChebyshevReaction(size_t i, ChebyshevReaction2& r)
{
    m_cheb_rates.replace(i, r.rate);
}

void GasKinetics::modifyBlowersMaselReaction(size_t i, BlowersMaselReaction& r)
{
    m_blowersmasel_rates.replace(i, r.rate);
}

void GasKinetics::init()
{
    BulkKinetics::init();
    m_logp_ref = log(thermo().refPressure()) - log(GasConstant);
}

void GasKinetics::invalidateCache()
{
    BulkKinetics::invalidateCache();
    m_pres += 0.13579;
}

}
