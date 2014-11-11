/**
 *  @file GasKinetics.cpp
 *
 * Homogeneous kinetics in ideal gases
 */

// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/GasKinetics.h"

using namespace std;

namespace Cantera
{
GasKinetics::GasKinetics(thermo_t* thermo) :
    BulkKinetics(thermo),
    m_nfall(0),
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_pres(0.0)
{
}

Kinetics* GasKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    GasKinetics* gK = new GasKinetics(*this);
    gK->assignShallowPointers(tpVector);
    return gK;
}

void GasKinetics::update_rates_T()
{
    doublereal T = thermo().temperature();
    doublereal P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());
    doublereal logT = log(T);

    if (T != m_temp) {
        if (!m_rfn.empty()) {
            m_rates.update(T, logT, &m_rfn[0]);
        }

        if (!m_rfn_low.empty()) {
            m_falloff_low_rates.update(T, logT, &m_rfn_low[0]);
            m_falloff_high_rates.update(T, logT, &m_rfn_high[0]);
        }
        if (!falloff_work.empty()) {
            m_falloffn.updateTemp(T, &falloff_work[0]);
        }
        updateKc();
        m_ROP_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.update(T, logT, &m_rfn[0]);
            m_ROP_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.update(T, logT, &m_rfn[0]);
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void GasKinetics::update_rates_C()
{
    thermo().getActivityConcentrations(&m_conc[0]);
    doublereal ctot = thermo().molarDensity();

    // 3-body reactions
    if (!concm_3b_values.empty()) {
        m_3b_concm.update(m_conc, ctot, &concm_3b_values[0]);
    }

    // Falloff reactions
    if (!concm_falloff_values.empty()) {
        m_falloff_concm.update(m_conc, ctot, &concm_falloff_values[0]);
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
    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(&m_grt[0], &m_rkcn[0]);

    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
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
    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(&m_grt[0], &m_rkcn[0]);

    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-m_rkcn[i]*rrt + m_dn[i]*m_logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void GasKinetics::processFalloffReactions()
{
    // use m_ropr for temporary storage of reduced pressure
    vector_fp& pr = m_ropr;

    for (size_t i = 0; i < m_nfall; i++) {
        pr[i] = concm_falloff_values[i] * m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
        AssertFinite(pr[i], "GasKinetics::processFalloffReactions",
                     "pr[" + int2str(i) + "] is not finite.");
    }

    double* work = (falloff_work.empty()) ? 0 : &falloff_work[0];
    m_falloffn.pr_to_falloff(&pr[0], work);

    for (size_t i = 0; i < m_nfall; i++) {
        if (m_rxntype[m_fallindx[i]] == FALLOFF_RXN) {
            pr[i] *= m_rfn_high[i];
        } else { // CHEMACT_RXN
            pr[i] *= m_rfn_low[i];
        }
    }

    scatter_copy(pr.begin(), pr.begin() + m_nfall,
                 m_ropf.begin(), m_fallindx.begin());
}

void GasKinetics::updateROP()
{
    update_rates_C();
    update_rates_T();

    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(&m_ropf[0], &concm_3b_values[0]);
    }

    if (m_nfall) {
        processFalloffReactions();
    }

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // for reverse rates computed from thermochemistry, multiply the forward
    // rates copied into m_ropr by the reciprocals of the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by concentration products
    m_reactantStoich.multiply(&m_conc[0], &m_ropf[0]);

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich.multiply(&m_conc[0], &m_ropr[0]);

    for (size_t j = 0; j != m_ii; ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    for (size_t i = 0; i < m_rfn.size(); i++) {
        AssertFinite(m_rfn[i], "GasKinetics::updateROP",
                     "m_rfn[" + int2str(i) + "] is not finite.");
        AssertFinite(m_ropf[i], "GasKinetics::updateROP",
                     "m_ropf[" + int2str(i) + "] is not finite.");
        AssertFinite(m_ropr[i], "GasKinetics::updateROP",
                     "m_ropr[" + int2str(i) + "] is not finite.");
    }

    m_ROP_ok = true;
}

void GasKinetics::getFwdRateConstants(doublereal* kfwd)
{
    update_rates_C();
    update_rates_T();

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(&m_ropf[0], &concm_3b_values[0]);
    }

    if (m_nfall) {
        processFalloffReactions();
    }

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    for (size_t i = 0; i < m_ii; i++) {
        kfwd[i] = m_ropf[i];
    }
}

void GasKinetics::addReaction(ReactionData& r)
{
    switch (r.reactionType) {
    case ELEMENTARY_RXN:
        addElementaryReaction(r);
        break;
    case THREE_BODY_RXN:
        addThreeBodyReaction(r);
        break;
    case FALLOFF_RXN:
    case CHEMACT_RXN:
        addFalloffReaction(r);
        break;
    case PLOG_RXN:
        addPlogReaction(r);
        break;
    case CHEBYSHEV_RXN:
        addChebyshevReaction(r);
        break;
    default:
        throw CanteraError("GasKinetics::addReaction", "Invalid reaction type specified");
    }

    // operations common to all reaction types
    BulkKinetics::addReaction(r);
}

void GasKinetics::addReaction(shared_ptr<Reaction> r)
{
    switch (r->reaction_type) {
    case ELEMENTARY_RXN:
        addElementaryReaction(dynamic_cast<ElementaryReaction&>(*r));
        break;
    case THREE_BODY_RXN:
        addThreeBodyReaction(dynamic_cast<ThirdBodyReaction&>(*r));
        break;
    case FALLOFF_RXN:
    case CHEMACT_RXN:
        addFalloffReaction(dynamic_cast<FalloffReaction&>(*r));
        break;
    case PLOG_RXN:
        addPlogReaction(dynamic_cast<PlogReaction&>(*r));
        break;
    case CHEBYSHEV_RXN:
        addChebyshevReaction(dynamic_cast<ChebyshevReaction&>(*r));
        break;
    default:
        throw CanteraError("GasKinetics::addReaction",
            "Unknown reaction type specified: " + int2str(r->reaction_type));
    }

    // operations common to all reaction types
    BulkKinetics::addReaction(r);
}

void GasKinetics::addFalloffReaction(ReactionData& r)
{
    // install high and low rate coeff calculators
    // and add constant terms to high and low rate coeff value vectors
    m_falloff_high_rates.install(m_nfall, r);
    m_rfn_high.push_back(r.rateCoeffParameters[0]);
    std::swap(r.rateCoeffParameters, r.auxRateCoeffParameters);
    m_falloff_low_rates.install(m_nfall, r);
    m_rfn_low.push_back(r.rateCoeffParameters[0]);

    // add this reaction number to the list of falloff reactions
    m_fallindx.push_back(nReactions());

    // install the enhanced third-body concentration calculator for this
    // reaction
    m_falloff_concm.install(m_nfall, r.thirdBodyEfficiencies,
                            r.default_3b_eff);

    // install the falloff function calculator for this reaction
    m_falloffn.install(m_nfall, r.falloffType, r.reactionType,
                       r.falloffParameters);

    // increment the falloff reaction counter
    ++m_nfall;
}

void GasKinetics::addThreeBodyReaction(ReactionData& r)
{
    m_rates.install(nReactions(), r);
    m_3b_concm.install(nReactions(), r.thirdBodyEfficiencies,
                       r.default_3b_eff);
}

void GasKinetics::addPlogReaction(ReactionData& r)
{
    m_plog_rates.install(nReactions(), r);
}

void GasKinetics::addChebyshevReaction(ReactionData& r)
{
    m_cheb_rates.install(nReactions(), r);
}

void GasKinetics::addFalloffReaction(FalloffReaction& r)
{
    // install high and low rate coeff calculators
    // and extend the high and low rate coeff value vectors
    m_falloff_high_rates.install(m_nfall, r.high_rate);
    m_rfn_high.push_back(0.0);
    m_falloff_low_rates.install(m_nfall, r.low_rate);
    m_rfn_low.push_back(0.0);

    // add this reaction number to the list of falloff reactions
    m_fallindx.push_back(nReactions());

    // install the enhanced third-body concentration calculator
    map<size_t, double> efficiencies;
    for (Composition::const_iterator iter = r.third_body.efficiencies.begin();
         iter != r.third_body.efficiencies.end();
         ++iter) {
        efficiencies[kineticsSpeciesIndex(iter->first)] = iter->second;
    }
    m_falloff_concm.install(m_nfall, efficiencies,
                            r.third_body.default_efficiency);

    // install the falloff function calculator for this reaction
    m_falloffn.install(m_nfall, r.falloff_type, r.reaction_type,
                       r.falloff_parameters);

    // increment the falloff reaction counter
    ++m_nfall;
}

void GasKinetics::addThreeBodyReaction(ThirdBodyReaction& r)
{
    m_rates.install(nReactions(), r.rate);
    map<size_t, double> efficiencies;
    for (Composition::const_iterator iter = r.third_body.efficiencies.begin();
         iter != r.third_body.efficiencies.end();
         ++iter) {
        efficiencies[kineticsSpeciesIndex(iter->first)] = iter->second;
    }
    m_3b_concm.install(nReactions(), efficiencies,
                       r.third_body.default_efficiency);
}

void GasKinetics::addPlogReaction(PlogReaction& r)
{
    m_plog_rates.install(nReactions(), r.rate);
}

void GasKinetics::addChebyshevReaction(ChebyshevReaction& r)
{
    m_cheb_rates.install(nReactions(), r.rate);
}

void GasKinetics::init()
{
    BulkKinetics::init();
    m_logp_ref = log(thermo().refPressure()) - log(GasConstant);
}

void GasKinetics::finalize()
{
    BulkKinetics::finalize();
    falloff_work.resize(m_falloffn.workSize());
    concm_3b_values.resize(m_3b_concm.workSize());
    concm_falloff_values.resize(m_falloff_concm.workSize());
}

bool GasKinetics::ready() const
{
    return m_finalized;
}

}
