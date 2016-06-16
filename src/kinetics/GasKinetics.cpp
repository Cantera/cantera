/**
 *  @file GasKinetics.cpp
 *
 * Homogeneous kinetics in ideal gases
 */

// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/GasKinetics.h"

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Enhanced3BConc.h"
#include "cantera/kinetics/ThirdBodyMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"

using namespace std;

namespace Cantera
{
GasKinetics::GasKinetics(thermo_t* thermo) :
    Kinetics(),
    m_nfall(0),
    m_ntedep(0),
    m_nvibrel(0),
    m_nirrev(0),
    m_nrev(0),
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_etemp(0.0),
    m_pres(0.0),
    m_finalized(false)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
    m_temp  = 0.0;
    m_etemp = 0.0;
}

GasKinetics::GasKinetics(const GasKinetics& right) :
    Kinetics(),
    m_nfall(0),
    m_ntedep(0),
    m_nvibrel(0),
    m_nirrev(0),
    m_nrev(0),
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_etemp(0.0),
    m_pres(0.0),
    m_finalized(false)
{
    m_temp  = 0.0;
    m_etemp = 0.0;
    *this   = right;
}

GasKinetics& GasKinetics::operator=(const GasKinetics& right)
{
    if (this == &right) {
        return *this;
    }

    Kinetics::operator=(right);

    m_nfall = right.m_nfall;
    m_ntedep = right.m_ntedep;
    m_nvibrel = right.m_nvibrel;
    m_fallindx = right.m_fallindx;
    m_tedepindx = m_tedepindx;
    m_vibrelindx = m_vibrelindx;
    m_falloff_low_rates = right.m_falloff_low_rates;
    m_falloff_high_rates = right.m_falloff_high_rates;
    m_rates = right.m_rates;
    m_tedep_rates = right.m_tedep_rates;
    m_vibrel_rates = right.m_vibrel_rates;
    m_index = right.m_index;
    m_falloffn = right.m_falloffn;
    m_3b_concm = right.m_3b_concm;
    m_falloff_concm = right.m_falloff_concm;
    m_irrev = right.m_irrev;
    m_plog_rates = right.m_plog_rates;
    m_cheb_rates = right.m_cheb_rates;

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
    m_reactantStrings = right.m_reactantStrings;
    m_productStrings = right.m_productStrings;

    m_logp_ref = right.m_logp_ref;
    m_logc_ref  = right.m_logc_ref;
    m_logStandConc = right.m_logStandConc;
    m_ropf = right.m_ropf;
    m_ropr = right.m_ropr;
    m_ropnet = right.m_ropnet;
    m_tedep_rfn = right.m_tedep_rfn;
    m_vibrel_rfn = right.m_vibrel_rfn;
    m_rfn_low = right.m_rfn_low;
    m_rfn_high = right.m_rfn_high;
    m_ROP_ok  = right.m_ROP_ok;
    m_temp = right.m_temp;
    m_etemp = right.m_etemp;
    m_rfn  = right.m_rfn;
    falloff_work = right.falloff_work;
    concm_3b_values = right.concm_3b_values;
    concm_falloff_values = right.concm_falloff_values;
    m_rkcn = right.m_rkcn;
    m_deltaE = right.m_deltaE;

    m_conc = right.m_conc;
    m_grt = right.m_grt;
    m_finalized = right.m_finalized;

    throw CanteraError("GasKinetics::operator=()",
                       "Unfinished implementation");

    return *this;
}

Kinetics* GasKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    GasKinetics* gK = new GasKinetics(*this);
    gK->assignShallowPointers(tpVector);
    return gK;
}

void GasKinetics::update_rates_Te()
{
    doublereal Te    = thermo().elec_temperature();
    m_logStandConc   = log(thermo().standardConcentration());
    doublereal logTe = log(Te);

    if (Te != m_etemp) {
        if (!m_tedep_rfn.empty()) {
            m_tedep_rates.update(Te, logTe, &m_tedep_rfn[0]);
        }

        updateKc();
        m_ROP_ok = false;
    }

    m_etemp = Te;
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
    m_rxnstoich.getRevReactionDelta(m_ii, &m_grt[0], &m_rkcn[0]);

    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_nrev; i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = std::min(exp(m_rkcn[irxn]*rrt - m_dn[irxn]*m_logStandConc),
                                BigNumber);
    }

    for (size_t i = 0; i != m_nirrev; ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void GasKinetics::getEquilibriumConstants(doublereal* kc)
{
    update_rates_T();
    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reactions
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], &m_rkcn[0]);

    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-m_rkcn[i]*rrt + m_dn[i]*m_logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void GasKinetics::getDeltaGibbs(doublereal* deltaG)
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

void GasKinetics::getDeltaEnthalpy(doublereal* deltaH)
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

void GasKinetics::getDeltaEntropy(doublereal* deltaS)
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

void GasKinetics::getDeltaSSGibbs(doublereal* deltaG)
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

void GasKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
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

void GasKinetics::getDeltaSSEntropy(doublereal* deltaS)
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

  void GasKinetics::getDeltaEPlasma(doublereal* deltaE)
{

  /*
   *  Get the energy release of each electron-impact reactions.
   */  
  for(size_t k = 0; k < m_ntedep; ++k) {
    deltaE[k] = m_deltaE[k];
  }
  
}

void GasKinetics::getNetProductionRates(doublereal* net)
{
    updateROP();
    m_rxnstoich.getNetProductionRates(m_kk, &m_ropnet[0], net);
}

void GasKinetics::getCreationRates(doublereal* cdot)
{
    updateROP();
    m_rxnstoich.getCreationRates(m_kk, &m_ropf[0], &m_ropr[0], cdot);
}

void GasKinetics::getDestructionRates(doublereal* ddot)
{
    updateROP();
    m_rxnstoich.getDestructionRates(m_kk, &m_ropf[0], &m_ropr[0], ddot);
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

    if (m_ntedep) {
      update_rates_Te();
      scatter_copy(m_tedep_rfn.begin(), m_tedep_rfn.end(),
		   m_ropf.begin(), m_tedepindx.begin());
    }

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

    /*
     * This routine is hardcoded to replace some of the values
     * of the ropf vector.
     */
    if (m_nfall) {
        processFalloffReactions();
    }

    if (m_ntedep) {
      update_rates_Te();
      scatter_copy(m_tedep_rfn.begin(), m_tedep_rfn.end(),
		   m_ropf.begin(), m_tedepindx.begin());
    }

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    for (size_t i = 0; i < m_ii; i++) {
        kfwd[i] = m_ropf[i];
    }
}

void GasKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
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
        // m_rkcn[] is zero for irreversible reactions
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

void GasKinetics::writeMech(const string& filename)
{

  string cppfile;
  size_t ns = thermo().nSpecies();

  cppfile = filename+".cpp";
  m_rxnstoich.writeMech(ns, cppfile, thermo().speciesThermo());
  
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
    case TEDEP_RXN:
        addTeDependentReaction(r);
	incrementTeDependentRxnCount();
        break;
    case VIBREL_RXN:
        addVibRelaxationReaction(r);
	incrementVibRelaxationRxnCount();
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
    installReagents(r);
    installGroups(reactionNumber(), r.rgroups, r.pgroups);
    incrementRxnCount();
    m_rxneqn.push_back(r.equation);
    m_reactantStrings.push_back(r.reactantString);
    m_productStrings.push_back(r.productString);
    m_rxntype.push_back(r.reactionType);
}

void GasKinetics::addFalloffReaction(ReactionData& r)
{
    // install high and low rate coeff calculators
    // and add constant terms to high and low rate coeff value vectors
    size_t iloc = m_falloff_high_rates.install(m_nfall, r);
    m_rfn_high.push_back(r.rateCoeffParameters[0]);
    std::swap(r.rateCoeffParameters, r.auxRateCoeffParameters);
    m_falloff_low_rates.install(m_nfall, r);
    m_rfn_low.push_back(r.rateCoeffParameters[0]);

    // add a dummy entry in m_rf, where computed falloff
    // rate coeff will be put
    m_rfn.push_back(0.0);

    // add this reaction number to the list of
    // falloff reactions
    m_fallindx.push_back(reactionNumber());

    // install the enhanced third-body concentration
    // calculator for this reaction
    m_falloff_concm.install(m_nfall, r.thirdBodyEfficiencies,
                            r.default_3b_eff);

    // install the falloff function calculator for
    // this reaction
    m_falloffn.install(m_nfall, r.falloffType, r.reactionType,
                       r.falloffParameters);

    // forward rxn order equals number of reactants, since rate
    // coeff is defined in terms of the high-pressure limit
    m_fwdOrder.push_back(r.reactants.size());

    // increment the falloff reaction counter
    ++m_nfall;
    registerReaction(reactionNumber(), r.reactionType, iloc);
}

void GasKinetics::addElementaryReaction(ReactionData& r)
{
    // install rate coeff calculator
    size_t iloc = m_rates.install(reactionNumber(), r);

    // add constant term to rate coeff value vector
    m_rfn.push_back(r.rateCoeffParameters[0]);

    // forward rxn order equals number of reactants
    m_fwdOrder.push_back(r.reactants.size());
    registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);
}

void GasKinetics::addThreeBodyReaction(ReactionData& r)
{
    // install rate coeff calculator
    size_t iloc = m_rates.install(reactionNumber(), r);

    // add constant term to rate coeff value vector
    m_rfn.push_back(r.rateCoeffParameters[0]);

    // forward rxn order equals number of reactants + 1
    m_fwdOrder.push_back(r.reactants.size() + 1);

    m_3b_concm.install(reactionNumber(), r.thirdBodyEfficiencies,
                       r.default_3b_eff);
    registerReaction(reactionNumber(), THREE_BODY_RXN, iloc);
}

void GasKinetics::addTeDependentReaction(ReactionData& r)
{
  // install rate coeff calculator
  size_t iloc = m_tedep_rates.install(reactionNumber(), r);

  // add this reaction number to the list of
  // te-dependent reactions
  m_tedepindx.push_back(reactionNumber());

  // add constant term to rate coeff value vector
  m_rfn.push_back(r.rateCoeffParameters[0]);
  m_tedep_rfn.push_back(r.rateCoeffParameters[0]);

  // forward rxn order equals number of reactants
  m_fwdOrder.push_back(r.reactants.size());

  // excitation energy
  m_deltaE.push_back(r.deltaE);

  registerReaction(reactionNumber(), TEDEP_RXN, iloc);

  ++m_ntedep;
  
}

void GasKinetics::addVibRelaxationReaction(ReactionData& r)
{
  // install rate coeff calculator
  size_t iloc = m_vibrel_rates.install(reactionNumber(), r);

  // add this reaction number to the list of
  // te-dependent reactions
  m_vibrelindx.push_back(reactionNumber());

  // add constant term to rate coeff value vector
  m_rfn.push_back(r.rateCoeffParameters[0]);
  m_vibrel_rfn.push_back(r.rateCoeffParameters[0]);

  // forward rxn order equals number of reactants
  m_fwdOrder.push_back(r.reactants.size());

  registerReaction(reactionNumber(), VIBREL_RXN, iloc);

  ++m_nvibrel;
  
}

void GasKinetics::addPlogReaction(ReactionData& r)
{
    // install rate coefficient calculator
    size_t iloc = m_plog_rates.install(reactionNumber(), r);

    // add a dummy entry in m_rfn, where computed rate coeff will be put
    m_rfn.push_back(0.0);

    m_fwdOrder.push_back(r.reactants.size());
    registerReaction(reactionNumber(), PLOG_RXN, iloc);
}

void GasKinetics::addChebyshevReaction(ReactionData& r)
{
    // install rate coefficient calculator
    size_t iloc = m_cheb_rates.install(reactionNumber(), r);

    // add a dummy entry in m_rfn, where computed rate coeff will be put
    m_rfn.push_back(0.0);

    m_fwdOrder.push_back(r.reactants.size());
    registerReaction(reactionNumber(), CHEBYSHEV_RXN, iloc);
}

void GasKinetics::installReagents(const ReactionData& r)
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
            ns = std::max<size_t>(ns, 1);
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
            ns = std::max<size_t>(ns, 1);
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

void GasKinetics::installGroups(size_t irxn,
                                const vector<grouplist_t>& r, const vector<grouplist_t>& p)
{
    if (!r.empty()) {
        writelog("installing groups for reaction "+int2str(reactionNumber()));
        m_rgroups[reactionNumber()] = r;
        m_pgroups[reactionNumber()] = p;
    }
}

void GasKinetics::init()
{
    m_kk = thermo().nSpecies();
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_conc.resize(m_kk);
    m_grt.resize(m_kk);
    m_logp_ref = log(thermo().refPressure()) - log(GasConstant);
}

void GasKinetics::finalize()
{
    if (!m_finalized) {
        falloff_work.resize(m_falloffn.workSize());
        concm_3b_values.resize(m_3b_concm.workSize());
        concm_falloff_values.resize(m_falloff_concm.workSize());
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

bool GasKinetics::ready() const
{
    return m_finalized;
}

}
