/**
 *  @file GasKinetics.cpp
 *
 * Homogeneous kinetics in ideal gases
 *
 */

// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/GasKinetics.h"

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Enhanced3BConc.h"
#include "cantera/kinetics/ThirdBodyMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"

//#include "../user/grirxnstoich.h"

#include <iostream>
using namespace std;


namespace Cantera
{

//====================================================================================================================
GasKineticsData::GasKineticsData() :
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_ROP_ok(false),
    m_temp(0.0)
{
}
//====================================================================================================================
GasKineticsData::GasKineticsData(const GasKineticsData& right) :
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_ROP_ok(false),
    m_temp(0.0)
{
    *this = right;
}
//====================================================================================================================
GasKineticsData::~GasKineticsData()
{
}
//====================================================================================================================
GasKineticsData&  GasKineticsData::operator=(const GasKineticsData& right)
{
    if (this == &right) {
        return *this;
    }

    m_logp_ref = right.m_logp_ref;
    m_logc_ref  = right.m_logc_ref;
    m_logStandConc = right.m_logStandConc;
    m_ropf = right.m_ropf;
    m_ropr = right.m_ropr;
    m_ropnet = right.m_ropnet;
    m_rfn_low = right.m_rfn_low;
    m_rfn_high = right.m_rfn_high;
    m_ROP_ok  = right.m_ROP_ok;
    m_temp = right.m_temp;
    m_rfn  = right.m_rfn;
    falloff_work = right.falloff_work;
    concm_3b_values = right.concm_3b_values;
    concm_falloff_values = right.concm_falloff_values;
    m_rkcn = right.m_rkcn;

    return *this;
}
//====================================================================================================================
/*
 * Construct an empty reaction mechanism.
 */
GasKinetics::
GasKinetics(thermo_t* thermo) :
    Kinetics(),
    m_nfall(0),
    m_nirrev(0),
    m_nrev(0),
    m_finalized(false)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
    m_kdata = new GasKineticsData();
    m_kdata->m_temp = 0.0;
    m_rxnstoich = new ReactionStoichMgr();
}

//====================================================================================================================
GasKinetics::GasKinetics(const GasKinetics& right) :
    Kinetics(),
    m_nfall(0),
    m_nirrev(0),
    m_nrev(0),
    m_finalized(false)
{
    m_kdata = new GasKineticsData();
    m_kdata->m_temp = 0.0;
    m_rxnstoich = new ReactionStoichMgr();
    *this = right;
}
//====================================================================================================================
GasKinetics::~GasKinetics()
{
    delete m_kdata;
    delete m_rxnstoich;
}
//====================================================================================================================
GasKinetics& GasKinetics::operator=(const GasKinetics& right)
{
    if (this == &right) {
        return *this;
    }

    Kinetics::operator=(right);

    m_nfall = right.m_nfall;
    m_fallindx = right.m_fallindx;
    m_falloff_low_rates = right.m_falloff_low_rates;
    m_falloff_high_rates = right.m_falloff_high_rates;
    m_rates = right.m_rates;
    m_index = right.m_index;
    m_falloffn = right.m_falloffn;
    m_3b_concm = right.m_3b_concm;
    m_falloff_concm = right.m_falloff_concm;
    m_irrev = right.m_irrev;

    *m_rxnstoich = *(right.m_rxnstoich);

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

    *m_kdata = *(right.m_kdata);

    m_conc = right.m_conc;
    m_grt = right.m_grt;
    m_finalized = right.m_finalized;

    throw CanteraError("GasKinetics::operator=()",
                       "Unfinished implementation");

    return *this;
}
//====================================================================================================================
// Duplication routine for objects which inherit from Kinetics
/*
 *  This virtual routine can be used to duplicate %Kinetics objects
 *  inherited from %Kinetics even if the application only has
 *  a pointer to %Kinetics to work with.
 *
 *  These routines are basically wrappers around the derived copy
 *  constructor.
 *
 * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
 *                  m_thermo vector within this object
 */
Kinetics* GasKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    GasKinetics* gK = new GasKinetics(*this);
    gK->assignShallowPointers(tpVector);
    return dynamic_cast<Kinetics*>(gK);
}
//====================================================================================================================
/**
 * Update temperature-dependent portions of reaction rates and
 * falloff functions.
 */
void GasKinetics::update_T()
{
}
//====================================================================================================================
void GasKinetics::
update_C() {}
//====================================================================================================================
void GasKinetics::
_update_rates_T()
{
    doublereal T = thermo().temperature();
    m_kdata->m_logStandConc = log(thermo().standardConcentration());
    doublereal logT = log(T);
    if (!m_kdata->m_rfn.empty()) {
        m_rates.update(T, logT, &m_kdata->m_rfn[0]);
    }
    if (!m_kdata->m_rfn_low.empty()) {
        m_falloff_low_rates.update(T, logT, &m_kdata->m_rfn_low[0]);
        m_falloff_high_rates.update(T, logT, &m_kdata->m_rfn_high[0]);
    }
    if (!m_kdata->falloff_work.empty()) {
        m_falloffn.updateTemp(T, &m_kdata->falloff_work[0]);
    }
    m_kdata->m_temp = T;
    updateKc();
    m_kdata->m_ROP_ok = false;
    //}
};

//====================================================================================================================
/**
 * Update properties that depend on concentrations. Currently only
 * the enhanced collision partner concentrations are updated here.
 */
void GasKinetics::
_update_rates_C()
{
    thermo().getActivityConcentrations(&m_conc[0]);
    doublereal ctot = thermo().molarDensity();
    if (!m_kdata->concm_3b_values.empty()) {
        m_3b_concm.update(m_conc, ctot, &m_kdata->concm_3b_values[0]);
    }
    if (!m_kdata->concm_falloff_values.empty()) {
        m_falloff_concm.update(m_conc, ctot,
                               &m_kdata->concm_falloff_values[0]);
    }
    m_kdata->m_ROP_ok = false;
}
//====================================================================================================================
/**
 * Update the equilibrium constants in molar units.
 */
void GasKinetics::updateKc()
{
    vector_fp& m_rkc = m_kdata->m_rkcn;

    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(m_rkc.begin(), m_rkc.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    m_rxnstoich->getRevReactionDelta(m_ii, &m_grt[0], &m_rkc[0]);

    doublereal logStandConc = m_kdata->m_logStandConc;
    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_nrev; i++) {
        size_t irxn = m_revindex[i];
        m_rkc[irxn] = exp(m_rkc[irxn]*rrt - m_dn[irxn]*logStandConc);
    }

    for (size_t i = 0; i != m_nirrev; ++i) {
        m_rkc[ m_irrev[i] ] = 0.0;
    }
}
//====================================================================================================================
/**
 * Get the equilibrium constants of all reactions, whether
 * reversible or not.
 */
void GasKinetics::getEquilibriumConstants(doublereal* kc)
{
    _update_rates_T();
    vector_fp& rkc = m_kdata->m_rkcn;
    //thermo().getGibbs_RT(m_grt.begin());
    thermo().getStandardChemPotentials(&m_grt[0]);
    fill(rkc.begin(), rkc.end(), 0.0);

    // compute Delta G^0 for all reactions
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], &rkc[0]);

    doublereal logStandConc = m_kdata->m_logStandConc;
    doublereal rrt = 1.0/(GasConstant * thermo().temperature());
    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-rkc[i]*rrt + m_dn[i]*logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_kdata->m_temp = 0.0;
}
//====================================================================================================================
/**
 *
 * getDeltaGibbs():
 *
 * Return the vector of values for the reaction gibbs free energy
 * change
 * These values depend upon the concentration
 * of the ideal gas.
 *
 *  units = J kmol-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaG);
}
//====================================================================================================================
/**
 *
 * getDeltaEnthalpy():
 *
 * Return the vector of values for the reactions change in
 * enthalpy.
 * These values depend upon the concentration
 * of the solution.
 *
 *  units = J kmol-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaH);
}
//====================================================================================================================
/*
 *
 * getDeltaEntropy():
 *
 * Return the vector of values for the reactions change in
 * entropy.
 * These values depend upon the concentration
 * of the solution.
 *
 *  units = J kmol-1 Kelvin-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaS);
}
//====================================================================================================================
/**
 *
 * getDeltaSSGibbs():
 *
 * Return the vector of values for the reaction
 * standard state gibbs free energy change.
 * These values don't depend upon the concentration
 * of the solution.
 *
 *  units = J kmol-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaG);
}
//====================================================================================================================
/**
 *
 * getDeltaSSEnthalpy():
 *
 * Return the vector of values for the change in the
 * standard state enthalpies of reaction.
 * These values don't depend upon the concentration
 * of the solution.
 *
 *  units = J kmol-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaH);
}
//====================================================================================================================
/*********************************************************************
 *
 * getDeltaSSEntropy():
 *
 * Return the vector of values for the change in the
 * standard state entropies for each reaction.
 * These values don't depend upon the concentration
 * of the solution.
 *
 *  units = J kmol-1 Kelvin-1
 */
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
    m_rxnstoich->getReactionDelta(m_ii, &m_grt[0], deltaS);
}

//====================================================================================================================
// Return the species net production rates
/*
 * Species net production rates [kmol/m^3/s]. Return the species
 * net production rates (creation - destruction) in array
 * wdot, which must be dimensioned at least as large as the
 * total number of species.
 *
 *  @param net  Array of species production rates.
 *             units kmol m-3 s-1
 */
void GasKinetics::getNetProductionRates(doublereal* net)
{
    updateROP();
    m_rxnstoich->getNetProductionRates(m_kk, &m_kdata->m_ropnet[0], net);
}
//====================================================================================================================
// Return the species creation rates
/*
 * Species creation rates [kmol/m^3]. Return the species
 * creation rates in array cdot, which must be
 * dimensioned at least as large as the total number of
 * species.
 *
 *  @param cdot  Array of species production rates.
 *              units kmol m-3 s-1
 */
void GasKinetics::getCreationRates(doublereal* cdot)
{
    updateROP();
    m_rxnstoich->getCreationRates(m_kk, &m_kdata->m_ropf[0], &m_kdata->m_ropr[0], cdot);
}
//====================================================================================================================
//   Return a vector of the species destruction rates
/*
 * Species destruction rates [kmol/m^3]. Return the species
 * destruction rates in array ddot, which must be
 * dimensioned at least as large as the total number of
 * species.
 *
 *
 *  @param ddot  Array of species destruction rates.
 *               units kmol m-3 s-1
 *
 */
void GasKinetics::getDestructionRates(doublereal* ddot)
{
    updateROP();
    m_rxnstoich->getDestructionRates(m_kk, &m_kdata->m_ropf[0], &m_kdata->m_ropr[0], ddot);
}
//====================================================================================================================
void GasKinetics::processFalloffReactions()
{
    const vector_fp& fc = m_kdata->concm_falloff_values;
    const vector_fp& m_rf_low = m_kdata->m_rfn_low;
    const vector_fp& m_rf_high = m_kdata->m_rfn_high;

    // use m_ropr for temporary storage of reduced pressure
    vector_fp& pr = m_kdata->m_ropr;

    vector_fp& ropf = m_kdata->m_ropf;

    for (size_t i = 0; i < m_nfall; i++) {
        pr[i] = fc[i] * m_rf_low[i] / m_rf_high[i];
    }

    double* falloff_work =
        (m_kdata->falloff_work.empty()) ? 0 : &m_kdata->falloff_work[0];
    m_falloffn.pr_to_falloff(&pr[0], falloff_work);

    for (size_t i = 0; i < m_nfall; i++) {
        pr[i] *= m_rf_high[i];
    }

    scatter_copy(pr.begin(), pr.begin() + m_nfall,
                 ropf.begin(), m_fallindx.begin());
}

//====================================================================================================================
void GasKinetics::updateROP()
{

    _update_rates_T();
    _update_rates_C();

    if (m_kdata->m_ROP_ok) {
        return;
    }

    const vector_fp& rf = m_kdata->m_rfn;
    const vector_fp& m_rkc = m_kdata->m_rkcn;
    vector_fp& ropf = m_kdata->m_ropf;
    vector_fp& ropr = m_kdata->m_ropr;
    vector_fp& ropnet = m_kdata->m_ropnet;

    // copy rate coefficients into ropf
    copy(rf.begin(), rf.end(), ropf.begin());

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!m_kdata->concm_3b_values.empty()) {
        m_3b_concm.multiply(&ropf[0], &m_kdata->concm_3b_values[0]);
    }

    if (m_nfall) {
        processFalloffReactions();
    }

    // multiply by perturbation factor
    multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    copy(ropf.begin(), ropf.end(), ropr.begin());

    // for reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(ropr.begin(), ropr.end(), m_rkc.begin());

    // multiply ropf by concentration products
    m_rxnstoich->multiplyReactants(&m_conc[0], &ropf[0]);
    //m_reactantStoich.multiply(m_conc.begin(), ropf.begin());

    // for reversible reactions, multiply ropr by concentration
    // products
    m_rxnstoich->multiplyRevProducts(&m_conc[0], &ropr[0]);
    //m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

    for (size_t j = 0; j != m_ii; ++j) {
        ropnet[j] = ropf[j] - ropr[j];
    }

    m_kdata->m_ROP_ok = true;
}
//====================================================================================================================
/**
 *
 * getFwdRateConstants():
 *
 * Update the rate of progress for the reactions.
 * This key routine makes sure that the rate of progress vectors
 * located in the solid kinetics data class are up to date.
 */
void GasKinetics::
getFwdRateConstants(doublereal* kfwd)
{
    _update_rates_T();
    _update_rates_C();

    // copy rate coefficients into ropf
    const vector_fp& rf = m_kdata->m_rfn;
    vector_fp& ropf = m_kdata->m_ropf;
    copy(rf.begin(), rf.end(), ropf.begin());

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!m_kdata->concm_3b_values.empty()) {
        m_3b_concm.multiply(&ropf[0], &m_kdata->concm_3b_values[0]);
    }

    /*
     * This routine is hardcoded to replace some of the values
     * of the ropf vector.
     */
    if (m_nfall) {
        processFalloffReactions();
    }

    // multiply by perturbation factor
    multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());

    for (size_t i = 0; i < m_ii; i++) {
        kfwd[i] = ropf[i];
    }
}
//====================================================================================================================
/**
 *
 * getRevRateConstants():
 *
 * Return a vector of the reverse reaction rate constants
 *
 * Length is the number of reactions. units depends
 * on many issues. Note, this routine will return rate constants
 * for irreversible reactions if the default for
 * doIrreversible is overridden.
 */
void GasKinetics::
getRevRateConstants(doublereal* krev, bool doIrreversible)
{
    /*
     * go get the forward rate constants. -> note, we don't
     * really care about speed or redundancy in these
     * informational routines.
     */
    getFwdRateConstants(krev);

    if (doIrreversible) {
        doublereal* tmpKc = &m_kdata->m_ropnet[0];
        getEquilibriumConstants(tmpKc);
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] /=  tmpKc[i];
        }
    } else {
        /*
         * m_rkc[] is zero for irreversibly reactions
         */
        const vector_fp& m_rkc = m_kdata->m_rkcn;
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] *= m_rkc[i];
        }
    }
}
//====================================================================================================================
void GasKinetics::
addReaction(ReactionData& r)
{
    switch (r.reactionType) {
    case ELEMENTARY_RXN:
        addElementaryReaction(r);
        break;
    case THREE_BODY_RXN:
        addThreeBodyReaction(r);
        break;
    case FALLOFF_RXN:
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
}

//====================================================================================================================
void GasKinetics::
addFalloffReaction(ReactionData& r)
{
    // install high and low rate coeff calculators
    // and add constant terms to high and low rate coeff value vectors
    size_t iloc = m_falloff_high_rates.install(m_nfall, r);
    m_kdata->m_rfn_high.push_back(r.rateCoeffParameters[0]);
    std::swap(r.rateCoeffParameters, r.auxRateCoeffParameters);
    m_falloff_low_rates.install(m_nfall, r);
    m_kdata->m_rfn_low.push_back(r.rateCoeffParameters[0]);

    // add a dummy entry in m_rf, where computed falloff
    // rate coeff will be put
    m_kdata->m_rfn.push_back(0.0);

    // add this reaction number to the list of
    // falloff reactions
    m_fallindx.push_back(reactionNumber());

    // install the enhanced third-body concentration
    // calculator for this reaction
    m_falloff_concm.install(m_nfall, r.thirdBodyEfficiencies,
                            r.default_3b_eff);

    // install the falloff function calculator for
    // this reaction
    m_falloffn.install(m_nfall, r.falloffType, r.falloffParameters);

    // forward rxn order equals number of reactants, since rate
    // coeff is defined in terms of the high-pressure limit
    m_fwdOrder.push_back(r.reactants.size());

    // increment the falloff reaction counter
    ++m_nfall;
    registerReaction(reactionNumber(), FALLOFF_RXN, iloc);
}
//====================================================================================================================

void GasKinetics::
addElementaryReaction(ReactionData& r)
{
    size_t iloc;

    // install rate coeff calculator
    iloc = m_rates.install(reactionNumber(), r);

    // add constant term to rate coeff value vector
    m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);

    // forward rxn order equals number of reactants
    m_fwdOrder.push_back(r.reactants.size());
    registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);
}

//====================================================================================================================
void GasKinetics::
addThreeBodyReaction(ReactionData& r)
{
    size_t iloc;
    // install rate coeff calculator
    iloc = m_rates.install(reactionNumber(), r);

    // add constant term to rate coeff value vector
    m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);

    // forward rxn order equals number of reactants + 1
    m_fwdOrder.push_back(r.reactants.size() + 1);

    m_3b_concm.install(reactionNumber(), r.thirdBodyEfficiencies,
                       r.default_3b_eff);
    registerReaction(reactionNumber(), THREE_BODY_RXN, iloc);
}
//====================================================================================================================

void GasKinetics::addPlogReaction(ReactionData& r)
{
    // @todo: Not yet implemented
}

void GasKinetics::addChebyshevReaction(ReactionData& r)
{
    // @todo: Not yet implemented
}

void GasKinetics::installReagents(const ReactionData& r)
{

    m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn
    m_kdata->m_ropr.push_back(0.0);
    m_kdata->m_ropnet.push_back(0.0);
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

    m_kdata->m_rkcn.push_back(0.0);

    m_rxnstoich->add(reactionNumber(), r);

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
//====================================================================================================================


void GasKinetics::installGroups(size_t irxn,
                                const vector<grouplist_t>& r, const vector<grouplist_t>& p)
{
    if (!r.empty()) {
        writelog("installing groups for reaction "+int2str(reactionNumber()));
        m_rgroups[reactionNumber()] = r;
        m_pgroups[reactionNumber()] = p;
    }
}

//====================================================================================================================
void GasKinetics::init()
{
    m_kk = thermo().nSpecies();
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_conc.resize(m_kk);
    m_grt.resize(m_kk);
    m_kdata->m_logp_ref = log(thermo().refPressure()) - log(GasConstant);
}
//====================================================================================================================
void GasKinetics::finalize()
{
    if (!m_finalized) {
        //            int i, j, nr, np;
        m_kdata->falloff_work.resize(
            static_cast<size_t>(m_falloffn.workSize()));
        m_kdata->concm_3b_values.resize(
            static_cast<size_t>(m_3b_concm.workSize()));
        m_kdata->concm_falloff_values.resize(
            static_cast<size_t>(m_falloff_concm.workSize()));

        //             for (i = 0; i < m_ii; i++) {
        //                 nr = m_reactants[i].size();
        //                 for (j = 0; j < nr; j++) {
        //                     m_rstoich[i][m_reactants[i][j]]++;
        //                 }
        //                 np = m_products[i].size();
        //                 for (j = 0; j < np; j++) {
        //                     m_pstoich[i][m_products[i][j]]++;
        //                 }
        //             }
        //m_rxnstoich->write("c.cpp");
        m_finalized = true;
    }
}
//====================================================================================================================
bool GasKinetics::ready() const
{
    return (m_finalized);
}
//====================================================================================================================
}
//======================================================================================================================
