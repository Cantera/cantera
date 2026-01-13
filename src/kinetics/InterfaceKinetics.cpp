/**
 *  @file InterfaceKinetics.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorFactory.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

// Constructor / destructor definitions required due to forward-declared unique_ptr
// members.
InterfaceKinetics::~InterfaceKinetics()
{
    delete m_integrator;
}

void InterfaceKinetics::resizeReactions()
{
    Kinetics::resizeReactions();

    // resize buffer
    m_rbuf0.resize(nReactions());
    m_rbuf1.resize(nReactions());

    for (auto& rates : m_rateHandlers) {
        rates->resize(nTotalSpecies(), nReactions(), nPhases());
        // @todo ensure that ReactionData are updated; calling rates->update
        //      blocks correct behavior in InterfaceKinetics::_update_rates_T
        //      and running updateROP() is premature
    }
}

void InterfaceKinetics::setElectricPotential(int n, double V)
{
    thermo(n).setElectricPotential(V);
    m_redo_rates = true;
}

void InterfaceKinetics::_update_rates_T()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    // Go find the temperature from the surface
    double T = thermo(0).temperature();
    m_redo_rates = true;
    if (T != m_temp || m_redo_rates) {
        //  Calculate the forward rate constant by calling m_rates and store it in m_rfn[]
        for (size_t n = 0; n < nPhases(); n++) {
            thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
        }

        // Use the stoichiometric manager to find deltaH for each reaction.
        getReactionDelta(m_grt.data(), m_dH.data());

        m_temp = T;
        m_ROP_ok = false;
        m_redo_rates = false;
    }

    // loop over interface MultiRate evaluators for each reaction type
    for (auto& rates : m_rateHandlers) {
        bool changed = rates->update(thermo(0), *this);
        if (changed) {
            rates->getRateConstants(m_rfn.data());
            m_ROP_ok = false;
            m_redo_rates = true;
        }
    }

    if (!m_ROP_ok) {
        updateKc();
    }
}

void InterfaceKinetics::_update_rates_phi()
{
    // Store electric potentials for each phase in the array m_phi[].
    for (size_t n = 0; n < nPhases(); n++) {
        if (thermo(n).electricPotential() != m_phi[n]) {
            m_phi[n] = thermo(n).electricPotential();
            m_redo_rates = true;
        }
    }
}

void InterfaceKinetics::_update_rates_C()
{
    for (size_t n = 0; n < nPhases(); n++) {
        const auto& tp = thermo(n);
        /*
         * We call the getActivityConcentrations function of each ThermoPhase
         * class that makes up this kinetics object to obtain the generalized
         * concentrations for species within that class. This is collected in
         * the vector m_conc. m_start[] are integer indices for that vector
         * denoting the start of the species for each phase.
         */
        tp.getActivityConcentrations(m_actConc.data() + m_start[n]);

        // Get regular concentrations too
        tp.getConcentrations(m_conc.data() + m_start[n]);
    }
    m_ROP_ok = false;
}

void InterfaceKinetics::updateKc()
{
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    if (m_revindex.size() > 0) {
        /*
         * Get the vector of standard state electrochemical potentials for
         * species in the Interfacial kinetics object and store it in m_mu0[]
         * and m_mu0_Kc[]
         */
        updateMu0();
        double rrt = 1.0 / thermo(0).RT();

        // compute Delta mu^0 for all reversible reactions
        getRevReactionDelta(m_mu0_Kc.data(), m_rkcn.data());

        for (size_t i = 0; i < m_revindex.size(); i++) {
            size_t irxn = m_revindex[i];
            if (irxn == npos || irxn >= nReactions()) {
                throw CanteraError("InterfaceKinetics::updateKc",
                                   "illegal value: irxn = {}", irxn);
            }
            // WARNING this may overflow HKM
            m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
        }
        for (size_t i = 0; i != m_irrev.size(); ++i) {
            m_rkcn[ m_irrev[i] ] = 0.0;
        }
    }
}

void InterfaceKinetics::updateMu0()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    // @todo  There is significant potential to further simplify calculations
    //      once the old framework is removed
    size_t ik = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
        for (size_t k = 0; k < thermo(n).nSpecies(); k++) {
            m_mu0_Kc[ik] = m_mu0[ik] + Faraday * m_phi[n] * thermo(n).charge(k);
            m_mu0_Kc[ik] -= thermo(0).RT() * thermo(n).logStandardConc(k);
            ik++;
        }
    }
}

void InterfaceKinetics::getEquilibriumConstants(double* kc)
{
    updateMu0();
    double rrt = 1.0 / thermo(0).RT();
    std::fill(kc, kc + nReactions(), 0.0);
    getReactionDelta(m_mu0_Kc.data(), kc);
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-kc[i]*rrt);
    }
}

void InterfaceKinetics::getFwdRateConstants(double* kfwd)
{
    updateROP();
    for (size_t i = 0; i < nReactions(); i++) {
        // base rate coefficient multiplied by perturbation factor
        kfwd[i] = m_rfn[i] * m_perturb[i];
    }
}

void InterfaceKinetics::getRevRateConstants(double* krev, bool doIrreversible)
{
    getFwdRateConstants(krev);
    if (doIrreversible) {
        getEquilibriumConstants(m_ropnet.data());
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] /= m_ropnet[i];
        }
    } else {
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

void InterfaceKinetics::updateROP()
{
    // evaluate rate constants and equilibrium constants at temperature and phi
    // (electric potential)
    _update_rates_T();
    // get updated activities (rates updated below)
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    for (size_t i = 0; i < nReactions(); i++) {
        // Scale the base forward rate coefficient by the perturbation factor
        m_ropf[i] = m_rfn[i] * m_perturb[i];
        // Multiply the scaled forward rate coefficient by the reciprocal of the
        // equilibrium constant
        m_ropr[i] = m_ropf[i] * m_rkcn[i];
    }

    for (auto& rates : m_rateHandlers) {
        rates->modifyRateConstants(m_ropf.data(), m_ropr.data());
    }

    // multiply ropf by the activity concentration reaction orders to obtain
    // the forward rates of progress.
    m_reactantStoich.multiply(m_actConc.data(), m_ropf.data());

    // For reversible reactions, multiply ropr by the activity concentration
    // products
    m_revProductStoich.multiply(m_actConc.data(), m_ropr.data());

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    // For reactions involving multiple phases, we must check that the phase
    // being consumed actually exists. This is particularly important for phases
    // that are stoichiometric phases containing one species with a unity
    // activity
    if (m_phaseExistsCheck) {
        for (size_t j = 0; j != nReactions(); ++j) {
            if ((m_ropr[j] > m_ropf[j]) && (m_ropr[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsProduct[j][p] && !m_phaseExists[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropr[j] = m_ropf[j];
                        if (m_ropf[j] > 0.0) {
                            for (size_t rp = 0; rp < nPhases(); rp++) {
                                if (m_rxnPhaseIsReactant[j][rp] && !m_phaseExists[rp]) {
                                    m_ropnet[j] = 0.0;
                                    m_ropr[j] = m_ropf[j] = 0.0;
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsReactant[j][p] && !m_phaseIsStable[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropr[j] = m_ropf[j];
                    }
                }
            } else if ((m_ropf[j] > m_ropr[j]) && (m_ropf[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsReactant[j][p] && !m_phaseExists[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropf[j] = m_ropr[j];
                        if (m_ropf[j] > 0.0) {
                            for (size_t rp = 0; rp < nPhases(); rp++) {
                                if (m_rxnPhaseIsProduct[j][rp] && !m_phaseExists[rp]) {
                                    m_ropnet[j] = 0.0;
                                    m_ropf[j] = m_ropr[j] = 0.0;
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsProduct[j][p] && !m_phaseIsStable[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropf[j] = m_ropr[j];
                    }
                }
            }
        }
    }
    m_ROP_ok = true;
}

void InterfaceKinetics::getDeltaGibbs(double* deltaG)
{
    // Get the chemical potentials of the species in the all of the phases used
    // in the kinetics mechanism
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(m_mu.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_mu.data(), m_rbuf.data());
    if (deltaG != 0 && (m_rbuf.data() != deltaG)) {
        for (size_t j = 0; j < nReactions(); ++j) {
            deltaG[j] = m_rbuf[j];
        }
    }
}

void InterfaceKinetics::getDeltaElectrochemPotentials(double* deltaM)
{
    // Get the chemical potentials of the species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getElectrochemPotentials(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaM);
}

void InterfaceKinetics::getDeltaEnthalpy(double* deltaH)
{
    // Get the partial molar enthalpy of all species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void InterfaceKinetics::getDeltaEntropy(double* deltaS)
{
    // Get the partial molar entropy of all species in all of the phases
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void InterfaceKinetics::getDeltaSSGibbs(double* deltaGSS)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity We define these here as the
    // chemical potentials of the pure species at the temperature and pressure
    // of the solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_mu0.data(), deltaGSS);
}

void InterfaceKinetics::getDeltaSSEnthalpy(double* deltaH)
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

void InterfaceKinetics::getDeltaSSEntropy(double* deltaS)
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

bool InterfaceKinetics::addReaction(shared_ptr<Reaction> r_base, bool resize)
{
    if (!m_surf) {
        init();
    }

    size_t i = nReactions();
    shared_ptr<ReactionRate> rate = r_base->rate();
    if (rate) {
        rate->setContext(*r_base, *this);
    }

    bool added = Kinetics::addReaction(r_base, resize);
    if (!added) {
        return false;
    }

    m_rxnPhaseIsReactant.emplace_back(nPhases(), false);
    m_rxnPhaseIsProduct.emplace_back(nPhases(), false);

    for (const auto& [name, stoich] : r_base->reactants) {
        size_t k = kineticsSpeciesIndex(name);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    for (const auto& [name, stoich] : r_base->products) {
        size_t k = kineticsSpeciesIndex(name);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }

    // Set index of rate to number of reaction within kinetics
    rate->setRateIndex(nReactions() - 1);

    string rtype = rate->subType();
    if (rtype == "") {
        rtype = rate->type();
    }

    // If necessary, add new interface MultiRate evaluator
    if (m_rateTypes.find(rtype) == m_rateTypes.end()) {
        m_rateTypes[rtype] = m_rateHandlers.size();
        m_rateHandlers.push_back(rate->newMultiRate());
        m_rateHandlers.back()->resize(m_kk, nReactions(), nPhases());
    }

    // Add reaction rate to evaluator
    size_t index = m_rateTypes[rtype];
    m_rateHandlers[index]->add(nReactions() - 1, *rate);

    // Set flag for coverage dependence to true
    if (rate->compositionDependent()) {
        m_has_coverage_dependence = true;
    }

    // Set flag for electrochemistry to true
    if (r_base->usesElectrochemistry(*this)) {
        m_has_electrochemistry = true;
    }

    return true;
}

void InterfaceKinetics::modifyReaction(size_t i, shared_ptr<Reaction> r_base)
{
    Kinetics::modifyReaction(i, r_base);

    shared_ptr<ReactionRate> rate = r_base->rate();
    rate->setRateIndex(i);
    rate->setContext(*r_base, *this);

    string rtype = rate->subType();
    if (rtype == "") {
        rtype = rate->type();
    }

    // Ensure that interface MultiRate evaluator is available
    if (!m_rateTypes.count(rtype)) {
        throw CanteraError("InterfaceKinetics::modifyReaction",
            "Interface evaluator not available for type '{}'.", rtype);
    }
    // Replace reaction rate evaluator
    size_t index = m_rateTypes[rate->type()];
    m_rateHandlers[index]->replace(i, *rate);

    // Invalidate cached data
    m_redo_rates = true;
    m_temp += 0.1;
}

void InterfaceKinetics::setMultiplier(size_t i, double f)
{
    Kinetics::setMultiplier(i, f);
    m_ROP_ok = false;
}

void InterfaceKinetics::addThermo(shared_ptr<ThermoPhase> thermo)
{
    Kinetics::addThermo(thermo);
    m_phaseExists.push_back(true);
    m_phaseIsStable.push_back(true);
}

void InterfaceKinetics::init()
{
    if (thermo(0).nDim() > 2) {
        throw CanteraError("InterfaceKinetics::init", "no interface phase is present.");
    }
}

void InterfaceKinetics::resizeSpecies()
{
    size_t kOld = m_kk;
    Kinetics::resizeSpecies();
    if (m_kk != kOld && nReactions()) {
        throw CanteraError("InterfaceKinetics::resizeSpecies", "Cannot add"
            " species to InterfaceKinetics after reactions have been added.");
    }
    m_actConc.resize(m_kk);
    m_conc.resize(m_kk);
    m_mu0.resize(m_kk);
    m_mu.resize(m_kk);
    m_mu0_Kc.resize(m_kk);
    m_grt.resize(m_kk);
    m_phi.resize(nPhases(), 0.0);
}

void InterfaceKinetics::buildNetwork()
{
    for (auto& phase : m_thermo) {
        if (!phase->root()) {
            throw CanteraError("InterfaceKinetics::buildNetwork",
                "Phase '{}' is not attached to a Solution.", phase->name());
        }
    }
    vector<shared_ptr<ReactorBase>> reservoirs;
    for (size_t i = 1; i < nPhases(); i++) {
        auto r = newReservoir(thermo(i).root(), false);
        reservoirs.push_back(r);
    }
    auto rsurf = newReactorSurface(thermo(0).root(), reservoirs, false);
    m_integrator = new ReactorNet(rsurf);
}

void InterfaceKinetics::advanceCoverages(double tstep, double rtol, double atol,
        double maxStepSize, size_t maxSteps, size_t maxErrTestFails)
{
    // Stash the state of adjacent phases, and set their T and P to match the surface
    vector<vector<double>> savedStates(nPhases());
    for (size_t i = 1; i < nPhases(); i++) {
        savedStates[i].resize(thermo(i).partialStateSize());
        thermo(i).savePartialState(savedStates[i]);
        thermo(i).setState_TP(thermo(0).temperature(), thermo(0).pressure());
    }

    if (!m_integrator) {
        buildNetwork();
    }

    m_integrator->setTolerances(rtol, atol);
    m_integrator->setMaxTimeStep(maxStepSize);
    m_integrator->setMaxSteps(maxSteps);
    m_integrator->setMaxErrTestFails(maxErrTestFails);
    m_integrator->setInitialTime(0.0);
    m_integrator->advance(tstep);

    // Restore adjacent phases to their original states
    for (size_t i = 1; i < nPhases(); i++) {
        thermo(i).restorePartialState(savedStates[i]);
    }
}

void InterfaceKinetics::solvePseudoSteadyStateProblem(int loglevel)
{
    // Stash the state of adjacent phases, and set their T and P to match the surface
    vector<vector<double>> savedStates(nPhases());
    for (size_t i = 1; i < nPhases(); i++) {
        savedStates[i].resize(thermo(i).partialStateSize());
        thermo(i).savePartialState(savedStates[i]);
        thermo(i).setState_TP(thermo(0).temperature(), thermo(0).pressure());
    }

    if (!m_integrator) {
        buildNetwork();
    }

    m_integrator->setVerbose(loglevel != 0);
    m_integrator->solveSteady(loglevel);

    // Restore adjacent phases to their original states
    for (size_t i = 1; i < nPhases(); i++) {
        thermo(i).restorePartialState(savedStates[i]);
    }
}

void InterfaceKinetics::setPhaseExistence(const size_t iphase, const int exists)
{
    checkPhaseIndex(iphase);
    if (exists) {
        if (!m_phaseExists[iphase]) {
            m_phaseExistsCheck--;
            m_phaseExistsCheck = std::max(m_phaseExistsCheck, 0);
            m_phaseExists[iphase] = true;
        }
        m_phaseIsStable[iphase] = true;
    } else {
        if (m_phaseExists[iphase]) {
            m_phaseExistsCheck++;
            m_phaseExists[iphase] = false;
        }
        m_phaseIsStable[iphase] = false;
    }
}

int InterfaceKinetics::phaseExistence(const size_t iphase) const
{
    checkPhaseIndex(iphase);
    return m_phaseExists[iphase];
}

int InterfaceKinetics::phaseStability(const size_t iphase) const
{
    checkPhaseIndex(iphase);
    return m_phaseIsStable[iphase];
}

void InterfaceKinetics::setPhaseStability(const size_t iphase, const int isStable)
{
    checkPhaseIndex(iphase);
    if (isStable) {
        m_phaseIsStable[iphase] = true;
    } else {
        m_phaseIsStable[iphase] = false;
    }
}

double InterfaceKinetics::interfaceCurrent(const size_t iphase)
{
    vector<double> charges(m_kk, 0.0);
    vector<double> netProdRates(m_kk, 0.0);
    double dotProduct = 0.0;

    thermo(iphase).getCharges(charges.data());
    getNetProductionRates(netProdRates.data());

    for (size_t k = 0; k < thermo(iphase).nSpecies(); k++)
    {
        dotProduct += charges[k] * netProdRates[m_start[iphase] + k];
    }

    return dotProduct * Faraday;
}

Eigen::SparseMatrix<double> InterfaceKinetics::fwdRatesOfProgress_ddCi()
{
    // check derivatives are valid
    assertDerivativesValid("InterfaceKinetics::fwdRatesOfProgress_ddCi");
    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_reactantStoich, rop_rates);
}

Eigen::SparseMatrix<double> InterfaceKinetics::revRatesOfProgress_ddCi()
{
    // check derivatives are valid
    assertDerivativesValid("InterfaceKinetics::revRatesOfProgress_ddCi");
    // reverse reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    applyEquilibriumConstants(rop_rates.data());
    return calculateCompositionDerivatives(m_revProductStoich, rop_rates);
}

Eigen::SparseMatrix<double> InterfaceKinetics::netRatesOfProgress_ddCi()
{
    // check derivatives are valid
    assertDerivativesValid("InterfaceKinetics::netRatesOfProgress_ddCi");
    // forward reaction rate coefficients
    vector<double>& rop_rates = m_rbuf0;
    getFwdRateConstants(rop_rates.data());
    Eigen::SparseMatrix<double> jac = calculateCompositionDerivatives(m_reactantStoich,
        rop_rates);

    // reverse reaction rate coefficients
    applyEquilibriumConstants(rop_rates.data());
    return jac - calculateCompositionDerivatives(m_revProductStoich, rop_rates);
}

void InterfaceKinetics::setDerivativeSettings(const AnyMap& settings)
{
    bool force = settings.empty();
    if (force || settings.hasKey("skip-coverage-dependence")) {
        m_jac_skip_coverage_dependence = settings.getBool("skip-coverage-dependence",
            false);
    }
    if (force || settings.hasKey("skip-electrochemistry")) {
        m_jac_skip_electrochemistry = settings.getBool("skip-electrochemistry",
            false);
    }
    if (force || settings.hasKey("rtol-delta")) {
        m_jac_rtol_delta = settings.getDouble("rtol-delta", 1e-8);
    }
}

void InterfaceKinetics::getDerivativeSettings(AnyMap& settings) const
{
    settings["skip-coverage-dependence"] = m_jac_skip_coverage_dependence;
    settings["skip-electrochemistry"] = m_jac_skip_electrochemistry;
    settings["rtol-delta"] = m_jac_rtol_delta;
}

Eigen::SparseMatrix<double> InterfaceKinetics::calculateCompositionDerivatives(
    StoichManagerN& stoich, const vector<double>& in)
{
    vector<double>& outV = m_rbuf1;
    // derivatives handled by StoichManagerN
    copy(in.begin(), in.end(), outV.begin());
    return stoich.derivatives(m_actConc.data(), outV.data());
}

void InterfaceKinetics::assertDerivativesValid(const string& name)
{
    if (!m_jac_skip_coverage_dependence && m_has_coverage_dependence) {
        throw NotImplementedError(name, "Coverage-dependent reactions not supported.");
    } else if (!m_jac_skip_electrochemistry && m_has_electrochemistry) {
        throw NotImplementedError(name, "Electrochemical reactions not supported.");
    }
}

void InterfaceKinetics::applyEquilibriumConstants(double* rop)
{
    // For reverse rates computed from thermochemistry, multiply the forward
    // rate coefficients by the reciprocals of the equilibrium constants
    for (size_t i = 0; i < nReactions(); ++i) {
        rop[i] *= m_rkcn[i];
    }
}

}
