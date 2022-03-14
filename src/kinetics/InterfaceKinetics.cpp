/**
 *  @file InterfaceKinetics.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

InterfaceKinetics::InterfaceKinetics(ThermoPhase* thermo) :
    m_redo_rates(false),
    m_surf(0),
    m_integrator(0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_logtemp(0.0),
    m_has_coverage_dependence(false),
    m_has_electrochem_rxns(false),
    m_has_exchange_current_density_formulation(false),
    m_phaseExistsCheck(false),
    m_ioFlag(0),
    m_nDim(2)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
}

InterfaceKinetics::~InterfaceKinetics()
{
    delete m_integrator;
}

void InterfaceKinetics::resizeReactions()
{
    Kinetics::resizeReactions();

    for (auto& rates : m_interfaceRates) {
        rates->resize(nTotalSpecies(), nReactions(), nPhases());
        // @todo ensure that ReactionData are updated; calling rates->update
        //      blocks correct behavior in InterfaceKinetics::_update_rates_T
        //      and running updateROP() is premature
    }
}

void InterfaceKinetics::setElectricPotential(int n, doublereal V)
{
    thermo(n).setElectricPotential(V);
    m_redo_rates = true;
}

void InterfaceKinetics::_update_rates_T()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();
    if (m_has_coverage_dependence) {
        m_surf->getCoverages(m_actConc.data());
        m_rates.update_C(m_actConc.data());
        m_redo_rates = true;
    }

    // Go find the temperature from the surface
    doublereal T = thermo(surfacePhaseIndex()).temperature();
    m_redo_rates = true;
    if (T != m_temp || m_redo_rates) {
        m_logtemp = log(T);

        //  Calculate the forward rate constant by calling m_rates and store it in m_rfn[]
        m_rates.update(T, m_logtemp, m_rfn.data());
        for (size_t n = 0; n < nPhases(); n++) {
            thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
        }

        // Use the stoichiometric manager to find deltaH for each reaction.
        getReactionDelta(m_grt.data(), m_dH.data());
        applyStickingCorrection(T, m_rfn.data());

        // If we need to do conversions between exchange current density
        // formulation and regular formulation (either way) do it here.
        if (m_has_exchange_current_density_formulation) {
            convertExchangeCurrentDensityFormulation(m_rfn.data());
        }
        if (m_has_electrochem_rxns) {
            applyVoltageKfwdCorrection(m_rfn.data());
        }
        m_temp = T;
        m_ROP_ok = false;
        m_redo_rates = false;
    }

    // loop over interface MultiRate evaluators for each reaction type
    for (auto& rates : m_interfaceRates) {
        bool changed = rates->update(thermo(surfacePhaseIndex()), *this);
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
        const ThermoPhase* tp = m_thermo[n];
        /*
         * We call the getActivityConcentrations function of each ThermoPhase
         * class that makes up this kinetics object to obtain the generalized
         * concentrations for species within that class. This is collected in
         * the vector m_conc. m_start[] are integer indices for that vector
         * denoting the start of the species for each phase.
         */
        tp->getActivityConcentrations(m_actConc.data() + m_start[n]);

        // Get regular concentrations too
        tp->getConcentrations(m_conc.data() + m_start[n]);
    }
    m_ROP_ok = false;
}

void InterfaceKinetics::getActivityConcentrations(doublereal* const conc)
{
    _update_rates_C();
    copy(m_actConc.begin(), m_actConc.end(), conc);
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
        doublereal rrt = 1.0 / thermo(reactionPhaseIndex()).RT();

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
    updateExchangeCurrentQuantities();
    size_t ik = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
        for (size_t k = 0; k < thermo(n).nSpecies(); k++) {
            m_mu0_Kc[ik] = m_mu0[ik] + Faraday * m_phi[n] * thermo(n).charge(k);
            m_mu0_Kc[ik] -= thermo(reactionPhaseIndex()).RT()
                            * thermo(n).logStandardConc(k);
            ik++;
        }
    }
}

void InterfaceKinetics::getEquilibriumConstants(doublereal* kc)
{
    updateMu0();
    doublereal rrt = 1.0 / thermo(reactionPhaseIndex()).RT();
    std::fill(kc, kc + nReactions(), 0.0);
    getReactionDelta(m_mu0_Kc.data(), kc);
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-kc[i]*rrt);
    }
}

void InterfaceKinetics::updateExchangeCurrentQuantities()
{
    // Calculate:
    //   - m_StandardConc[]
    //   - m_ProdStanConcReac[]
    //   - m_deltaG0[]
    //   - m_mu0[]

    // First collect vectors of the standard Gibbs free energies of the
    // species and the standard concentrations
    //   - m_mu0
    //   - m_StandardConc
    size_t ik = 0;

    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_StandardConc[ik] = thermo(n).standardConcentration(k);
            ik++;
        }
    }

    getReactionDelta(m_mu0.data(), m_deltaG0.data());

    //  Calculate the product of the standard concentrations of the reactants
    for (size_t i = 0; i < nReactions(); i++) {
        m_ProdStanConcReac[i] = 1.0;
    }
    m_reactantStoich.multiply(m_StandardConc.data(), m_ProdStanConcReac.data());
}

void InterfaceKinetics::applyVoltageKfwdCorrection(doublereal* const kf)
{
    // Compute the electrical potential energy of each species
    size_t ik = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_pot[ik] = Faraday * thermo(n).charge(k) * m_phi[n];
            ik++;
        }
    }

    // Compute the change in electrical potential energy for each reaction. This
    // will only be non-zero if a potential difference is present.
    getReactionDelta(m_pot.data(), deltaElectricEnergy_.data());

    // Modify the reaction rates. Only modify those with a non-zero activation
    // energy. Below we decrease the activation energy below zero but in some
    // debug modes we print out a warning message about this.

    // NOTE, there is some discussion about this point. Should we decrease the
    // activation energy below zero? I don't think this has been decided in any
    // definitive way. The treatment below is numerically more stable, however.
    for (size_t i = 0; i < m_beta.size(); i++) {
        size_t irxn = m_ctrxn[i];

        // Add the voltage correction to the forward reaction rate constants.
        double eamod = m_beta[i] * deltaElectricEnergy_[irxn];
        if (eamod != 0.0) {
            kf[irxn] *= exp(-eamod/thermo(reactionPhaseIndex()).RT());
        }
    }
}

void InterfaceKinetics::convertExchangeCurrentDensityFormulation(doublereal* const kfwd)
{
    updateExchangeCurrentQuantities();
    // Loop over all reactions which are defined to have a voltage transfer
    // coefficient that affects the activity energy for the reaction
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        size_t irxn = m_ctrxn[i];

        // Determine whether the reaction rate constant is in an exchange
        // current density formulation format.
        int iECDFormulation = m_ctrxn_ecdf[i];
        if (iECDFormulation) {
            // We need to have the straight chemical reaction rate constant to
            // come out of this calculation.
            double tmp = exp(- m_beta[i] * m_deltaG0[irxn]
                                / thermo(reactionPhaseIndex()).RT());
            tmp *= 1.0 / m_ProdStanConcReac[irxn] / Faraday;
            kfwd[irxn] *= tmp;
        }
    }
}

void InterfaceKinetics::getFwdRateConstants(doublereal* kfwd)
{
    updateROP();
    for (size_t i = 0; i < nReactions(); i++) {
        // base rate coefficient multiplied by perturbation factor
        kfwd[i] = m_rfn[i] * m_perturb[i];
    }
}

void InterfaceKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
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

void InterfaceKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the all of the phases used
    // in the kinetics mechanism
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(m_mu.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_mu.data(), m_deltaG.data());
    if (deltaG != 0 && (m_deltaG.data() != deltaG)) {
        for (size_t j = 0; j < nReactions(); ++j) {
            deltaG[j] = m_deltaG[j];
        }
    }
}

void InterfaceKinetics::getDeltaElectrochemPotentials(doublereal* deltaM)
{
    // Get the chemical potentials of the species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getElectrochemPotentials(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaM);
}

void InterfaceKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void InterfaceKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in all of the phases
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void InterfaceKinetics::getDeltaSSGibbs(doublereal* deltaGSS)
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

void InterfaceKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    // Get the standard state enthalpies of the species. This is the array of
    // chemical potentials at unit activity We define these here as the
    // enthalpies of the pure species at the temperature and pressure of the
    // solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEnthalpy_RT(m_grt.data() + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= thermo(reactionPhaseIndex()).RT();
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void InterfaceKinetics::getDeltaSSEntropy(doublereal* deltaS)
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
    bool added = Kinetics::addReaction(r_base, resize);
    if (!added) {
        return false;
    }

    if (r_base->reversible) {
        m_revindex.push_back(i);
    } else {
        m_irrev.push_back(i);
    }

    m_rxnPhaseIsReactant.emplace_back(nPhases(), false);
    m_rxnPhaseIsProduct.emplace_back(nPhases(), false);

    for (const auto& sp : r_base->reactants) {
        size_t k = kineticsSpeciesIndex(sp.first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    for (const auto& sp : r_base->products) {
        size_t k = kineticsSpeciesIndex(sp.first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }

    if (!(r_base->usesLegacy())) {
        // Set index of rate to number of reaction within kinetics
        shared_ptr<ReactionRate> rate = r_base->rate();
        rate->setRateIndex(nReactions() - 1);
        rate->setContext(*r_base, *this);

        // If necessary, add new interface MultiRate evaluator
        if (m_interfaceTypes.find(rate->type()) == m_interfaceTypes.end()) {
            m_interfaceTypes[rate->type()] = m_interfaceRates.size();
            m_interfaceRates.push_back(rate->newMultiRate());
            m_interfaceRates.back()->resize(m_kk, nReactions(), nPhases());
        }

        // Add reaction rate to evaluator
        size_t index = m_interfaceTypes[rate->type()];
        m_interfaceRates[index]->add(nReactions() - 1, *rate);

    } else if (r_base->reaction_type == SURFACE_RXN) {
        InterfaceReaction2& r = dynamic_cast<InterfaceReaction2&>(*r_base);
        SurfaceArrhenius rate = buildSurfaceArrhenius(i, r, false);
        m_rates.install(i, rate);

        // Turn on the global flag indicating surface coverage dependence
        if (!r.coverage_deps.empty()) {
            m_has_coverage_dependence = true;
        }
        ElectrochemicalReaction2* re = dynamic_cast<ElectrochemicalReaction2*>(&r);
        if (re) {
            m_has_electrochem_rxns = true;
            m_beta.push_back(re->beta);
            m_ctrxn.push_back(i);
            if (re->exchange_current_density_formulation) {
                m_has_exchange_current_density_formulation = true;
                m_ctrxn_ecdf.push_back(1);
            } else {
                m_ctrxn_ecdf.push_back(0);
            }
        }
    } else {
        throw NotImplementedError("InterfaceKinetics::addReaction");
    }
    deltaElectricEnergy_.push_back(0.0);
    m_deltaG0.push_back(0.0);
    m_deltaG.push_back(0.0);
    m_ProdStanConcReac.push_back(0.0);

    return true;
}

void InterfaceKinetics::modifyReaction(size_t i, shared_ptr<Reaction> r_base)
{
    Kinetics::modifyReaction(i, r_base);
    if (!(r_base->usesLegacy())) {
        shared_ptr<ReactionRate> rate = r_base->rate();
        rate->setRateIndex(i);
        rate->setContext(*r_base, *this);

        const auto& rtype = rate->type();
        // Ensure that interface MultiRate evaluator is available
        if (!m_interfaceTypes.count(rtype)) {
            throw CanteraError("InterfaceKinetics::modifyReaction",
                "Interface evaluator not available for type '{}'.", rtype);
        }
        // Replace reaction rate evaluator
        size_t index = m_interfaceTypes[rate->type()];
        m_interfaceRates[index]->replace(i, *rate);

    } else if (r_base->reaction_type == SURFACE_RXN) {
        InterfaceReaction2& r = dynamic_cast<InterfaceReaction2&>(*r_base);
        SurfaceArrhenius rate = buildSurfaceArrhenius(i, r, true);
        m_rates.replace(i, rate);
    } else {
        throw NotImplementedError("InterfaceKinetics::modifyReaction");
    }
    // Invalidate cached data
    m_redo_rates = true;
    m_temp += 0.1;
}

SurfaceArrhenius InterfaceKinetics::buildSurfaceArrhenius(
    size_t i, InterfaceReaction2& r, bool replace)
{
    if (r.is_sticking_coefficient) {
        // Identify the interface phase
        size_t iInterface = npos;
        size_t min_dim = 4;
        for (size_t n = 0; n < nPhases(); n++) {
            if (thermo(n).nDim() < min_dim) {
                iInterface = n;
                min_dim = thermo(n).nDim();
            }
        }

        std::string sticking_species = r.sticking_species;
        if (sticking_species == "") {
            // Identify the sticking species if not explicitly given
            bool foundStick = false;
            for (const auto& sp : r.reactants) {
                size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(sp.first));
                if (iPhase != iInterface) {
                    // Non-interface species. There should be exactly one of these
                    if (foundStick) {
                        throw InputFileError("InterfaceKinetics::buildSurfaceArrhenius",
                            r.input, "Multiple non-interface species ('{}' and '{}')\n"
                            "found in sticking reaction: '{}'.\nSticking species "
                            "must be explicitly specified.",
                            sticking_species, sp.first, r.equation());
                    }
                    foundStick = true;
                    sticking_species = sp.first;
                }
            }
            if (!foundStick) {
                throw InputFileError("InterfaceKinetics::buildSurfaceArrhenius",
                    r.input, "No non-interface species found "
                    "in sticking reaction: '{}'", r.equation());
            }
        }

        double surface_order = 0.0;
        double multiplier = 1.0;
        // Adjust the A-factor
        for (const auto& sp : r.reactants) {
            size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(sp.first));
            const ThermoPhase& p = thermo(iPhase);
            size_t k = p.speciesIndex(sp.first);
            if (sp.first == sticking_species) {
                multiplier *= sqrt(GasConstant/(2*Pi*p.molecularWeight(k)));
            } else {
                // Non-sticking species. Convert from coverages used in the
                // sticking probability expression to the concentration units
                // used in the mass action rate expression. For surface phases,
                // the dependence on the site density is incorporated when the
                // rate constant is evaluated, since we don't assume that the
                // site density is known at this time.
                double order = getValue(r.orders, sp.first, sp.second);
                if (&p == m_surf) {
                    multiplier *= pow(m_surf->size(k), order);
                    surface_order += order;
                } else {
                    multiplier *= pow(p.standardConcentration(k), -order);
                }
            }
        }

        if (!replace) {
            m_stickingData.emplace_back(StickData{i, surface_order, multiplier,
                                                  r.use_motz_wise_correction});
        } else {
            // Modifying an existing sticking reaction.
            for (auto& item : m_stickingData) {
                if (item.index == i) {
                    item.order = surface_order;
                    item.multiplier = multiplier;
                    item.use_motz_wise = r.use_motz_wise_correction;
                    break;
                }
            }
        }
    }

    SurfaceArrhenius rate(r.rate.preExponentialFactor(),
                          r.rate.temperatureExponent(),
                          r.rate.activationEnergy_R());

    // Set up coverage dependencies
    for (const auto& sp : r.coverage_deps) {
        size_t k = thermo(reactionPhaseIndex()).speciesIndex(sp.first);
        rate.addCoverageDependence(k, sp.second.a, sp.second.m, sp.second.E);
    }
    return rate;
}

void InterfaceKinetics::setIOFlag(int ioFlag)
{
    m_ioFlag = ioFlag;
    if (m_integrator) {
        m_integrator->setIOFlag(ioFlag);
    }
}

void InterfaceKinetics::addPhase(ThermoPhase& thermo)
{
    Kinetics::addPhase(thermo);
    m_phaseExists.push_back(true);
    m_phaseIsStable.push_back(true);
}

void InterfaceKinetics::init()
{
    size_t ks = reactionPhaseIndex();
    if (ks == npos) {
        throw CanteraError("InterfaceKinetics::init",
                           "no surface phase is present.");
    }

    // Check to see that the interface routine has a dimension of 2
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != m_nDim) {
        throw CanteraError("InterfaceKinetics::init",
                           "expected interface dimension = 2, but got dimension = {}",
                           m_surf->nDim());
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
    m_StandardConc.resize(m_kk, 0.0);
    m_mu0.resize(m_kk);
    m_mu.resize(m_kk);
    m_mu0_Kc.resize(m_kk);
    m_grt.resize(m_kk);
    m_pot.resize(m_kk, 0.0);
    m_phi.resize(nPhases(), 0.0);
}

doublereal InterfaceKinetics::electrochem_beta(size_t irxn) const
{
    warn_deprecated("InterfaceKinetics::electrochem_beta",
                    "This function only works for the legacy framework. "
                    "To be removed after Cantera 2.6.");

    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        if (m_ctrxn[i] == irxn) {
            return m_beta[i];
        }
    }
    return 0.0;
}

void InterfaceKinetics::advanceCoverages(doublereal tstep, doublereal rtol,
                                         doublereal atol, doublereal maxStepSize,
                                         size_t maxSteps, size_t maxErrTestFails)
{
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k{this};
        m_integrator = new ImplicitSurfChem(k);
    }
    m_integrator->setTolerances(rtol, atol);
    m_integrator->setMaxStepSize(maxStepSize);
    m_integrator->setMaxSteps(maxSteps);
    m_integrator->setMaxErrTestFails(maxErrTestFails);
    m_integrator->integrate(0.0, tstep);
    delete m_integrator;
    m_integrator = 0;
}

void InterfaceKinetics::solvePseudoSteadyStateProblem(
    int ifuncOverride, doublereal timeScaleOverride)
{
    // create our own solver object
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k{this};
        m_integrator = new ImplicitSurfChem(k);
        m_integrator->initialize();
    }
    m_integrator->setIOFlag(m_ioFlag);
    // New direct method to go here
    m_integrator->solvePseudoSteadyStateProblem(ifuncOverride, timeScaleOverride);
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

void InterfaceKinetics::applyStickingCorrection(double T, double* kf)
{
    if (m_stickingData.empty()) {
        return;
    }

    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    vector_fp& factors = cached.value;

    double n0 = m_surf->siteDensity();
    if (!cached.validate(n0)) {
        factors.resize(m_stickingData.size());
        for (size_t n = 0; n < m_stickingData.size(); n++) {
            factors[n] = pow(n0, -m_stickingData[n].order);
        }
    }

    for (size_t n = 0; n < m_stickingData.size(); n++) {
        const StickData& item = m_stickingData[n];
        if (item.use_motz_wise) {
            kf[item.index] /= 1 - 0.5 * kf[item.index];
        }
        kf[item.index] *= factors[n] * sqrt(T) * item.multiplier;
    }
}

}
