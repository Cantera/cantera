/**
 *  @file InterfaceKinetics.cpp
 */

// Copyright 2002  California Institute of Technology

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/thermo/SurfPhase.h"

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/RateCoeffMgr.h"

#include "cantera/kinetics/ImplicitSurfChem.h"

using namespace std;

namespace Cantera
{
InterfaceKinetics::InterfaceKinetics(thermo_t* thermo) :
    Kinetics(),
    m_redo_rates(false),
    m_nirrev(0),
    m_nrev(0),
    m_surf(0),
    m_integrator(0),
    m_beta(0),
    m_ctrxn(0),
    m_ctrxn_ecdf(0),
    m_StandardConc(0),
    m_deltaG0(0),
    m_ProdStanConcReac(0),
    m_logp0(0.0),
    m_logc0(0.0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_logtemp(0.0),
    m_finalized(false),
    m_has_coverage_dependence(false),
    m_has_electrochem_rxns(false),
    m_has_exchange_current_density_formulation(false),
    m_phaseExistsCheck(false),
    m_phaseExists(0),
    m_phaseIsStable(0),
    m_rxnPhaseIsReactant(0),
    m_rxnPhaseIsProduct(0),
    m_ioFlag(0)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
}

InterfaceKinetics::~InterfaceKinetics()
{
    delete m_integrator;
}

InterfaceKinetics::InterfaceKinetics(const InterfaceKinetics& right) :
    Kinetics(),
    m_redo_rates(false),
    m_nirrev(0),
    m_nrev(0),
    m_surf(0),
    m_integrator(0),
    m_beta(0),
    m_ctrxn(0),
    m_ctrxn_ecdf(0),
    m_StandardConc(0),
    m_deltaG0(0),
    m_ProdStanConcReac(0),
    m_logp0(0.0),
    m_logc0(0.0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_logtemp(0.0),
    m_finalized(false),
    m_has_coverage_dependence(false),
    m_has_electrochem_rxns(false),
    m_has_exchange_current_density_formulation(false),
    m_phaseExistsCheck(false),
    m_phaseExists(0),
    m_phaseIsStable(0),
    m_rxnPhaseIsReactant(0),
    m_rxnPhaseIsProduct(0),
    m_ioFlag(0)
{
    /*
     * Call the assignment operator
     */
    *this = operator=(right);
}

InterfaceKinetics& InterfaceKinetics::
operator=(const InterfaceKinetics& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Kinetics::operator=(right);

    m_grt                  = right.m_grt;
    m_revindex             = right.m_revindex;
    m_rates                = right.m_rates;
    m_redo_rates           = right.m_redo_rates;
    m_index                = right.m_index;
    m_irrev                = right.m_irrev;
    m_rxnstoich            = right.m_rxnstoich;
    m_nirrev               = right.m_nirrev;
    m_nrev                 = right.m_nrev;
    m_rrxn                 = right.m_rrxn;
    m_prxn                 = right.m_prxn;
    m_rxneqn               = right.m_rxneqn;
    m_conc                 = right.m_conc;
    m_mu0                  = right.m_mu0;
    m_phi                  = right.m_phi;
    m_pot                  = right.m_pot;
    m_rwork                = right.m_rwork;
    m_E                    = right.m_E;
    m_surf                 = right.m_surf;  //DANGER - shallow copy
    m_integrator           = right.m_integrator;  //DANGER - shallow copy
    m_beta                 = right.m_beta;
    m_ctrxn                = right.m_ctrxn;
    m_ctrxn_ecdf           = right.m_ctrxn_ecdf;
    m_StandardConc         = right.m_StandardConc;
    m_deltaG0              = right.m_deltaG0;
    m_ProdStanConcReac     = right.m_ProdStanConcReac;
    m_logp0                = right.m_logp0;
    m_logc0                = right.m_logc0;
    m_ropf                 = right.m_ropf;
    m_ropr                 = right.m_ropr;
    m_ropnet               = right.m_ropnet;
    m_ROP_ok               = right.m_ROP_ok;
    m_temp                 = right.m_temp;
    m_logtemp              = right.m_logtemp;
    m_rfn                  = right.m_rfn;
    m_rkcn                 = right.m_rkcn;
    m_finalized            = right.m_finalized;
    m_has_coverage_dependence = right.m_has_coverage_dependence;
    m_has_electrochem_rxns = right.m_has_electrochem_rxns;
    m_has_exchange_current_density_formulation = right.m_has_exchange_current_density_formulation;
    m_phaseExistsCheck     = right.m_phaseExistsCheck;
    m_phaseExists          = right.m_phaseExists;
    m_phaseIsStable        = right.m_phaseIsStable;
    m_rxnPhaseIsReactant   = right.m_rxnPhaseIsReactant;
    m_rxnPhaseIsProduct    = right.m_rxnPhaseIsProduct;
    m_ioFlag               = right.m_ioFlag;

    return *this;
}

int InterfaceKinetics::type() const
{
    return cInterfaceKinetics;
}

Kinetics* InterfaceKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    InterfaceKinetics* iK = new InterfaceKinetics(*this);
    iK->assignShallowPointers(tpVector);
    return iK;
}

void InterfaceKinetics::setElectricPotential(int n, doublereal V)
{
    thermo(n).setElectricPotential(V);
    m_redo_rates = true;
}

void InterfaceKinetics::_update_rates_T()
{
    _update_rates_phi();
    if (m_has_coverage_dependence) {
        m_surf->getCoverages(DATA_PTR(m_conc));
        m_rates.update_C(DATA_PTR(m_conc));
        m_redo_rates = true;
    }
    doublereal T = thermo(surfacePhaseIndex()).temperature();
    m_redo_rates = true;
    if (T != m_temp || m_redo_rates) {
        m_logtemp = log(T);
        m_rates.update(T, m_logtemp, DATA_PTR(m_rfn));
        if (m_has_exchange_current_density_formulation) {
            applyExchangeCurrentDensityFormulation(DATA_PTR(m_rfn));
        }
        if (m_has_electrochem_rxns) {
            applyButlerVolmerCorrection(DATA_PTR(m_rfn));
        }
        m_temp = T;
        updateKc();
        m_ROP_ok = false;
        m_redo_rates = false;
    }
}

void InterfaceKinetics::_update_rates_phi()
{
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
        /*
         * We call the getActivityConcentrations function of each
         * ThermoPhase class that makes up this kinetics object to
         * obtain the generalized concentrations for species within that
         * class. This is collected in the vector m_conc. m_start[]
         * are integer indices for that vector denoting the start of the
         * species for each phase.
         */
        thermo(n).getActivityConcentrations(DATA_PTR(m_conc) + m_start[n]);
    }
    m_ROP_ok = false;
}

void InterfaceKinetics::getActivityConcentrations(doublereal* const conc)
{
    _update_rates_C();
    copy(m_conc.begin(), m_conc.end(), conc);
}

void InterfaceKinetics::updateKc()
{
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    //static vector_fp mu(nTotalSpecies());
    if (m_nrev > 0) {
        /*
         * Get the vector of standard state electrochemical potentials for species in the Interfacial
         * kinetics object and store it in m_mu0[]
         */
        size_t nsp, ik = 0;
        doublereal rt = GasConstant*thermo(0).temperature();
        doublereal rrt = 1.0 / rt;
        size_t np = nPhases();
        for (size_t n = 0; n < np; n++) {
            thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                m_mu0[ik] -= rt * thermo(n).logStandardConc(k);
                m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                ik++;
            }
        }

        // compute Delta mu^0 for all reversible reactions
        m_rxnstoich.getRevReactionDelta(m_ii, DATA_PTR(m_mu0), DATA_PTR(m_rkcn));

        for (size_t i = 0; i < m_nrev; i++) {
            size_t irxn = m_revindex[i];
            if (irxn == npos || irxn >= nReactions()) {
                throw CanteraError("InterfaceKinetics",
                                   "illegal value: irxn = "+int2str(irxn));
            }
            // WARNING this may overflow HKM
            m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
        }
        for (size_t i = 0; i != m_nirrev; ++i) {
            m_rkcn[ m_irrev[i] ] = 0.0;
        }
    }
}

void InterfaceKinetics::checkPartialEquil()
{
    vector_fp dmu(nTotalSpecies(), 0.0);
    vector_fp rmu(std::max<size_t>(nReactions(), 1), 0.0);
    if (m_nrev > 0) {
        doublereal rt = GasConstant*thermo(0).temperature();
        cout << "T = " << thermo(0).temperature() << " " << rt << endl;
        size_t nsp, ik=0;
        //doublereal rt = GasConstant*thermo(0).temperature();
        //            doublereal rrt = 1.0/rt;
        doublereal delta;
        for (size_t n = 0; n < nPhases(); n++) {
            thermo(n).getChemPotentials(DATA_PTR(dmu) + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                delta = Faraday * m_phi[n] * thermo(n).charge(k);
                //cout << thermo(n).speciesName(k) << "   " << (delta+dmu[ik])/rt << " " << dmu[ik]/rt << endl;
                dmu[ik] += delta;
                ik++;
            }
        }

        // compute Delta mu^ for all reversible reactions
        m_rxnstoich.getRevReactionDelta(m_ii, DATA_PTR(dmu), DATA_PTR(rmu));
        updateROP();
        for (size_t i = 0; i < m_nrev; i++) {
            size_t irxn = m_revindex[i];
            cout << "Reaction " << reactionString(irxn)
                 << "  " << rmu[irxn]/rt << endl;
            printf("%12.6e  %12.6e  %12.6e  %12.6e \n",
                   m_ropf[irxn], m_ropr[irxn], m_ropnet[irxn],
                   m_ropnet[irxn]/(m_ropf[irxn] + m_ropr[irxn]));
        }
    }
}

void InterfaceKinetics::getFwdRatesOfProgress(doublereal* fwdROP)
{
    updateROP();
    std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
}

void InterfaceKinetics::getRevRatesOfProgress(doublereal* revROP)
{
    updateROP();
    std::copy(m_ropr.begin(), m_ropr.end(), revROP);
}

void InterfaceKinetics::getNetRatesOfProgress(doublereal* netROP)
{
    updateROP();
    std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
}

void InterfaceKinetics::getEquilibriumConstants(doublereal* kc)
{
    size_t ik=0;
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
            m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
            ik++;
        }
    }

    fill(kc, kc + m_ii, 0.0);

    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_mu0), kc);

    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-kc[i]*rrt);
    }
}

void InterfaceKinetics::getExchangeCurrentQuantities()
{
    /*
     * First collect vectors of the standard Gibbs free energies of the
     * species and the standard concentrations
     *   - m_mu0
     *   - m_logStandardConc
     */
    size_t ik = 0;

    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_StandardConc[ik] = thermo(n).standardConcentration(k);
            ik++;
        }
    }

    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_mu0), DATA_PTR(m_deltaG0));

    for (size_t i = 0; i < m_ii; i++) {
        m_ProdStanConcReac[i] = 1.0;
    }

    m_rxnstoich.multiplyReactants(DATA_PTR(m_StandardConc), DATA_PTR(m_ProdStanConcReac));

}

void InterfaceKinetics::getCreationRates(doublereal* cdot)
{
    updateROP();
    m_rxnstoich.getCreationRates(m_kk, &m_ropf[0], &m_ropr[0], cdot);
}

void InterfaceKinetics::getDestructionRates(doublereal* ddot)
{
    updateROP();
    m_rxnstoich.getDestructionRates(m_kk, &m_ropf[0], &m_ropr[0], ddot);
}

void InterfaceKinetics::getNetProductionRates(doublereal* net)
{
    updateROP();
    m_rxnstoich.getNetProductionRates(m_kk, &m_ropnet[0], net);
}

void InterfaceKinetics::applyButlerVolmerCorrection(doublereal* const kf)
{
    // compute the electrical potential energy of each species
    size_t ik = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_pot[ik] = Faraday*thermo(n).charge(k)*m_phi[n];
            ik++;
        }
    }

    // Compute the change in electrical potential energy for each
    // reaction. This will only be non-zero if a potential
    // difference is present.
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_pot), DATA_PTR(m_rwork));

    // Modify the reaction rates. Only modify those with a
    // non-zero activation energy. Below we decrease the
    // activation energy below zero but in some debug modes
    // we print out a warning message about this.
    /*
     *   NOTE, there is some discussion about this point.
     *   Should we decrease the activation energy below zero?
     *   I don't think this has been decided in any definitive way.
     *   The treatment below is numerically more stable, however.
     */
    doublereal eamod;
#ifdef DEBUG_KIN_MODE
    doublereal ea;
#endif
    for (size_t i = 0; i < m_beta.size(); i++) {
        size_t irxn = m_ctrxn[i];
        eamod = m_beta[i]*m_rwork[irxn];
        //  if (eamod != 0.0 && m_E[irxn] != 0.0) {
        if (eamod != 0.0) {
#ifdef DEBUG_KIN_MODE
            ea = GasConstant * m_E[irxn];
            if (eamod + ea < 0.0) {
                writelog("Warning: act energy mod too large!\n");
                writelog("  Delta phi = "+fp2str(m_rwork[irxn]/Faraday)+"\n");
                writelog("  Delta Ea = "+fp2str(eamod)+"\n");
                writelog("  Ea = "+fp2str(ea)+"\n");
                for (n = 0; n < np; n++) {
                    writelog("Phase "+int2str(n)+": phi = "
                             +fp2str(m_phi[n])+"\n");
                }
            }
#endif
            doublereal rt = GasConstant*thermo(0).temperature();
            doublereal rrt = 1.0/rt;
            kf[irxn] *= exp(-eamod*rrt);
        }
    }
}

void InterfaceKinetics::applyExchangeCurrentDensityFormulation(doublereal* const kfwd)
{
    getExchangeCurrentQuantities();
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        size_t irxn = m_ctrxn[i];
        int iECDFormulation =  m_ctrxn_ecdf[i];
        if (iECDFormulation) {
            double tmp = exp(- m_beta[i] * m_deltaG0[irxn] * rrt);
            double tmp2 = m_ProdStanConcReac[irxn];
            tmp *= 1.0  / tmp2 / Faraday;
            kfwd[irxn] *= tmp;
        }
    }
}

void InterfaceKinetics::getFwdRateConstants(doublereal* kfwd)
{

    updateROP();

    // copy rate coefficients into kfwd
    copy(m_rfn.begin(), m_rfn.end(), kfwd);

    // multiply by perturbation factor
    multiply_each(kfwd, kfwd + nReactions(), m_perturb.begin());

}

void InterfaceKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
{
    getFwdRateConstants(krev);
    if (doIrreversible) {
        getEquilibriumConstants(&m_ropnet[0]);
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] /= m_ropnet[i];
        }
    } else {
        multiply_each(krev, krev + nReactions(), m_rkcn.begin());
    }
}

void InterfaceKinetics::getActivationEnergies(doublereal* E)
{
    copy(m_E.begin(), m_E.end(), E);
}

void InterfaceKinetics::updateROP()
{
    _update_rates_T();
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // for reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by concentration products
    m_rxnstoich.multiplyReactants(DATA_PTR(m_conc), DATA_PTR(m_ropf));
    //m_reactantStoich.multiply(m_conc.begin(), ropf.begin());

    // for reversible reactions, multiply ropr by concentration
    // products
    m_rxnstoich.multiplyRevProducts(DATA_PTR(m_conc),
                                    DATA_PTR(m_ropr));
    //m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

    // do global reactions
    //m_globalReactantStoich.power(m_conc.begin(), ropf.begin());

    for (size_t j = 0; j != m_ii; ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    /*
     *  For reactions involving multiple phases, we must check that the phase
     *  being consumed actually exists. This is particularly important for
     *  phases that are stoichiometric phases containing one species with a unity activity
     */
    if (m_phaseExistsCheck) {
        for (size_t j = 0; j != m_ii; ++j) {
            if ((m_ropr[j] >  m_ropf[j]) && (m_ropr[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsProduct[j][p]) {
                        if (! m_phaseExists[p]) {
                            m_ropnet[j] = 0.0;
                            m_ropr[j] = m_ropf[j];
                            if (m_ropf[j] > 0.0) {
                                for (size_t rp = 0; rp < nPhases(); rp++) {
                                    if (m_rxnPhaseIsReactant[j][rp]) {
                                        if (! m_phaseExists[rp]) {
                                            m_ropnet[j] = 0.0;
                                            m_ropr[j] = m_ropf[j] = 0.0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsReactant[j][p]) {
                        if (! m_phaseIsStable[p]) {
                            m_ropnet[j] = 0.0;
                            m_ropr[j] = m_ropf[j];
                        }
                    }
                }
            } else if ((m_ropf[j] > m_ropr[j]) && (m_ropf[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsReactant[j][p]) {
                        if (! m_phaseExists[p]) {
                            m_ropnet[j] = 0.0;
                            m_ropf[j] = m_ropr[j];
                            if (m_ropf[j] > 0.0) {
                                for (size_t rp = 0; rp < nPhases(); rp++) {
                                    if (m_rxnPhaseIsProduct[j][rp]) {
                                        if (! m_phaseExists[rp]) {
                                            m_ropnet[j] = 0.0;
                                            m_ropf[j] = m_ropr[j] = 0.0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsProduct[j][p]) {
                        if (! m_phaseIsStable[p]) {
                            m_ropnet[j] = 0.0;
                            m_ropf[j] = m_ropr[j];
                        }
                    }
                }
            }
        }
    }

    m_ROP_ok = true;
}

#ifdef KINETICS_WITH_INTERMEDIATE_ZEROED_PHASES
//=================================================================================================
InterfaceKinetics::adjustRatesForIntermediatePhases()
{
    doublereal sFac = 1.0;

    getCreatingRates(DATA_PTR(m_speciestmpP));
    getDestructionRates(DATA_PTR(m_speciestmpD));

    for (iphase = 0; iphase < nphases; iphase++) {
        if (m_intermediatePhases(iphase)) {
            for (isp = 0; isp < nspecies; isp++) {
                if (m_speciesTmpD[ispI] > m_speciesTmpP[I]) {
                    sFac = m_speciesTmpD[ispI]/ m_speciesTmpP[I];
                }
                // Loop over reactions that are reactants for the species in the phase
                // reducing their rates.


            }
        }

    }

}
#endif

void InterfaceKinetics::getDeltaGibbs(doublereal* deltaG)
{
    /*
     * Get the chemical potentials of the species in the
     * ideal gas solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getChemPotentials(DATA_PTR(m_grt) + m_start[n]);
    }
    //for (n = 0; n < m_grt.size(); n++) {
    //    cout << n << "G_RT = " << m_grt[n] << endl;
    //}

    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaG);
}

void InterfaceKinetics::getDeltaElectrochemPotentials(doublereal* deltaM)
{
    /*
     * Get the chemical potentials of the species in the
     * ideal gas solution.
     */
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getElectrochemPotentials(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaM);
}

void InterfaceKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    /*
     * Get the partial molar enthalpy of all species in the
     * ideal gas.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaH);
}

void InterfaceKinetics::getDeltaEntropy(doublereal* deltaS)
{
    /*
     * Get the partial molar entropy of all species in all of
     * the phases
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaS);
}

void InterfaceKinetics::getDeltaSSGibbs(doublereal* deltaG)
{
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaG);
}

void InterfaceKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    /*
     *  Get the standard state enthalpies of the species.
     *  This is the array of chemical potentials at unit activity
     *  We define these here as the enthalpies of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEnthalpy_RT(DATA_PTR(m_grt) + m_start[n]);
    }
    doublereal RT = thermo().temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= RT;
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaH);
}

void InterfaceKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    /*
     *  Get the standard state entropy of the species.
     *  We define these here as the entropies of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEntropy_R(DATA_PTR(m_grt) + m_start[n]);
    }
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= R;
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaS);
}

void InterfaceKinetics::addReaction(ReactionData& r)
{
    /*
     * Install the rate coefficient for the current reaction
     * in the appropriate data structure.
     */
    addElementaryReaction(r);
    /*
     * Add the reactants and products for  m_ropnet;the current reaction
     * to the various stoichiometric coefficient arrays.
     */
    installReagents(r);
    /*
     * Save the reaction and product groups, which are
     * part of the ReactionData class, in this class.
     * They aren't used for anything but reaction path
     * analysis.
     */
    //installGroups(reactionNumber(), r.rgroups, r.pgroups);
    /*
     * Increase the internal number of reactions, m_ii, by one.
     * increase the size of m_perturb by one as well.
     */
    incrementRxnCount();
    m_rxneqn.push_back(r.equation);

    m_rxnPhaseIsReactant.push_back(std::vector<bool>(nPhases(), false));
    m_rxnPhaseIsProduct.push_back(std::vector<bool>(nPhases(), false));

    size_t i = m_ii - 1;
    const std::vector<size_t>& vr = reactants(i);
    for (size_t ik = 0; ik < vr.size(); ik++) {
        size_t k = vr[ik];
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    const std::vector<size_t>& vp = products(i);
    for (size_t ik = 0; ik < vp.size(); ik++) {
        size_t k = vp[ik];
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }
}

void InterfaceKinetics::addElementaryReaction(ReactionData& r)
{
    // install rate coeff calculator
    vector_fp& rp = r.rateCoeffParameters;
    size_t ncov = r.cov.size();
    if (ncov > 3) {
        m_has_coverage_dependence = true;
    }
    for (size_t m = 0; m < ncov; m++) {
        rp.push_back(r.cov[m]);
    }

    /*
     * Temporarily change the reaction rate coefficient type to surface arrhenius.
     * This is what is expected. We'll handle exchange current types below by hand.
     */
    int reactionRateCoeffType_orig = r.rateCoeffType;
    if (r.rateCoeffType == EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE) {
        r.rateCoeffType = SURF_ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }
    if (r.rateCoeffType == ARRHENIUS_REACTION_RATECOEFF_TYPE) {
        r.rateCoeffType = SURF_ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }
    /*
     * Install the reaction rate into the vector of reactions handled by this class
     */
    size_t iloc = m_rates.install(reactionNumber(), r);

    /*
     * Change the reaction rate coefficient type back to its original value
     */
    r.rateCoeffType = reactionRateCoeffType_orig;

    // store activation energy
    m_E.push_back(r.rateCoeffParameters[2]);

    if (r.beta > 0.0) {
        m_has_electrochem_rxns = true;
        m_beta.push_back(r.beta);
        m_ctrxn.push_back(reactionNumber());
        if (r.rateCoeffType == EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE) {
            m_has_exchange_current_density_formulation = true;
            m_ctrxn_ecdf.push_back(1);
        } else {
            m_ctrxn_ecdf.push_back(0);
        }
    }

    // add constant term to rate coeff value vector
    m_rfn.push_back(r.rateCoeffParameters[0]);
    registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);
}

void InterfaceKinetics::setIOFlag(int ioFlag)
{
    m_ioFlag = ioFlag;
    if (m_integrator) {
        m_integrator->setIOFlag(ioFlag);
    }
}

//     void InterfaceKinetics::
//     addGlobalReaction(const ReactionData& r) {

//         int iloc;
//         // install rate coeff calculator
//         vector_fp rp = r.rateCoeffParameters;
//         int ncov = r.cov.size();
//         for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);
//         iloc = m_rates.install( reactionNumber(),
//             r.rateCoeffType, rp.size(),
//             rp.begin() );
//         // store activation energy
//         m_E.push_back(r.rateCoeffParameters[2]);
//         // add constant term to rate coeff value vector
//         m_rfn.push_back(r.rateCoeffParameters[0]);

//         int nr = r.order.size();
//         vector_fp ordr(nr);
//         for (int n = 0; n < nr; n++) {
//             ordr[n] = r.order[n] - r.rstoich[n];
//         }
//         m_globalReactantStoich.add( reactionNumber(),
//             r.reactants, ordr);

//         registerReaction( reactionNumber(), GLOBAL_RXN, iloc);
//     }

void InterfaceKinetics::installReagents(const ReactionData& r)
{

    size_t n, ns, m;
    doublereal nsFlt;
    /*
     * extend temporary storage by one for this rxn.
     */
    m_ropf.push_back(0.0);
    m_ropr.push_back(0.0);
    m_ropnet.push_back(0.0);
    m_rkcn.push_back(0.0);

    /*
     * Obtain the current reaction index for the reaction that we
     * are adding. The first reaction is labeled 0.
     */
    size_t rnum = reactionNumber();

    // vectors rk and pk are lists of species numbers, with
    // repeated entries for species with stoichiometric
    // coefficients > 1. This allows the reaction to be defined
    // with unity reaction order for each reactant, and so the
    // faster method 'multiply' can be used to compute the rate of
    // progress instead of 'power'.

    std::vector<size_t> rk;
    size_t nr = r.reactants.size();
    for (n = 0; n < nr; n++) {
        nsFlt = r.rstoich[n];
        ns = (size_t) nsFlt;
        if ((doublereal) ns != nsFlt) {
            if (ns < 1) {
                ns = 1;
            }
        }
        /*
         * Add to m_rrxn. m_rrxn is a vector of maps. m_rrxn has a length
         * equal to the total number of species for each species, there
         * exists a map, with the reaction number being the key, and the
         * reactant stoichiometric coefficient being the value.
         */
        m_rrxn[r.reactants[n]][rnum] = nsFlt;
        for (m = 0; m < ns; m++) {
            rk.push_back(r.reactants[n]);
        }
    }
    /*
     * Now that we have rk[], we add it into the vector<vector_int> m_reactants
     * in the rnum index spot. Thus m_reactants[rnum] yields a vector
     * of reactants for the rnum'th reaction
     */
    m_reactants.push_back(rk);
    std::vector<size_t> pk;
    size_t np = r.products.size();
    for (n = 0; n < np; n++) {
        nsFlt = r.pstoich[n];
        ns = (size_t) nsFlt;
        if ((doublereal) ns != nsFlt) {
            if (ns < 1) {
                ns = 1;
            }
        }
        /*
         * Add to m_prxn. m_prxn is a vector of maps. m_prxn has a length
         * equal to the total number of species for each species, there
         * exists a map, with the reaction number being the key, and the
         * product stoichiometric coefficient being the value.
         */
        m_prxn[r.products[n]][rnum] = nsFlt;
        for (m = 0; m < ns; m++) {
            pk.push_back(r.products[n]);
        }
    }
    /*
     * Now that we have pk[], we add it into the vector<vector_int> m_products
     * in the rnum index spot. Thus m_products[rnum] yields a vector
     * of products for the rnum'th reaction
     */
    m_products.push_back(pk);
    /*
     * Add this reaction to the stoichiometric coefficient manager. This
     * calculates rates of species production from reaction rates of
     * progress.
     */
    m_rxnstoich.add(reactionNumber(), r);
    /*
     * register reaction in lists of reversible and irreversible rxns.
     */
    if (r.reversible) {
        m_revindex.push_back(reactionNumber());
        m_nrev++;
    } else {
        m_irrev.push_back(reactionNumber());
        m_nirrev++;
    }
}

void InterfaceKinetics::addPhase(thermo_t& thermo)
{
    Kinetics::addPhase(thermo);
    m_phaseExists.push_back(true);
    m_phaseIsStable.push_back(true);
}

void InterfaceKinetics::init()
{
    m_kk = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        m_kk += thermo(n).nSpecies();
    }
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_conc.resize(m_kk);
    m_mu0.resize(m_kk);
    m_grt.resize(m_kk);
    m_pot.resize(m_kk, 0.0);
    m_phi.resize(nPhases(), 0.0);
}

void InterfaceKinetics::finalize()
{
    Kinetics::finalize();
    size_t safe_reaction_size = std::max<size_t>(nReactions(), 1);
    m_rwork.resize(safe_reaction_size);
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("InterfaceKinetics::finalize",
                                           "no surface phase is present.");
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 2)
        throw CanteraError("InterfaceKinetics::finalize",
                           "expected interface dimension = 2, but got dimension = "
                           +int2str(m_surf->nDim()));

    m_StandardConc.resize(m_kk, 0.0);
    m_deltaG0.resize(safe_reaction_size, 0.0);
    m_ProdStanConcReac.resize(safe_reaction_size, 0.0);

    if (m_thermo.size() != m_phaseExists.size()) {
        throw CanteraError("InterfaceKinetics::finalize", "internal error");
    }

    // Guarantee that these arrays can be converted to double* even in the
    // special case where there are no reactions defined.
    if (!m_ii) {
        m_perturb.resize(1, 1.0);
        m_ropf.resize(1, 0.0);
        m_ropr.resize(1, 0.0);
        m_ropnet.resize(1, 0.0);
        m_rkcn.resize(1, 0.0);
    }

    m_finalized = true;
}

doublereal InterfaceKinetics::electrochem_beta(size_t irxn) const
{
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        if (m_ctrxn[i] == irxn) {
            return m_beta[i];
        }
    }
    return 0.0;
}

bool InterfaceKinetics::ready() const
{
    return m_finalized;
}

void InterfaceKinetics::
advanceCoverages(doublereal tstep)
{
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k;
        k.push_back(this);
        m_integrator = new ImplicitSurfChem(k);
        m_integrator->initialize();
    }
    m_integrator->integrate(0.0, tstep);
    delete m_integrator;
    m_integrator = 0;
}

void InterfaceKinetics::
solvePseudoSteadyStateProblem(int ifuncOverride, doublereal timeScaleOverride)
{
    // create our own solver object
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k;
        k.push_back(this);
        m_integrator = new ImplicitSurfChem(k);
        m_integrator->initialize();
    }
    m_integrator->setIOFlag(m_ioFlag);
    /*
     * New direct method to go here
     */
    m_integrator->solvePseudoSteadyStateProblem(ifuncOverride, timeScaleOverride);
}

void InterfaceKinetics::setPhaseExistence(const size_t iphase, const int exists)
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:setPhaseExistence", "out of bounds");
    }
    if (exists) {
        if (!m_phaseExists[iphase]) {
            m_phaseExistsCheck--;
            if (m_phaseExistsCheck < 0) {
                m_phaseExistsCheck = 0;
            }
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
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:phaseExistence()", "out of bounds");
    }
    return m_phaseExists[iphase];
}

int InterfaceKinetics::phaseStability(const size_t iphase) const
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:phaseStability()", "out of bounds");
    }
    return m_phaseIsStable[iphase];
}

void InterfaceKinetics::setPhaseStability(const size_t iphase, const int isStable)
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:setPhaseStability", "out of bounds");
    }
    if (isStable) {
        m_phaseIsStable[iphase] = true;
    } else {
        m_phaseIsStable[iphase] = false;
    }
}

void EdgeKinetics::finalize()
{
    m_rwork.resize(std::max<size_t>(nReactions(), 1));
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("EdgeKinetics::finalize",
                                           "no edge phase is present.");
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 1)
        throw CanteraError("EdgeKinetics::finalize",
                           "expected interface dimension = 1, but got dimension = "
                           +int2str(m_surf->nDim()));

    // Guarantee that these arrays can be converted to double* even in the
    // special case where there are no reactions defined.
    if (!m_ii) {
        m_perturb.resize(1, 1.0);
        m_ropf.resize(1, 0.0);
        m_ropr.resize(1, 0.0);
        m_ropnet.resize(1, 0.0);
        m_rkcn.resize(1, 0.0);
    }

    m_finalized = true;
}

}
