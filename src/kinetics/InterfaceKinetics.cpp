/**
 *  @file InterfaceKinetics.cpp
 */

// Copyright 2002  California Institute of Technology

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/thermo/SurfPhase.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

InterfaceKinetics::InterfaceKinetics(thermo_t* thermo) :
    m_redo_rates(false),
    m_nirrev(0),
    m_nrev(0),
    m_surf(0),
    m_integrator(0),
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
    m_ioFlag(0)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
}

InterfaceKinetics::~InterfaceKinetics()
{
    delete m_integrator;

    for (size_t i = 0; i <  m_ctrxn_ROPOrdersList_.size(); i++) {
        delete m_ctrxn_ROPOrdersList_[i];
    }
    for (size_t i = 0; i <  m_ctrxn_FwdOrdersList_.size(); i++) {
        delete m_ctrxn_FwdOrdersList_[i];
    }
}

InterfaceKinetics::InterfaceKinetics(const InterfaceKinetics& right)
{
    /*
     * Call the assignment operator
     */
    operator=(right);
}

InterfaceKinetics& InterfaceKinetics::operator=(const InterfaceKinetics& right)
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
    m_irrev                = right.m_irrev;
    m_nirrev               = right.m_nirrev;
    m_nrev                 = right.m_nrev;
    m_conc                 = right.m_conc;
    m_actConc              = right.m_actConc;
    m_mu0                  = right.m_mu0;
    m_mu                   = right.m_mu;
    m_mu0_Kc               = right.m_mu0_Kc;
    m_phi                  = right.m_phi;
    m_pot                  = right.m_pot;
    deltaElectricEnergy_   = right.deltaElectricEnergy_;
    m_E                    = right.m_E;
    m_surf                 = right.m_surf;  //DANGER - shallow copy
    m_integrator           = right.m_integrator;  //DANGER - shallow copy
    m_beta                 = right.m_beta;
    m_ctrxn                = right.m_ctrxn;
    m_ctrxn_BVform         = right.m_ctrxn_BVform;
    m_ctrxn_ecdf           = right.m_ctrxn_ecdf;
    m_StandardConc         = right.m_StandardConc;
    m_deltaG0              = right.m_deltaG0;
    m_deltaG               = right.m_deltaG;
    m_ProdStanConcReac     = right.m_ProdStanConcReac;
    m_logp0                = right.m_logp0;
    m_logc0                = right.m_logc0;
    m_ROP_ok               = right.m_ROP_ok;
    m_temp                 = right.m_temp;
    m_logtemp              = right.m_logtemp;
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

    for (size_t i = 0; i <  m_ctrxn_ROPOrdersList_.size(); i++) {
        delete m_ctrxn_ROPOrdersList_[i];
    }
    m_ctrxn_ROPOrdersList_  = right.m_ctrxn_ROPOrdersList_;
    for (size_t i = 0; i <  m_ctrxn_ROPOrdersList_.size(); i++) {
        RxnOrders* ro = right.m_ctrxn_ROPOrdersList_[i];
        m_ctrxn_ROPOrdersList_[i] = new RxnOrders(*ro);
    }

    for (size_t i = 0; i <  m_ctrxn_FwdOrdersList_.size(); i++) {
        delete m_ctrxn_FwdOrdersList_[i];
    }
    m_ctrxn_FwdOrdersList_  = right.m_ctrxn_FwdOrdersList_;
    for (size_t i = 0; i <  m_ctrxn_FwdOrdersList_.size(); i++) {
        RxnOrders* ro = right.m_ctrxn_FwdOrdersList_[i];
        m_ctrxn_FwdOrdersList_[i] = new RxnOrders(*ro);
    }

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
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();
    if (m_has_coverage_dependence) {
        m_surf->getCoverages(DATA_PTR(m_actConc));
        m_rates.update_C(DATA_PTR(m_actConc));
        m_redo_rates = true;
    }

    // Go find the temperature from the surface
    doublereal T = thermo(surfacePhaseIndex()).temperature();
    m_redo_rates = true;
    if (T != m_temp || m_redo_rates) {
        m_logtemp = log(T);

        //  Calculate the forward rate constant by calling m_rates and store it in m_rfn[]
        m_rates.update(T, m_logtemp, DATA_PTR(m_rfn));
        applyStickingCorrection(&m_rfn[0]);

        //  If we need to do conversions between exchange current density formulation and regular formulation
        //  (either way) do it here.
        if (m_has_exchange_current_density_formulation) {
            convertExchangeCurrentDensityFormulation(DATA_PTR(m_rfn));
        }
        if (m_has_electrochem_rxns) {
            applyVoltageKfwdCorrection(DATA_PTR(m_rfn));
        }
        m_temp = T;
        updateKc();
        m_ROP_ok = false;
        m_redo_rates = false;
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

// Updates the internal variables m_actConc and m_conc
void InterfaceKinetics::_update_rates_C()
{
    for (size_t n = 0; n < nPhases(); n++) {
        const ThermoPhase* tp = m_thermo[n];
        /*
         * We call the getActivityConcentrations function of each
         * ThermoPhase class that makes up this kinetics object to
         * obtain the generalized concentrations for species within that
         * class. This is collected in the vector m_conc. m_start[]
         * are integer indices for that vector denoting the start of the
         * species for each phase.
         */
        tp->getActivityConcentrations(DATA_PTR(m_actConc) + m_start[n]);

        // Get regular concentrations too
        tp->getConcentrations(DATA_PTR(m_conc) + m_start[n]);
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

    if (m_nrev > 0) {
        /*
         * Get the vector of standard state electrochemical potentials for species in the Interfacial
         * kinetics object and store it in m_mu0[] and m_mu0_Kc[]
         */
        updateMu0();
        doublereal rrt = 1.0 / (GasConstant * thermo(0).temperature());

        // compute Delta mu^0 for all reversible reactions
        getRevReactionDelta(DATA_PTR(m_mu0_Kc), DATA_PTR(m_rkcn));

        for (size_t i = 0; i < m_nrev; i++) {
            size_t irxn = m_revindex[i];
            if (irxn == npos || irxn >= nReactions()) {
                throw CanteraError("InterfaceKinetics", "illegal value: irxn = "+int2str(irxn));
            }
            // WARNING this may overflow HKM
            m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
        }
        for (size_t i = 0; i != m_nirrev; ++i) {
            m_rkcn[ m_irrev[i] ] = 0.0;
        }
    }
}

void InterfaceKinetics::updateMu0()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    updateExchangeCurrentQuantities();
    /*
     * Get the vector of standard state electrochemical potentials for species in the Interfacial
     * kinetics object and store it in m_mu0[] and in m_mu0_Kc[]
     */
    size_t nsp, ik = 0;
    doublereal rt = GasConstant * thermo(0).temperature();
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
        nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_mu0_Kc[ik] = m_mu0[ik] + Faraday * m_phi[n] * thermo(n).charge(k);
            m_mu0_Kc[ik] -= rt * thermo(n).logStandardConc(k);
            ik++;
        }
    }
}

void InterfaceKinetics::checkPartialEquil()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    vector_fp dmu(nTotalSpecies(), 0.0);
    vector_fp rmu(std::max<size_t>(nReactions(), 1), 0.0);
    if (m_nrev > 0) {
        doublereal rt = GasConstant*thermo(0).temperature();
        cout << "T = " << thermo(0).temperature() << " " << rt << endl;
        size_t nsp, ik=0;
        doublereal delta;
        for (size_t n = 0; n < nPhases(); n++) {
            thermo(n).getChemPotentials(DATA_PTR(dmu) + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                delta = Faraday * m_phi[n] * thermo(n).charge(k);
                dmu[ik] += delta;
                ik++;
            }
        }

        // compute Delta mu^ for all reversible reactions
        getRevReactionDelta(DATA_PTR(dmu), DATA_PTR(rmu));
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

void InterfaceKinetics::getEquilibriumConstants(doublereal* kc)
{
    updateMu0();
    doublereal rrt = 1.0 / (GasConstant * thermo(0).temperature());

    std::fill(kc, kc + m_ii, 0.0);

    getReactionDelta(DATA_PTR(m_mu0_Kc), kc);

    for (size_t i = 0; i < m_ii; i++) {
        kc[i] = exp(-kc[i]*rrt);
    }
}

void InterfaceKinetics::updateExchangeCurrentQuantities()
{
    /*
     *   Calculate:
     *     - m_StandardConc[]
     *     - m_ProdStandConcReac[]
     *     - m_deltaG0[]
     *     - m_mu0[]
     */

    /*
     * First collect vectors of the standard Gibbs free energies of the
     * species and the standard concentrations
     *   - m_mu0
     *   - m_StandardConc
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

    getReactionDelta(DATA_PTR(m_mu0), DATA_PTR(m_deltaG0));

    //  Calculate the product of the standard concentrations of the reactants
    for (size_t i = 0; i < m_ii; i++) {
        m_ProdStanConcReac[i] = 1.0;
    }
    m_reactantStoich.multiply(DATA_PTR(m_StandardConc), DATA_PTR(m_ProdStanConcReac));
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

    // Compute the change in electrical potential energy for each
    // reaction. This will only be non-zero if a potential
    // difference is present.
    getReactionDelta(DATA_PTR(m_pot), DATA_PTR(deltaElectricEnergy_));

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

        //   If we calculate the BV form directly, we don't add the voltage correction to the
        //   forward reaction rate constants.
        if (m_ctrxn_BVform[i] == 0) {
            eamod = m_beta[i] * deltaElectricEnergy_[irxn];
            if (eamod != 0.0) {
#ifdef DEBUG_KIN_MODE
                ea = GasConstant * m_E[irxn];
                if (eamod + ea < 0.0) {
                    writelog("Warning: act energy mod too large!\n");
                    writelog("  Delta phi = "+fp2str(deltaElectricEnergy_[irxn]/Faraday)+"\n");
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
}

void InterfaceKinetics::convertExchangeCurrentDensityFormulation(doublereal* const kfwd)
{
    updateExchangeCurrentQuantities();
    doublereal rt = GasConstant * thermo(0).temperature();
    doublereal rrt = 1.0/rt;

    //  Loop over all reactions which are defined to have a voltage transfer coefficient that
    //  affects the activity energy for the reaction
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        size_t irxn = m_ctrxn[i];

        // Determine whether the reaction rate constant is in an exchange current density formulation format.
        int iECDFormulation =  m_ctrxn_ecdf[i];
        if (iECDFormulation) {
            // If the BV form is to be converted into the normal form then we go through this process.
            // If it isn't to be converted, then we don't go through this process.
            //
            //   We need to have the straight chemical reaction rate constant to come out of this calculation.
            if (m_ctrxn_BVform[i] == 0) {
                //
                //  Calculate the term and modify the forward reaction
                //
                double tmp = exp(- m_beta[i] * m_deltaG0[irxn] * rrt);
                double tmp2 = m_ProdStanConcReac[irxn];
                tmp *= 1.0  / tmp2 / Faraday;
                kfwd[irxn] *= tmp;
            }
            //  If BVform is nonzero we don't need to do anything.

        } else {
            // kfwd[] is the chemical reaction rate constant
            //
            // If we are to calculate the BV form directly, then we will do the reverse.
            // We will calculate the exchange current density formulation here and
            // substitute it.
            if (m_ctrxn_BVform[i] != 0) {

                //  Calculate the term and modify the forward reaction rate constant so that
                //  it's in the exchange current density formulation format
                double tmp = exp(m_beta[i] * m_deltaG0[irxn] * rrt);
                double tmp2 = m_ProdStanConcReac[irxn];
                tmp *= Faraday * tmp2;
                kfwd[irxn] *= tmp;
            }
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

void InterfaceKinetics::updateROP()
{
    // evaluate rate constants and equilibrium constants at temperature and phi (electric potential)
    _update_rates_T();
    // get updated activities (rates updated below)
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    // Copy the reaction rate coefficients, m_rfn, into m_ropf
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());

    // Multiply by the perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());
    //
    // Copy the forward rate constants to the reverse rate constants
    //
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());

    // For reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by the activity concentration reaction orders to obtain
    // the forward rates of progress.
    m_reactantStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(m_ropf));

    // For reversible reactions, multiply ropr by the activity concentration products
    m_revProductStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(m_ropr));

    //  Fix up these calculations for cases where the above formalism doesn't hold
    double OCV = 0.0;
    for (size_t jrxn = 0; jrxn != m_ii; ++jrxn) {
        int reactionType = m_rxntype[jrxn];
        if (reactionType == BUTLERVOLMER_RXN) {
            //
            // OK, the reaction rate constant contains the current density rate constant calculation
            // the rxnstoich calculation contained the dependence of the current density on the activity concentrations
            // We finish up with the ROP calculation
            //
            // Calculate the overpotential of the reaction
            //
            //    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseRS_];
            double nStoichElectrons=1;
            //*nStoich = nStoichElectrons;



            getDeltaGibbs(0);

            if (nStoichElectrons != 0.0) {
                OCV = m_deltaG[jrxn]/Faraday/ nStoichElectrons;
            }

            /*

            double exp1 = nu * nStoich * beta / rtdf
            double exp2 = -nu * nStoich * Faraday * (1.0 - beta) / (GasConstant * temp);
            double val = io * (exp(exp1) - exp(exp2));

            doublereal BVterm = exp(exp1  ) - exp(exp2);
            m_ropnet[j] = m_ropf[j] * BVterm
            m_ropf[j] =
                //
            m_ropr[j] =  m_ropnet[j] - m_ropf[j];
            */
        }
    }

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

void InterfaceKinetics::getDeltaGibbs(doublereal* deltaG)
{
    /*
     * Get the chemical potentials of the species in the all of the phases used in the
     * kinetics mechanism
     */
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(DATA_PTR(m_mu) + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(DATA_PTR(m_mu), DATA_PTR(m_deltaG));
    if (deltaG != 0 && (DATA_PTR(m_deltaG) != deltaG)) {
        for (size_t j = 0; j < m_ii; ++j) {
            deltaG[j] = m_deltaG[j];
        }
    }
}

void InterfaceKinetics::getDeltaElectrochemPotentials(doublereal* deltaM)
{
    /*
     * Get the chemical potentials of the species
     */
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getElectrochemPotentials(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    getReactionDelta(DATA_PTR(m_grt), deltaM);
}

void InterfaceKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    /*
     * Get the partial molar enthalpy of all species
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    getReactionDelta(DATA_PTR(m_grt), deltaH);
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
    getReactionDelta(DATA_PTR(m_grt), deltaS);
}

void InterfaceKinetics::getDeltaSSGibbs(doublereal* deltaGSS)
{
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    getReactionDelta(DATA_PTR(m_mu0), deltaGSS);
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
    doublereal RT = thermo(0).temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= RT;
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    getReactionDelta(DATA_PTR(m_grt), deltaH);
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
    getReactionDelta(DATA_PTR(m_grt), deltaS);
}

void InterfaceKinetics::addReaction(ReactionData& r)
{
    int reactionType = r.reactionType;

    // Install rate coeff calculator
    if (r.cov.size() > 3) {
        m_has_coverage_dependence = true;
    }
    for (size_t m = 0; m < r.cov.size(); m++) {
        r.rateCoeffParameters.push_back(r.cov[m]);
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
    m_rates.install(m_ii, r);

    /*
     * Change the reaction rate coefficient type back to its original value
     */
    r.rateCoeffType = reactionRateCoeffType_orig;

    // Store activation energy
    m_E.push_back(r.rateCoeffParameters[2]);

    if (r.beta > 0.0) {
        m_has_electrochem_rxns = true;
        m_beta.push_back(r.beta);
        m_ctrxn.push_back(m_ii);
        if (r.rateCoeffType == EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE) {
            m_has_exchange_current_density_formulation = true;
            m_ctrxn_ecdf.push_back(1);
        } else {
            m_ctrxn_ecdf.push_back(0);
        }
        m_ctrxn_resistivity_.push_back(r.filmResistivity);

        if (reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
            reactionType == BUTLERVOLMER_RXN ||
            reactionType == SURFACEAFFINITY_RXN ||
            reactionType == GLOBAL_RXN) {
            //   Specify alternative forms of the electrochemical reaction
            if (r.reactionType == BUTLERVOLMER_RXN) {
                m_ctrxn_BVform.push_back(1);
            } else if (r.reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {
                m_ctrxn_BVform.push_back(2);
            } else {
                // set the default to be the normal forward / reverse calculation method
                m_ctrxn_BVform.push_back(0);
            }
            if (r.forwardFullOrder_.size() > 0) {
                RxnOrders* ro = new RxnOrders();
                ro->fill(r.forwardFullOrder_);
                m_ctrxn_ROPOrdersList_.push_back(ro);
                m_ctrxn_FwdOrdersList_.push_back(0);

                // Fill in the Fwd Orders dependence here for B-V reactions
                if (r.reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
                    r.reactionType == BUTLERVOLMER_RXN) {
                    vector_fp fwdFullorders(m_kk, 0.0);
                    determineFwdOrdersBV(r, fwdFullorders);
                    RxnOrders* ro = new RxnOrders();
                    ro->fill(fwdFullorders);
                    m_ctrxn_FwdOrdersList_[m_ii] = ro;
                }
            } else {
                m_ctrxn_ROPOrdersList_.push_back(0);
                m_ctrxn_FwdOrdersList_.push_back(0);
            }

        } else {
            m_ctrxn_BVform.push_back(0);
            m_ctrxn_ROPOrdersList_.push_back(0);
            m_ctrxn_FwdOrdersList_.push_back(0);
            if (r.filmResistivity > 0.0) {
                throw CanteraError("InterfaceKinetics::addReaction()",
                                   "film resistivity set for elementary reaction");
            }
        }
    }

    if (r.reversible) {
        m_revindex.push_back(nReactions());
        m_nrev++;
    } else {
        m_irrev.push_back(nReactions());
        m_nirrev++;
    }
    Kinetics::addReaction(r);

    m_rxnPhaseIsReactant.push_back(std::vector<bool>(nPhases(), false));
    m_rxnPhaseIsProduct.push_back(std::vector<bool>(nPhases(), false));

    size_t i = m_ii - 1;
    for (size_t ik = 0; ik < r.reactants.size(); ik++) {
        size_t k = r.reactants[ik];
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    for (size_t ik = 0; ik < r.products.size(); ik++) {
        size_t k = r.products[ik];
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }
}

bool InterfaceKinetics::addReaction(shared_ptr<Reaction> r_base)
{
    size_t i = nReactions();
    bool added = Kinetics::addReaction(r_base);
    if (!added) {
        return false;
    }

    InterfaceReaction& r = dynamic_cast<InterfaceReaction&>(*r_base);
    SurfaceArrhenius rate = buildSurfaceArrhenius(i, r);
    m_rates.install(i, rate);

    // Turn on the global flag indicating surface coverage dependence
    if (!r.coverage_deps.empty()) {
        m_has_coverage_dependence = true;
    }

    // Store activation energy
    m_E.push_back(rate.activationEnergy_R());

    ElectrochemicalReaction* re = dynamic_cast<ElectrochemicalReaction*>(&r);
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
        m_ctrxn_resistivity_.push_back(re->film_resistivity);

        if (r.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
            r.reaction_type == BUTLERVOLMER_RXN ||
            r.reaction_type == SURFACEAFFINITY_RXN ||
            r.reaction_type == GLOBAL_RXN) {
            //   Specify alternative forms of the electrochemical reaction
            if (r.reaction_type == BUTLERVOLMER_RXN) {
                m_ctrxn_BVform.push_back(1);
            } else if (r.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {
                m_ctrxn_BVform.push_back(2);
            } else {
                // set the default to be the normal forward / reverse calculation method
                m_ctrxn_BVform.push_back(0);
            }
            if (!r.orders.empty()) {
                vector_fp orders(nTotalSpecies(), 0.0);
                for (Composition::const_iterator iter = r.orders.begin();
                     iter != r.orders.end();
                     ++iter) {
                    orders[kineticsSpeciesIndex(iter->first)] = iter->second;
                }
                RxnOrders* ro = new RxnOrders();
                ro->fill(orders);
                m_ctrxn_ROPOrdersList_.push_back(ro);
                m_ctrxn_FwdOrdersList_.push_back(0);

                // Fill in the Fwd Orders dependence here for B-V reactions
                if (r.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
                    r.reaction_type == BUTLERVOLMER_RXN) {
                    vector_fp fwdFullorders(m_kk, 0.0);
                    determineFwdOrdersBV(*re, fwdFullorders);
                    RxnOrders* ro = new RxnOrders();
                    ro->fill(fwdFullorders);
                    m_ctrxn_FwdOrdersList_[i] = ro;
                }
            } else {
                m_ctrxn_ROPOrdersList_.push_back(0);
                m_ctrxn_FwdOrdersList_.push_back(0);
            }

        } else {
            m_ctrxn_BVform.push_back(0);
            m_ctrxn_ROPOrdersList_.push_back(0);
            m_ctrxn_FwdOrdersList_.push_back(0);
            if (re->film_resistivity > 0.0) {
                throw CanteraError("InterfaceKinetics::addReaction()",
                                   "film resistivity set for elementary reaction");
            }
        }
    }

    if (r.reversible) {
        m_revindex.push_back(i);
        m_nrev++;
    } else {
        m_irrev.push_back(i);
        m_nirrev++;
    }

    m_rxnPhaseIsReactant.push_back(std::vector<bool>(nPhases(), false));
    m_rxnPhaseIsProduct.push_back(std::vector<bool>(nPhases(), false));

    for (Composition::const_iterator iter = r.reactants.begin();
         iter != r.reactants.end();
         ++iter) {
        size_t k = kineticsSpeciesIndex(iter->first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    for (Composition::const_iterator iter = r.products.begin();
         iter != r.products.end();
         ++iter) {
        size_t k = kineticsSpeciesIndex(iter->first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }
    return true;
}

void InterfaceKinetics::modifyReaction(size_t i, shared_ptr<Reaction> r_base)
{
    Kinetics::modifyReaction(i, r_base);
    InterfaceReaction& r = dynamic_cast<InterfaceReaction&>(*r_base);
    SurfaceArrhenius rate = buildSurfaceArrhenius(npos, r);
    m_rates.replace(i, rate);

    // Invalidate cached data
    m_redo_rates = true;
    m_temp += 0.1;
}

SurfaceArrhenius InterfaceKinetics::buildSurfaceArrhenius(
    size_t i, InterfaceReaction& r)
{
    double A_rate = r.rate.preExponentialFactor();
    double b_rate = r.rate.temperatureExponent();

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

        b_rate += 0.5;
        std::string sticking_species = r.sticking_species;
        if (sticking_species == "") {
            // Identify the sticking species if not explicitly given
            bool foundStick = false;
            for (Composition::const_iterator iter = r.reactants.begin();
                 iter != r.reactants.end();
                 ++iter) {
                size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(iter->first));
                if (iPhase != iInterface) {
                    // Non-interface species. There should be exactly one of these
                    if (foundStick) {
                        throw CanteraError("InterfaceKinetics::addReaction",
                            "Multiple non-interface species found"
                            "in sticking reaction: '" + r.equation() + "'");
                    }
                    foundStick = true;
                    sticking_species = iter->first;
                }
            }
            if (!foundStick) {
                throw CanteraError("InterfaceKinetics::addReaction",
                    "No non-interface species found"
                    "in sticking reaction: '" + r.equation() + "'");
            }
        }

        double surface_order = 0.0;
        // Adjust the A-factor
        for (Composition::const_iterator iter = r.reactants.begin();
             iter != r.reactants.end();
             ++iter) {
            size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(iter->first));
            const ThermoPhase& p = thermo(iPhase);
            const ThermoPhase& surf = thermo(surfacePhaseIndex());
            size_t k = p.speciesIndex(iter->first);
            if (iter->first == sticking_species) {
                A_rate *= sqrt(GasConstant/(2*Pi*p.molecularWeight(k)));
            } else {
                // Non-sticking species. Convert from coverages used in the
                // sticking probability expression to the concentration units
                // used in the mass action rate expression. For surface phases,
                // the dependence on the site density is incorporated when the
                // rate constant is evaluated, since we don't assume that the
                // site density is known at this time.
                double order = getValue(r.orders, iter->first, iter->second);
                if (&p == &surf) {
                    A_rate *= pow(p.size(k), order);
                    surface_order += order;
                } else {
                    A_rate *= pow(p.standardConcentration(k), -order);
                }
            }
        }
        if (i != npos) {
            m_sticking_orders.push_back(make_pair(i, surface_order));
        }
    }

    SurfaceArrhenius rate(A_rate, b_rate, r.rate.activationEnergy_R());

    // Set up coverage dependencies
    for (map<string, CoverageDependency>::const_iterator iter = r.coverage_deps.begin();
         iter != r.coverage_deps.end();
         ++iter) {
        size_t k = thermo(reactionPhaseIndex()).speciesIndex(iter->first);
        rate.addCoverageDependence(k, iter->second.a, iter->second.m, iter->second.E);
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
    m_actConc.resize(m_kk);
    m_conc.resize(m_kk);
    m_mu0.resize(m_kk);
    m_mu.resize(m_kk);
    m_mu0_Kc.resize(m_kk);
    m_grt.resize(m_kk);
    m_pot.resize(m_kk, 0.0);
    m_phi.resize(nPhases(), 0.0);
}

void InterfaceKinetics::finalize()
{
    Kinetics::finalize();
    size_t safe_reaction_size = std::max<size_t>(m_ii, 1);
    deltaElectricEnergy_.resize(safe_reaction_size);
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("InterfaceKinetics::finalize",
                                           "no surface phase is present.");

    // Check to see that the interface routine has a dimension of 2
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 2) {
        throw CanteraError("InterfaceKinetics::finalize",
                           "expected interface dimension = 2, but got dimension = "
                           +int2str(m_surf->nDim()));
    }
    m_StandardConc.resize(m_kk, 0.0);
    m_deltaG0.resize(safe_reaction_size, 0.0);
    m_deltaG.resize(safe_reaction_size, 0.0);

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

void InterfaceKinetics::advanceCoverages(doublereal tstep)
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

void InterfaceKinetics::solvePseudoSteadyStateProblem(
    int ifuncOverride, doublereal timeScaleOverride)
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

void InterfaceKinetics::determineFwdOrdersBV(ReactionData& rdata, std::vector<doublereal>& fwdFullorders)
{
    //   Start out with the full ROP orders vector.
    //   This vector will have the BV exchange current density orders in it.
    fwdFullorders = rdata.forwardFullOrder_;

    //   forward and reverse beta values
    double betaf = rdata.beta;

    //   Loop over the reactants doing away with the BV terms.
    //   This should leave the reactant terms only, even if they are non-mass action.
    for (size_t j = 0; j < rdata.reactants.size(); j++) {
        size_t kkin =  rdata.reactants[j];
        double oo = rdata.rstoich[j];
        fwdFullorders[kkin] += betaf * oo;
        // just to make sure roundoff doesn't leave a term that should be zero (haven't checked this out yet)
        if (abs(fwdFullorders[kkin]) < 0.00001) {
            fwdFullorders[kkin] = 0.0;
        }
    }

    //   Loop over the products doing away with the BV terms.
    //   This should leave the reactant terms only, even if they are non-mass action.
    for (size_t j = 0; j < rdata.products.size(); j++) {
        size_t kkin =  rdata.products[j];
        double oo = rdata.pstoich[j];
        fwdFullorders[kkin] -= betaf * oo;
        if (abs(fwdFullorders[kkin]) < 0.00001) {
            fwdFullorders[kkin] = 0.0;
        }
    }
}

void InterfaceKinetics::determineFwdOrdersBV(ElectrochemicalReaction& r, std::vector<doublereal>& fwdFullOrders)
{
    // Start out with the full ROP orders vector.
    // This vector will have the BV exchange current density orders in it.
    fwdFullOrders.assign(nTotalSpecies(), 0.0);
    for (Composition::const_iterator iter = r.orders.begin();
         iter != r.orders.end();
         ++iter) {
        fwdFullOrders[kineticsSpeciesIndex(iter->first)] = iter->second;
    }

    //   forward and reverse beta values
    double betaf = r.beta;

    // Loop over the reactants doing away with the BV terms.
    // This should leave the reactant terms only, even if they are non-mass action.
    for (Composition::const_iterator iter = r.reactants.begin();
         iter != r.reactants.end();
         ++iter) {
        size_t k = kineticsSpeciesIndex(iter->first);
        fwdFullOrders[k] += betaf * iter->second;
        // just to make sure roundoff doesn't leave a term that should be zero (haven't checked this out yet)
        if (abs(fwdFullOrders[k]) < 0.00001) {
            fwdFullOrders[k] = 0.0;
        }
    }

    // Loop over the products doing away with the BV terms.
    // This should leave the reactant terms only, even if they are non-mass action.
    for (Composition::const_iterator iter = r.products.begin();
         iter != r.products.end();
         ++iter) {
        size_t k = kineticsSpeciesIndex(iter->first);
        fwdFullOrders[k] -= betaf * iter->second;
        // just to make sure roundoff doesn't leave a term that should be zero (haven't checked this out yet)
        if (abs(fwdFullOrders[k]) < 0.00001) {
            fwdFullOrders[k] = 0.0;
        }
    }
}

void InterfaceKinetics::applyStickingCorrection(double* kf)
{
    if (m_sticking_orders.empty()) {
        return;
    }

    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    vector_fp& factors = cached.value;

    SurfPhase& surf = dynamic_cast<SurfPhase&>(thermo(reactionPhaseIndex()));
    double n0 = surf.siteDensity();
    if (!cached.validate(n0)) {
        factors.resize(m_sticking_orders.size());
        for (size_t n = 0; n < m_sticking_orders.size(); n++) {
            factors[n] = pow(n0, -m_sticking_orders[n].second);
        }
    }

    for (size_t n = 0; n < m_sticking_orders.size(); n++) {
        kf[m_sticking_orders[n].first] *= factors[n];
    }
}


void EdgeKinetics::finalize()
{
    //  Note we can't call the Interface::finalize() routine because we need to check for a dimension of 1 below.
    //  Therefore, we have to malloc room in arrays that would normally be
    //  handled by the InterfaceKinetics::finalize() call.
    Kinetics::finalize();

    size_t safe_reaction_size = std::max<size_t>(m_ii, 1);
    deltaElectricEnergy_.resize(safe_reaction_size);
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("EdgeKinetics::finalize",
                                           "no surface phase is present.");

    // Check to see edge phase has a dimension of 1
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 1) {
        throw CanteraError("EdgeKinetics::finalize",
                           "expected interface dimension = 1, but got dimension = "
                           +int2str(m_surf->nDim()));
    }
    m_StandardConc.resize(m_kk, 0.0);
    m_deltaG0.resize(safe_reaction_size, 0.0);
    m_deltaG.resize(safe_reaction_size, 0.0);

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

RxnOrders::RxnOrders(const RxnOrders& right) :
    kinSpeciesIDs_(right.kinSpeciesIDs_),
    kinSpeciesOrders_(right.kinSpeciesOrders_)
{
}

RxnOrders& RxnOrders::operator=(const RxnOrders& right)
{
    if (this == &right) {
        return *this;
    }
    kinSpeciesIDs_    = right.kinSpeciesIDs_;
    kinSpeciesOrders_ = right.kinSpeciesOrders_;
    return *this;
}

int RxnOrders::fill(const std::vector<doublereal>& fullForwardOrders)
{
    int nzeroes = 0;
    kinSpeciesIDs_.clear();
    kinSpeciesOrders_.clear();
    for (size_t k = 0; k < fullForwardOrders.size(); ++k) {
        if (fullForwardOrders[k] != 0.0) {
            kinSpeciesIDs_.push_back(k);
            kinSpeciesOrders_.push_back(fullForwardOrders[k]);
            ++nzeroes;
        }
    }
    return nzeroes;
}

}
