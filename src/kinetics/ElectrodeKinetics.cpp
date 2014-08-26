/**
 *  @file ElectrodeKinetics.cpp
 */

#include "cantera/kinetics/ElectrodeKinetics.h"
#include "cantera/thermo/SurfPhase.h"

using namespace std;

namespace Cantera
{
//============================================================================================================================
ElectrodeKinetics::ElectrodeKinetics(thermo_t* thermo) :
    InterfaceKinetics(thermo),
    metalPhaseRS_(npos),
    electronPhaseRS_(npos),
    solnPhaseRS_(npos),
    kElectronRS_(npos)
{
 
}
//============================================================================================================================
ElectrodeKinetics::~ElectrodeKinetics()
{
  
}
//============================================================================================================================
ElectrodeKinetics::ElectrodeKinetics(const ElectrodeKinetics& right) :
    InterfaceKinetics()

{
    /*
     * Call the assignment operator
     */
    ElectrodeKinetics::operator=(right);
}
//============================================================================================================================
ElectrodeKinetics& ElectrodeKinetics::operator=(const ElectrodeKinetics& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    InterfaceKinetics::operator=(right);

    metalPhaseRS_   = right.metalPhaseRS_;
    electronPhaseRS_ = right.electronPhaseRS_;
    solnPhaseRS_ = right.solnPhaseRS_;
    kElectronRS_ = right.kElectronRS_;
   
    return *this;
}
//============================================================================================================================
int ElectrodeKinetics::type() const
{
    return cInterfaceKinetics;
}
//============================================================================================================================
Kinetics* ElectrodeKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    ElectrodeKinetics* iK = new ElectrodeKinetics(*this);
    iK->assignShallowPointers(tpVector);
    return iK;
}
//============================================================================================================================
//====================================================================================================================
//  Identify the metal phase and the electrons species
void ElectrodeKinetics::identifyMetalPhase()
{
    metalPhaseRS_ = npos;
    kElectronRS_ = -1;
    size_t np = nPhases();
    //
    // Identify the metal phase as the phase with the electron species (element index of 1 for element E
    // Should probably also stipulate a charge of -1.
    //
    for (size_t iph = 0; iph < np; iph++) {
        ThermoPhase* tp = m_thermo[iph];
        size_t nSpecies = tp->nSpecies();
        size_t nElements = tp->nElements();
        size_t eElectron = tp->elementIndex("E");
        if (eElectron != npos) {
            for (size_t k = 0; k < nSpecies; k++) {
                if (tp->nAtoms(k,eElectron) == 1) {
                    int ifound = 1;
                    for (size_t e = 0; e < nElements; e++) {
                        if (tp->nAtoms(k,e) != 0.0) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        metalPhaseRS_ = iph;
                        kElectronRS_ = m_start[iph] + k;
                    }
                }
            }
        }
       //
       //  Identify the solution phase as a 3D phase, with nonzero phase charge change 
       //  in at least one reaction
       //
       if (iph != metalPhaseRS_) {
            for (size_t i = 0; i < m_ii; i++) {
                RxnMolChange* rmc = rmcVector[i];
                if (rmc->m_phaseChargeChange[iph] != 0) {
                    if (rmc->m_phaseDims[iph] == 3) {
                       solnPhaseRS_ = iph;
                       break;
                    }
                }
             }
        }

    }
}
//============================================================================================================================
// virtual from InterfaceKinetics
void ElectrodeKinetics::updateROP()
{
    // evaluate rate constants and equilibrium constants at temperature and phi (electric potential)
    _update_rates_T();
    // get updated activities (rates updated below)
    _update_rates_C();
    
    double TT = m_surf->temperature();
    double rtdf = GasConstant * TT / Faraday;

    if (m_ROP_ok) {
        return;
    }
    //
    // Copy the reaction rate coefficients, m_rfn, into m_ropf
    //
    copy(m_rfn.begin(), m_rfn.end(), m_ropf.begin());
    //
    // Multiply by the perturbation factor
    //
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());
    //
    // Copy the forward rate constants to the reverse rate constants
    //
    copy(m_ropf.begin(), m_ropf.end(), m_ropr.begin());



    //
    // For reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    //
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());
    //
    // multiply ropf by the actyivity concentration reaction orders to obtain
    // the forward rates of progress. 
    //
    m_rxnstoich.multiplyReactants(DATA_PTR(m_actConc), DATA_PTR(m_ropf));
    //
    // For reversible reactions, multiply ropr by the activity concentration products
    //
    m_rxnstoich.multiplyRevProducts(DATA_PTR(m_actConc), DATA_PTR(m_ropr));
    //
    //  Fix up these calculations for cases where the above formalism doesn't hold
    //
    double OCV = 0.0;
    for (size_t iBeta = 0; iBeta < m_beta.size(); iBeta++) {
        size_t irxn = m_ctrxn[iBeta];

	int reactionType = reactionTypes_[irxn];
	if (reactionType == BUTLERVOLMER_RXN) {
	    //
	    //   Get the beta value
	    //
	    double beta = m_beta[iBeta];
	    //
	    // OK, the reaction rate constant contains the current density rate constant calculation
	    // the rxnstoich calculation contained the dependence of the current density on the activity concentrations
	    // We finish up with the ROP calculation
	    //

	    //
	    //   Get the phase mole change structure
	    //
	    RxnMolChange* rmc = rmcVector[irxn];
	    //
	    //   Calculate the stoichiometric eletrons for the reaction
	    //   This is the number of electrons that are the net products of the reaction
	    //
	    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseRS_];
	    //
	    //   Calculate the open circuit voltage of the reaction
	    //
	    getDeltaGibbs(0);
	    if (nStoichElectrons != 0.0) {
		OCV = m_deltaG[irxn]/Faraday/ nStoichElectrons;
	    } else {
                OCV = 0.0;
	    }
	    //
	    //   Calculate the voltage of the electrode.
	    //
	    double voltage = m_phi[metalPhaseRS_] - m_phi[solnPhaseRS_];
	    //
	    //   Calculate the overpotential
	    //
	    double nu = voltage - OCV;
	    //
	    //   Calculate the exchange current density
	    //   m_ropf contains the exchange current reaction rate
	    //
	    double io = m_ropf[irxn] * nStoichElectrons;
	    
	    double exp1 = nu * nStoichElectrons * beta / rtdf;
	    double exp2 = - nu * nStoichElectrons * (1.0 - beta) / (rtdf);
	    m_ropnet[irxn] = io * (exp(exp1) - exp(exp2));

	    // Need to resurrect the forwards rate constant. 
	    //m_ropf[irxn] = ;
	    m_ropr[irxn] =  m_ropnet[irxn] - m_ropf[irxn];
	    
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
//==================================================================================================================
//
//  When the BV form is used we still need to go backwards to calculate the forward rate of progress.
//  This routine does that
//
double ElectrodeKinetics::calcForwardROP_BV(size_t irxn, size_t iBeta)
{
    doublereal rt = GasConstant * thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    //
    //    Calculate gather the exchange current reaction rate constant (where does n_s appear?)
    //
    double iorc = m_rfn[irxn] * m_perturb[irxn];
    //
    //    Determine whether the reaction rate constant is in an exchange current density formulation format.
    //
    int iECDFormulation = m_ctrxn_ecdf[iBeta];

    if (!iECDFormulation) {
	throw CanteraError("", "not handled yet");
    }
    //
    //  Calculate the forward chemical and modify the forward reaction rate coefficient 
    //
    double tmp = exp(- m_beta[iBeta] * m_deltaG0[irxn] * rrt);
    double tmp2 = m_ProdStanConcReac[irxn];
    tmp *= 1.0  / tmp2 / Faraday;
    //
    //  Calculate the chemical reaction rate constant
    //
    double kf = iorc * tmp;
    //
    //  Calculate the electrochemical factor
    // 
    double eamod = m_beta[iBeta] * deltaElectricEnergy_[irxn];
    kf *= exp(- eamod * rrt);
    //
    // Calculate the forward rate of progress
    //   -> get the pointer for the orders
    //
    const RxnOrders* ro_fwd = m_ctrxn_FwdOrdersList_[iBeta];
    if (ro_fwd == 0) {
	throw CanteraError("ElectrodeKinetics::calcForwardROP_BV()", "forward orders pointer is zero ?!?");
    }
    tmp = 1.0;
    const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
    const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
    for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	size_t k = kinSpeciesIDs[j];
	double oo = kinSpeciesOrders[j];
	tmp *= pow(m_actConc[k], oo);
    }
    double ropf = kf * tmp;

    return ropf;
}
//==================================================================================================================



}
