/**
 *  @file ElectrodeKinetics.cpp
 */

#include "cantera/kinetics/ElectrodeKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"

#include <cstdio>

using namespace std;

namespace Cantera
{
//============================================================================================================================
ElectrodeKinetics::ElectrodeKinetics(thermo_t* thermo) :
    InterfaceKinetics(thermo),
    metalPhaseIndex_(npos),
    solnPhaseIndex_(npos),
    kElectronIndex_(npos)
{
    warn_deprecated("class ElectrodeKinetics",
        "To be removed after Cantera 2.2.");
}
//============================================================================================================================
ElectrodeKinetics::~ElectrodeKinetics()
{
    for (size_t i = 0; i < rmcVector.size(); i++) {
        delete rmcVector[i];
    }
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

    metalPhaseIndex_   = right.metalPhaseIndex_;
    solnPhaseIndex_ = right.solnPhaseIndex_;
    kElectronIndex_ = right.kElectronIndex_;

    for (size_t i = 0; i <  rmcVector.size(); i++) {
        delete rmcVector[i];
    }
    rmcVector.resize(m_ii, 0);
    for (size_t i = 0; i < m_ii; i++) {
        if (right.rmcVector[i]) {
            rmcVector[i] = new RxnMolChange(*(right.rmcVector[i]));
        }
    }

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
//  Identify the metal phase and the electron species
void ElectrodeKinetics::identifyMetalPhase()
{
    metalPhaseIndex_ = npos;
    kElectronIndex_ = npos;
    solnPhaseIndex_ = npos;
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
                        metalPhaseIndex_ = iph;
                        kElectronIndex_ = m_start[iph] + k;
                    }
                }
            }
        }
       //
       //  Identify the solution phase as a 3D phase, with nonzero phase charge change 
       //  in at least one reaction
       //
	/*
	 *  Haven't filled in reactions yet when this is called, unlike previous treatment.
       if (iph != metalPhaseIndex_) {
            for (size_t i = 0; i < m_ii; i++) {
                RxnMolChange* rmc = rmcVector[i];
                if (rmc->m_phaseChargeChange[iph] != 0) {
                    if (rmc->m_phaseDims[iph] == 3) {
                       solnPhaseIndex_ = iph;
                       break;
                    }
                }
             }
        }
	*/
	//
	// New method is to find the first multispecies 3D phase with charged species as the solution phase
	//
	if (iph != metalPhaseIndex_) {
	    ThermoPhase& tp =*( m_thermo[iph]);
	    size_t nsp = tp.nSpecies();
	    size_t nd = tp.nDim();
	    if (nd == 3 && nsp > 1) {
		for (size_t k = 0; k < nsp; k++) {
		    if (tp.charge(k) != 0.0) {
			solnPhaseIndex_ = iph;
                        string ss = tp.name();
			// cout << "solution phase = "<< ss << endl;
			break;
		    }
		}
	    }
        }

    }
    //
    //  Right now, if we don't find an electron phase, we will not error exit. Some functions will
    //  be turned off and the object will behave as an InterfaceKinetics object.  This is needed
    //  because downstream electrode objects have internal reaction surfaces that don't have
    //  electrons.
    //
    /*
    if (metalPhaseIndex_ == npos) {
	throw CanteraError("ElectrodeKinetics::identifyMetalPhase()",
			   "Can't find electron phase -> treating this as an error right now");
    }
    if (solnPhaseIndex_ == npos) {
	throw CanteraError("ElectrodeKinetics::identifyMetalPhase()",
			   "Can't find solution phase -> treating this as an error right now");
    }
    */
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
    // multiply ropf by the activity concentration reaction orders to obtain
    // the forward rates of progress. 
    //
    m_reactantStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(m_ropf));
    //
    // For reversible reactions, multiply ropr by the activity concentration products
    //
    m_revProductStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(m_ropr));
    //
    //  Fix up these calculations for cases where the above formalism doesn't hold
    //
    double OCV = 0.0;
    for (size_t iBeta = 0; iBeta < m_beta.size(); iBeta++) {
        size_t irxn = m_ctrxn[iBeta];

	int reactionType = m_rxntype[irxn];
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
	    int iECDFormulation =  m_ctrxn_ecdf[iBeta];
            if (iECDFormulation == 0) {
		throw CanteraError(" ElectrodeKinetics::updateROP()",
				   "Straight kfwrd with  BUTLERVOLMER_RXN not handled yet");
	    }
	    //
	    //   Get the phase mole change structure
	    //
	    RxnMolChange* rmc = rmcVector[irxn];
	    //
	    //   Calculate the stoichiometric eletrons for the reaction
	    //   This is the number of electrons that are the net products of the reaction
	    //
            AssertThrow(metalPhaseIndex_ != npos, "ElectrodeKinetics::updateROP()");

	    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseIndex_];
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
	    double voltage = m_phi[metalPhaseIndex_] - m_phi[solnPhaseIndex_];
	    //
	    //   Calculate the overpotential
	    //
	    double nu = voltage - OCV;

	    //
	    //  Find the product of the standard concentrations for ROP orders that we used above
	    //
	    const RxnOrders* ro_rop = m_ctrxn_ROPOrdersList_[iBeta];
	    if (ro_rop == 0) {
		throw CanteraError("ElectrodeKinetics::", "ROP orders pointer is zero ?!?");
	    }
	    double tmp2 = 1.0;
	    const std::vector<size_t>& kinSpeciesIDs = ro_rop->kinSpeciesIDs_;
	    const std::vector<doublereal>& kinSpeciesOrders = ro_rop->kinSpeciesOrders_;
	    for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
		size_t k = kinSpeciesIDs[j];
		double oo = kinSpeciesOrders[j];
		tmp2 *= pow(m_StandardConc[k], oo);
	    }
	    //
	    //  Now have to divide this to get rid of standard concentrations. We should
	    //  have used just the activities in the m_rxnstoich.multiplyReactants(DATA_PTR(m_actConc), DATA_PTR(m_ropf));
	    //  calculation above!
	    //  That is because the exchange current density rate constants have the correct units in the first place.
	    //
	    m_ropf[irxn] /= tmp2;
	    //
	    //   Calculate the exchange current density
	    //   m_ropf contains the exchange current reaction rate
	    //
	    double ioc = m_ropf[irxn] * nStoichElectrons;
	    //
	    //   Add in the film resistance here
	    // 
	    double resist = m_ctrxn_resistivity_[iBeta];
	    double exp1 = nu * nStoichElectrons * beta / rtdf;
	    double exp2 = - nu * nStoichElectrons * (1.0 - beta) / (rtdf);
	    double io = ioc * (exp(exp1) - exp(exp2));

	    if (resist != 0.0) {
		io = solveCurrentRes(nu, nStoichElectrons, ioc, beta, TT, resist, 0);
	    }

            m_ropnet[irxn] = io / (Faraday * nStoichElectrons);
	    //
	    // Need to resurrect the forwards rate of progress -> there is some need to
	    // calculate each direction individually
	    //
	    m_ropf[irxn] = calcForwardROP_BV(irxn, iBeta, ioc, nStoichElectrons, nu, io);
	    //
	    // Calculate the reverse rate of progress from the difference
	    //
	    m_ropr[irxn] =  m_ropf[irxn] - m_ropnet[irxn];
	    
        } else if (reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {
	    //
	    //   Get the beta value
	    //
	    double beta = m_beta[iBeta];
	    //
	    // OK, the reaction rate constant contains the current density rate constant calculation
	    // the rxnstoich calculation contained the dependence of the current density on the activity concentrations
	    // We finish up with the ROP calculation
	    //
	    int iECDFormulation =  m_ctrxn_ecdf[iBeta];
            if (iECDFormulation == 0) {
		throw CanteraError("ElectrodeKinetics::updateROP()",
				   "Straight kfwrd with BUTLERVOLMER_NOACTIVITYCOEFFS_RXN not handled yet");
	    }
	    //
	    //   Get the phase mole change structure
	    //
	    RxnMolChange* rmc = rmcVector[irxn];
	    //
	    //   Calculate the stoichiometric eletrons for the reaction
	    //   This is the number of electrons that are the net products of the reaction
	    //
	    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseIndex_];
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
	    double voltage = m_phi[metalPhaseIndex_] - m_phi[solnPhaseIndex_];
	    //
	    //   Calculate the overpotential
	    //
	    double nu = voltage - OCV;
	    //
	    //  Unfortunately, we really need to recalculate everything from almost scratch
	    //  for this case, since it widely diverges from the thermo norm.
	    //
	    //   Start with the exchange current reaction rate constant, which should
	    //   be located in m_rfn[].
	    //
	    double ioc = m_rfn[irxn] *  nStoichElectrons *  m_perturb[irxn];
	    //
	    //   Now we need th mole fraction vector and we need the RxnOrders vector.
	    //
	    const RxnOrders* ro_fwd = m_ctrxn_ROPOrdersList_[iBeta];
	    if (ro_fwd == 0) {
		throw CanteraError("ElectrodeKinetics::calcForwardROP_BV()", "forward orders pointer is zero ?!?");
	    }
	    double tmp = 1.0;
	    double mfS = 0.0;
	    const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
	    const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
	    for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
		size_t ks = kinSpeciesIDs[j];
		thermo_t& th = speciesPhase(ks);
		size_t n = speciesPhaseIndex(ks);
                size_t klocal = ks - m_start[n];
		mfS = th.moleFraction(klocal);

		double oo = kinSpeciesOrders[j];
		tmp *= pow(mfS, oo);
	    }
	    ioc *= tmp;
	    //
	    //   Add in the film resistance here, later
	    //
	    double resist = m_ctrxn_resistivity_[iBeta];
	    double exp1 = nu * nStoichElectrons * beta / rtdf;
	    double exp2 = - nu * nStoichElectrons * (1.0 - beta) / (rtdf);
	    double io = ioc * (exp(exp1) - exp(exp2));
	    if (resist != 0.0) {
		io = solveCurrentRes(nu, nStoichElectrons, ioc, beta, TT, resist, 0);
	    }

            m_ropnet[irxn] = io / (Faraday * nStoichElectrons);
	    //
	    // Need to resurrect the forwards rate of progress -> there is some need to
	    // calculate each direction individually
	    //
	    m_ropf[irxn] = calcForwardROP_BV_NoAct(irxn, iBeta, ioc, nStoichElectrons, nu, io);
	    //
	    // Calculate the reverse rate of progress from the difference
	    //
	    m_ropr[irxn] =  m_ropf[irxn] - m_ropnet[irxn];
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
//  This version of takes the electrons out of the reaction rate expression
//  (note: with proper specification of the phase, this shouldn't make a numerical difference (power of 1).
//         But it certainly is a complication and unneeded work)
//   (TODO:  probably can take stoichiometric solids out of the reaction order expression as well.
//           They all contribute powers of 1 as well)
//
void ElectrodeKinetics::determineFwdOrdersBV(ReactionData& rdata, std::vector<doublereal>& fwdFullorders)
{
    //
    //   Start out with the full ROP orders vector.
    //   This vector will have the BV exchange current density orders in it.
    //
    fwdFullorders = rdata.forwardFullOrder_;
    //
    //   forward and reverse beta values
    //
    double betaf = rdata.beta;
    //double betar = 1.0 - betaf;
    //
    //   Loop over the reactants doing away the BV terms.
    //   This should leave the reactant terms only, even if they are non-mass action.
    //
    for (size_t j = 0; j < rdata.reactants.size(); j++) {
	size_t kkin =  rdata.reactants[j];
	double oo = rdata.rstoich[j];
	if (kkin != kElectronIndex_) {
	    fwdFullorders[kkin] += betaf * oo;
	    if (abs(fwdFullorders[kkin]) < 0.00001) {
		fwdFullorders[kkin] = 0.0;
	    }
	} else {
	    fwdFullorders[kkin] = 0.0;
	}
    }
    for (size_t j = 0; j < rdata.products.size(); j++) {
	size_t kkin =  rdata.products[j];
	double oo = rdata.pstoich[j];
	if (kkin != kElectronIndex_) {
	    fwdFullorders[kkin] -= betaf * oo;
	    if (abs(fwdFullorders[kkin]) < 0.00001) {
		fwdFullorders[kkin] = 0.0;
	    }
	} else {
	    fwdFullorders[kkin] = 0.0;
	}
    }
}
//==================================================================================================================
//
//  When the BV form is used we still need to go backwards to calculate the forward rate of progress.
//  This routine does that
//
double ElectrodeKinetics::calcForwardROP_BV(size_t irxn, size_t iBeta, double ioc, double nStoich, double nu, doublereal ioNet)
{ 
    double ropf;
    doublereal rt = GasConstant * thermo(0).temperature();
    //
    //    Calculate gather the exchange current reaction rate constant (where does n_s appear?)
    //
    doublereal beta =  m_beta[iBeta];

#ifdef DEBUG_MODE
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
    const RxnOrders* ro_fwd = m_ctrxn_FwdOrdersList_[iBeta];
    if (ro_fwd == 0) {
	throw CanteraError("ElectrodeKinetics::calcForwardROP_BV()", "forward orders pointer is zero ?!?");
    }
    double tmp = exp(- m_beta[iBeta] * m_deltaG0[irxn] / rt);
    double tmp2 = 1.0;
    const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
    const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
    for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	size_t k = kinSpeciesIDs[j];
	double oo = kinSpeciesOrders[j];
	tmp2 *= pow(m_StandardConc[k], oo);
    }

    //double tmp2 = m_ProdStanConcReac[irxn];
    tmp *= 1.0  / tmp2 / Faraday;
    //
    //  Calculate the chemical reaction rate constant
    //
    double iorc = m_rfn[irxn] * m_perturb[irxn];
    double kf = iorc * tmp;
    //
    //  Calculate the electrochemical factor
    // 
    double eamod = m_beta[iBeta] * deltaElectricEnergy_[irxn];
    kf *= exp(- eamod / rt);
    //
    // Calculate the forward rate of progress
    //   -> get the pointer for the orders
    //
    tmp = 1.0;

    for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	size_t k = kinSpeciesIDs[j];
	double oo = kinSpeciesOrders[j];
	tmp *= pow(m_actConc[k], oo);
    }
    ropf = kf * tmp;
#endif
    //
    //  Now calculate ropf in a separate but equivalent way.
    //  totally equivalent way if resistivity is zero, should be equal (HKM -> Proved exactly in one case)
    //
    double iof = ioc;
    double resistivity = m_ctrxn_resistivity_[iBeta];
    if (fabs(resistivity * ioNet) > fabs(nu)) {
	ioNet = nu / resistivity;
    }
    if (nStoich > 0.0) {
	double exp1 = nStoich * Faraday * beta * (nu - resistivity * ioNet)/ (rt);
	iof *= exp(exp1);
    } else {
#ifdef DEBUG_MODE
	if (ioc > 0) {
	    throw CanteraError(" ", "ioc should be less than zero here");
	}
#endif
	double exp2 = -nu * nStoich * Faraday * (1.0 - beta) / (rt);
        iof = ioc * ( - exp(exp2));
    }
    ropf = iof / ( Faraday * nStoich);

    return ropf;
}
//==================================================================================================================
//
//  When the BV form is used we still need to go backwards to calculate the forward rate of progress.
//  This routine does that
//
double ElectrodeKinetics::calcForwardROP_BV_NoAct(size_t irxn, size_t iBeta, double ioc, double nStoich, double nu,
						  doublereal ioNet)
{
    doublereal TT = thermo(0).temperature();
    doublereal rt = GasConstant * TT;
    //doublereal rrt = 1.0/rt;
    doublereal beta =  m_beta[iBeta];

/*
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
    //   (we don't use standard concentrations at all here);
    //
    double tmp = exp(- m_beta[iBeta] * m_deltaG0[irxn] * rrt);
    double tmp2 = 1.0;
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

	size_t ks = kinSpeciesIDs[j];
	thermo_t& th = speciesPhase(ks);
	size_t n = speciesPhaseIndex(ks);
	size_t klocal = ks - m_start[n];
	double mfS = th.moleFraction(klocal);
	double oo = kinSpeciesOrders[j];
	tmp *= pow(mfS, oo);
    }
    double ropf = kf * tmp;
*/
/*
    if (nStoich > 0) {
	double ropf = ioc / ( Faraday * nStoich);
	double exp1 = nu * nStoich * Faraday * beta / (rt);
        ropf *= exp(exp1);
    } else {
	double ropf = ioc / ( Faraday * nStoich);
	double exp1 = nu * nStoich * Faraday * beta / (rt);
	ropf *= exp(exp1);
    }
*/
    //
    //  With all of the thermo issues, I'm thinking this is the best we can do
    //   (it certainly maintains the forward and reverse rates of progress as being positive)
    //
    double iof = ioc;
    double resistivity = m_ctrxn_resistivity_[iBeta];
    if (fabs(resistivity * ioNet) > fabs(nu)) {
	ioNet = nu / resistivity;
    }
    if (nStoich > 0) {
	double exp1 = nStoich * Faraday * beta * (nu - resistivity * ioNet)/ (rt);
	iof *= exp(exp1);
    } else {
#ifdef DEBUG_MODE
	if (ioc > 0) {
	    throw CanteraError(" ", "ioc should be less than zero here");
	}
#endif
	double exp2 = -nu * nStoich * Faraday * (1.0 - beta) / (rt);
        iof = ioc * ( - exp(exp2));
    }
    double ropf = iof / ( Faraday * nStoich);
    return ropf;
}
//==================================================================================================================
double ElectrodeKinetics::openCircuitVoltage(size_t irxn)
{
    //
    //    Calculate deltaG for all reactions
    //
    getDeltaGibbs(0);
    //
    //    Look up the net number of electrons that are products.
    //
    RxnMolChange*  rmc = rmcVector[irxn];
    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseIndex_];
    double OCV = 0.0;
    if (nStoichElectrons != 0.0) {
        OCV = m_deltaG[irxn] / Faraday / nStoichElectrons;
    }
    return OCV;
}
//==================================================================================================================
//
//  Returns the local exchange current density formulation parameters
//
bool ElectrodeKinetics::
getExchangeCurrentDensityFormulation(size_t irxn,
				     doublereal& nStoichElectrons, doublereal& OCV, doublereal& io,
				     doublereal& overPotential, doublereal& beta,
				     doublereal& resistivity)
{
    size_t iBeta = npos;
    beta = 0.0;
    //
    //  Add logic to handle other reaction types -> return 0 if formulation isn't compatible
    //

    // evaluate rate constants and equilibrium constants at temperature and phi (electric potential)
    _update_rates_T();
    // get updated activities (rates updated below)
    _update_rates_C();

    updateExchangeCurrentQuantities();

    RxnMolChange*   rmc = rmcVector[irxn];
    // could also get this from reactant and product stoichiometry, maybe
    if (metalPhaseIndex_ == npos) {
        nStoichElectrons = 0;
        OCV = 0.0;
        return false;
    } else {
        nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseIndex_];
    }
 

    getDeltaGibbs(0);

    if (nStoichElectrons != 0.0) {
        OCV = m_deltaG[irxn] / Faraday / nStoichElectrons;
    }

    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        if (m_ctrxn[i] == irxn) {
            iBeta = i;
	    break;
        }
    }
    beta = m_beta[iBeta];
 
    doublereal rt = GasConstant*thermo(0).temperature();
  

    double mG0 =  m_deltaG0[irxn];
    int reactionType = m_rxntype[irxn];

    //
    // Start with the forward reaction rate
    //
    double iO = m_rfn[irxn] * m_perturb[irxn];
    int iECDFormulation =  m_ctrxn_ecdf[iBeta];
    if (! iECDFormulation) {
        iO = m_rfn[irxn] * Faraday * nStoichElectrons;
	if (beta > 0.0) {
	    double fac = exp(mG0 / (rt));
	    iO *= pow(fac, beta);
	    //  Need this step because m_rfn includes the inverse of this term, while the formulas
	    //  only use the chemical reaction rate constant.
	    fac = exp( beta * deltaElectricEnergy_[irxn] / (rt));
	    iO *= fac;
	}
    } else {
	iO *= nStoichElectrons;
    }
   
    double omb = 1.0 - beta;
    if (reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {  
	const RxnOrders* ro_fwd = m_ctrxn_ROPOrdersList_[iBeta];
	if (ro_fwd == 0) {
	    throw CanteraError("ElectrodeKinetics::calcForwardROP_BV()", "forward orders pointer is zero ?!?");
	}
	double tmp = 1.0;
	const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
	const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
	for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	    size_t ks = kinSpeciesIDs[j];
	    thermo_t& th = speciesPhase(ks);
	    size_t n = speciesPhaseIndex(ks);
	    size_t klocal = ks - m_start[n];
	    double mfS = th.moleFraction(klocal);
	
	    double oo = kinSpeciesOrders[j];
	    tmp *= pow(mfS, oo);
	}
	iO *= tmp;
    } else if (reactionType == BUTLERVOLMER_RXN) {  
	const RxnOrders* ro_fwd = m_ctrxn_ROPOrdersList_[iBeta];
	if (ro_fwd == 0) {
	    throw CanteraError("ElectrodeKinetics::calcForwardROP_BV()", "forward orders pointer is zero ?!?");
	}
	double tmp = 1.0;
	const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
	const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
	for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	    size_t ks = kinSpeciesIDs[j];
	
	    double oo = kinSpeciesOrders[j];
	    tmp *= pow((m_actConc[ks]/m_StandardConc[ks]), oo);
	}
	iO *= tmp;
    } else {
	for (size_t k = 0; k < m_kk; k++) {
	    doublereal reactCoeff = reactantStoichCoeff(k, irxn);
	    doublereal prodCoeff =  productStoichCoeff(k, irxn);

	    if (reactCoeff != 0.0) {
		iO *= pow(m_actConc[k], reactCoeff*omb);
		iO *= pow(m_StandardConc[k], reactCoeff*beta);
	    }
	    if (prodCoeff != 0.0) {
		iO *= pow(m_actConc[k], prodCoeff*beta);
		iO /= pow(m_StandardConc[k], prodCoeff*omb);
	    }
	}
    }
    io = iO;
    resistivity = m_ctrxn_resistivity_[iBeta];

    double phiMetal = m_thermo[metalPhaseIndex_]->electricPotential();
    double phiSoln = m_thermo[solnPhaseIndex_]->electricPotential();
    double E = phiMetal - phiSoln;
    overPotential = E - OCV;

    return true;
}
//====================================================================================================================
double ElectrodeKinetics::calcCurrentDensity(double nu, double nStoich, double ioc, double beta, double temp, 
					     doublereal resistivity) const
{
     double exp1 = nu * nStoich * Faraday * beta / (GasConstant * temp);
     double exp2 = -nu * nStoich * Faraday * (1.0 - beta) / (GasConstant * temp);
     double val = ioc * (exp(exp1) - exp(exp2));
     if (resistivity > 0.0) {
	 val = solveCurrentRes(nu, nStoich, ioc, beta, temp, resistivity, 0);
     }
     return val;
}
//==================================================================================================================
void ElectrodeKinetics::init()
{
     InterfaceKinetics::init();
     identifyMetalPhase();
}

void ElectrodeKinetics::finalize()
{
    InterfaceKinetics::finalize();
    // Malloc and calculate all of the quantities that go into the extra description of reactions
    rmcVector.resize(m_ii, 0);
    for (size_t i = 0; i < m_ii; i++) {
          rmcVector[i] = new RxnMolChange(this, static_cast<int>(i));
    }
}

//==================================================================================================================

double ElectrodeKinetics::solveCurrentRes(double nu, double nStoich, doublereal ioc, doublereal beta, doublereal temp,
					  doublereal resistivity, int iprob) const
{
    // int nits = 0;
    doublereal f, dfdi, deltai, eexp1, eexp2, exp1, exp2, icurr, deltai_damp;
    doublereal nFRT = nStoich * Faraday /  (GasConstant * temp);
    if (iprob == 0) {
	eexp1 = exp(nu * nFRT * beta);
	eexp2 = exp(-nu * nFRT * (1.0 - beta)) ;

    } else {
	eexp1 = exp(nu * nFRT * beta);
	eexp2 = 0.0;
    }
    icurr = ioc * (eexp1 - eexp2);
    double icurrDamp = icurr;
    if (fabs(resistivity * icurr) > 0.9 * fabs(nu)) {
	icurrDamp = 0.9 * nu / resistivity;
    }
    if (iprob == 0) {
	eexp1 = exp(  nFRT * beta         * (nu - resistivity * icurrDamp));
	eexp2 = exp(- nFRT * (1.0 - beta) * (nu - resistivity * icurrDamp));
    } else {
	eexp1 = exp(  nFRT * beta         * (nu - resistivity * icurrDamp));
	eexp2 = 0.0;
    }
    icurr = ioc * (eexp1 - eexp2);
    if (fabs(resistivity * icurr) > 0.99 * fabs(nu)) {
	icurr = 0.99 * nu / resistivity;
    } 

    do {
	// nits++;
	if (iprob == 0) {
	    exp1 =   nFRT * beta         * (nu - resistivity * icurr);
	    exp2 = - nFRT * (1.0 - beta) * (nu - resistivity * icurr);
	    eexp1 = exp(exp1);
	    eexp2 = exp(exp2);
	    f = icurr - ioc * (eexp1 - eexp2);
	    dfdi = 1.0 - ioc * eexp1 * ( - beta * nFRT * resistivity ) + 
		   ioc * eexp2 * ( (1.0 - beta) * nFRT * resistivity );
	} else {
	    exp1 =   nFRT * beta         * (nu - resistivity * icurr);
	    eexp1 =  exp(exp1);
	    f = icurr - ioc * (eexp1);
	    dfdi = 1.0 - ioc * eexp1 * ( - beta * nFRT * resistivity );
	}
	deltai = - f / dfdi;
	if (fabs(deltai) > 0.1 * fabs(icurr)) {
	    deltai_damp = 0.1 * deltai;
	    if (fabs(deltai_damp) > 0.1 * fabs(icurr)) {
		deltai_damp = 0.1 * icurr * (deltai_damp / fabs(deltai_damp)); 
	    }
	} else if (fabs(deltai) > 0.01 * fabs(icurr))  {
	    deltai_damp = 0.3 * deltai;
	} else if (fabs(deltai) > 0.001 * fabs(icurr))  {
	    deltai_damp = 0.5 * deltai;
	} else {
	    deltai_damp = deltai;
	}
	icurr += deltai_damp;
	if (fabs(resistivity * icurr) > fabs(nu)) {
	    icurr = 0.999 * nu / resistivity;
	}

    } while((fabs(deltai/icurr)> 1.0E-14) && (fabs(deltai) > 1.0E-20));

    // printf("   its = %d\n", nits);

    return icurr;
}
//==================================================================================================================
}
