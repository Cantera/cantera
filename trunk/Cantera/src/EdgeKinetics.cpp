/**
 *  @file EdgeKinetics.cpp 
 *
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "EdgeKinetics.h"
#include "SurfPhase.h"

#include "ReactionData.h"
//#include "StoichManager.h"
#include "RateCoeffMgr.h"

#include <iostream>
using namespace std;


namespace Cantera {

    //////////////////////////////////////////////////////////////////

    /**
     * Construct an empty EdgeKinetics reaction mechanism.
     */    
    EdgeKinetics::
    EdgeKinetics() :
        Kinetics(),
        m_kk(0), 
        m_redo_rates(false),
        m_nirrev(0), 
        m_nrev(0),
        m_finalized(false),
        m_has_electrochem_rxns(false)
    {
        m_kdata = new EdgeKineticsData;
        m_kdata->m_temp = 0.0;
    }

    /**
     * Destructor
     */
    EdgeKinetics::
    ~EdgeKinetics(){
        delete m_kdata;
    }


    /**
     * Update properties that depend on temperature
     *
     */ 
    void EdgeKinetics::
    _update_rates_T() {
        _update_rates_phi();
        doublereal T = thermo(surfacePhaseIndex()).temperature();
        if (T != m_kdata->m_temp || m_redo_rates) {
            m_kdata->m_logtemp = log(T);
            m_rates.update(T, m_kdata->m_logtemp, DATA_PTR(m_kdata->m_rfn));
            if (m_has_electrochem_rxns)
                applyButlerVolmerCorrection(DATA_PTR(m_kdata->m_rfn));
            m_kdata->m_temp = T;
            updateKc();
            m_kdata->m_ROP_ok = false;
            m_redo_rates = false;
        }
    }

    void EdgeKinetics::
    _update_rates_phi() {
        int np = nPhases();
        for (int n = 0; n < np; n++) {
            if (thermo(n).electricPotential() != m_phi[n]) {
                m_phi[n] = thermo(n).electricPotential();
                m_redo_rates = true;
            }
        }
    }


    /**
     * Update properties that depend on concentrations. This method
     * fills out the array of generalized concentrations by calling
     * method getActivityConcentrations for each phase, which classes
     * representing phases should overload to return the appropriate
     * quantities.
     */ 
    void EdgeKinetics::
    _update_rates_C() {
        int n;

        //m_rates.update(m_kdata->m_temp, 
        //    m_kdata->m_logtemp, m_kdata->m_rfn.begin());

        int np = nPhases();
        for (n = 0; n < np; n++) {
            thermo(n).getActivityConcentrations(DATA_PTR(m_conc) + m_start[n]);
        }
        m_kdata->m_ROP_ok = false;
    }


    /**
     * Update the equilibrium constants in molar units for all
     * reversible reactions. Irreversible reactions have their 
     * equilibrium constant set to zero.
     */
    void EdgeKinetics::updateKc() {
        int i, irxn;
        vector_fp& m_rkc = m_kdata->m_rkcn;
        fill(m_rkc.begin(), m_rkc.end(), 0.0);

        if (m_nrev > 0) {

            int n, nsp, k, ik=0;
            doublereal rt = GasConstant*thermo(0).temperature();
            doublereal rrt = 1.0/rt;
            int np = nPhases();
            for (n = 0; n < np; n++) {
                thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
                nsp = thermo(n).nSpecies();
                for (k = 0; k < nsp; k++) {
                    m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
                    m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                    ik++;
                }
            }

            // compute Delta mu^0 for all reversible reactions
            m_reactantStoich.decrementReactions(DATA_PTR(m_mu0), 
                DATA_PTR(m_rkc)); 
            m_revProductStoich.incrementReactions(DATA_PTR(m_mu0), 
                DATA_PTR(m_rkc));

            for (i = 0; i < m_nrev; i++) {
                irxn = m_revindex[i];
                m_rkc[irxn] = exp(m_rkc[irxn]*rrt);
            }
            for (i = 0; i != m_nirrev; ++i) {
                m_rkc[ m_irrev[i] ] = 0.0;
            }
        }
    }



    void EdgeKinetics::checkPartialEquil() {
        int i, irxn;
        vector_fp dmu(nTotalSpecies(), 0.0);
        vector_fp rmu(nReactions(), 0.0);
        if (m_nrev > 0) {

            int n, nsp, k, ik=0;
            doublereal rt = GasConstant*thermo(0).temperature();
            doublereal rrt = 1.0/rt;
            int np = nPhases();
            for (n = 0; n < np; n++) {
                thermo(n).getChemPotentials(DATA_PTR(dmu) + m_start[n]);
                nsp = thermo(n).nSpecies();
                for (k = 0; k < nsp; k++) {
                    dmu[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                    cout << thermo(n).speciesName(k) << "   " << dmu[ik] << endl;
                    ik++;
                }
            }

            // compute Delta mu^ for all reversible reactions
            m_reactantStoich.decrementReactions(DATA_PTR(dmu), DATA_PTR(rmu)); 
            m_revProductStoich.incrementReactions(DATA_PTR(dmu), DATA_PTR(rmu));

            for (i = 0; i < m_nrev; i++) {
                irxn = m_revindex[i];
                cout << "Reaction " << irxn << "  " << exp(rmu[irxn]*rrt) << endl;
            }
        }
    }


    /**
     * Get the equilibrium constants of all reactions, whether
     * reversible or not.
     */
    void EdgeKinetics::getEquilibriumConstants(doublereal* kc) {
        int i;

        int n, nsp, k, ik=0;
        doublereal rt = GasConstant*thermo(0).temperature();
        doublereal rrt = 1.0/rt;
        int np = nPhases();
        for (n = 0; n < np; n++) {
            thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (k = 0; k < nsp; k++) {
                m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
                m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                ik++;
            }
        }

        fill(kc, kc + m_ii, 0.0);

        m_reactantStoich.decrementReactions(DATA_PTR(m_mu0), kc); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_mu0), kc);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_mu0), kc);

        for (i = 0; i < m_ii; i++) {
            kc[i] = exp(-kc[i]*rrt);
        }
    }


    /**
     * For reactions that transfer charge across a potential difference,
     * the activation energies are modified by the potential difference.
     * (see, for example, Baird and Falkner, "Electrochemical Methods"). 
     * This method applies this correction.
     */
    void EdgeKinetics::applyButlerVolmerCorrection(doublereal* kf) {
        int i;

        int n, nsp, k, ik=0;
        doublereal rt = GasConstant*thermo(0).temperature();
        doublereal rrt = 1.0/rt;
        int np = nPhases();

        // compute the electrical potential energy of each species
        for (n = 0; n < np; n++) {
            nsp = thermo(n).nSpecies();
            for (k = 0; k < nsp; k++) {
                m_pot[ik] = Faraday*thermo(n).charge(k)*m_phi[n];
                ik++;
            }
        }

        // compute the change in electrical potential energy for each
        // reaction. This will only be non-zero if a potential
        // difference is present.
        fill(DATA_PTR(m_rwork), DATA_PTR(m_rwork) + m_ii, 0.0);
        m_reactantStoich.decrementReactions(DATA_PTR(m_pot), DATA_PTR(m_rwork)); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_pot), DATA_PTR(m_rwork));
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_pot), DATA_PTR(m_rwork));

        // modify the reaction rates. Only modify those with a
        // non-zero activation energy, and do not decrease the
        // activation energy below zero.
        doublereal ea, eamod;

        int nct = m_beta.size();
        int irxn;
        for (i = 0; i < nct; i++) {
            irxn = m_ctrxn[i];
            eamod = m_beta[i]*m_rwork[irxn];
            //cout << "i, beta = " << i << " " << m_beta[i] << endl;
            if (eamod != 0.0 && m_E[i] != 0.0) {
                ea = GasConstant * m_E[i];
                if (eamod + ea < 0.0) {
                    writelog("Warning: act energy mod too large");
                    eamod = -ea;
                }
                kf[irxn] *= exp(-eamod*rrt);
            }
        }
    }


    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void EdgeKinetics::updateROP() {

        _update_rates_T();
        _update_rates_C();

        if (m_kdata->m_ROP_ok) return;

        const vector_fp& rf = m_kdata->m_rfn;
        const vector_fp& m_rkc = m_kdata->m_rkcn;
        array_fp& ropf = m_kdata->m_ropf;
        array_fp& ropr = m_kdata->m_ropr;
        array_fp& ropnet = m_kdata->m_ropnet;

        // copy rate coefficients into ropf
        copy(rf.begin(), rf.end(), ropf.begin());

        // multiply by perturbation factor
        multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());
           
        // copy the forward rates to the reverse rates                
        copy(ropf.begin(), ropf.end(), ropr.begin());
        
        // for reverse rates computed from thermochemistry, multiply
        // the forward rates copied into m_ropr by the reciprocals of
        // the equilibrium constants
        multiply_each(ropr.begin(), ropr.end(), m_rkc.begin());

        // multiply ropf by concentration products
        m_reactantStoich.multiply(DATA_PTR(m_conc), DATA_PTR(ropf)); 
        
        // for reversible reactions, multiply ropr by concentration
        // products
        m_revProductStoich.multiply(DATA_PTR(m_conc), DATA_PTR(ropr));

        // do global reactions
        //m_globalReactantStoich.power(DATA_PTR(m_conc), ropf.begin());

        for (int j = 0; j != m_ii; ++j) {
            ropnet[j] = ropf[j] - ropr[j];
        }        

        m_kdata->m_ROP_ok = true;
    }


    /**
     * Add a single reaction to the mechanism. This routine
     * must be called after init() and before finalize().
     * This function branches on the types of reactions allowed
     * by the interfaceKinetics manager in order to install
     * the reaction correctly in the manager.
     * The manager allows the following reaction types
     *  Elementary
     *  Surface
     *  Global  
     * There is no difference between elementary and surface 
     * reactions.
     */
    void EdgeKinetics::
    addReaction(const ReactionData& r) {

        int nr = r.reactants.size();

        // a global reaction is idnetified as one with 
        // a reactant stoichiometric coefficient not equal 
        // to the molecularity for some reactant
        bool isglobal = false;
        for (int n = 0; n < nr; n++) {
            if (r.rstoich[n] != int(r.order[n])) {
                isglobal = true; break;
            }
        }
        if (isglobal)
            addGlobalReaction(r);
        else
            addElementaryReaction(r);

        installReagents( r );
        installGroups(reactionNumber(), r.rgroups, r.pgroups);
        incrementRxnCount();
        m_rxneqn.push_back(r.equation);
    }


    void EdgeKinetics::
    addElementaryReaction(const ReactionData& r) {
        int iloc;

        // install rate coeff calculator
        vector_fp rp = r.rateCoeffParameters;

        // coverage dependence
        int ncov = r.cov.size();
        for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);

        iloc = m_rates.install( reactionNumber(), r.rateCoeffType, rp.size(), 
            DATA_PTR(rp) );

        // store activation energy
        if (r.beta > 0.0) {
            m_has_electrochem_rxns = true;
            m_E.push_back(r.rateCoeffParameters[2]);
            m_beta.push_back(r.beta);
            m_ctrxn.push_back(reactionNumber());
        }

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);
                
        registerReaction( reactionNumber(), ELEMENTARY_RXN, iloc);
    }


    void EdgeKinetics::
    addGlobalReaction(const ReactionData& r) {
            
        int iloc;
        // install rate coeff calculator
        vector_fp rp = r.rateCoeffParameters;
        int ncov = r.cov.size();
        for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, rp.size(),
            DATA_PTR(rp) );

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);

        int nr = r.order.size();
        vector_fp ordr(nr);
        for (int n = 0; n < nr; n++) {
            ordr[n] = r.order[n] - r.rstoich[n];
        }
        m_globalReactantStoich.add( reactionNumber(),
            r.reactants, ordr);

        registerReaction( reactionNumber(), GLOBAL_RXN, iloc);
    }

        
    void EdgeKinetics::installReagents(const ReactionData& r) {
            
        m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn
        m_kdata->m_ropr.push_back(0.0);
        m_kdata->m_ropnet.push_back(0.0);
        int n, ns, m;
	doublereal nsFlt;
        int rnum = reactionNumber();

        vector_int rk;
        int nr = r.reactants.size();
        for (n = 0; n < nr; n++) {
            nsFlt = r.rstoich[n];
	    ns = (int) nsFlt;
	    if ((doublereal) ns != nsFlt) {
	      if (ns < 1) ns = 1;
	    }
            m_rrxn[r.reactants[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                rk.push_back(r.reactants[n]);
            }
        }
        m_reactants.push_back(rk);
        
        vector_int pk;
        int np = r.products.size();
        for (n = 0; n < np; n++) {
            nsFlt = r.pstoich[n];
	    ns = (int) nsFlt;
	    if ((doublereal) ns != nsFlt) {
	      if (ns < 1) ns = 1;
	    }
            m_prxn[r.products[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                pk.push_back(r.products[n]);
            }
        }
        m_products.push_back(pk);

        m_kdata->m_rkcn.push_back(0.0);

        m_reactantStoich.add( reactionNumber(), rk);

        if (r.reversible) {
            m_revProductStoich.add(reactionNumber(), pk);
            //m_dn.push_back(pk.size() - rk.size());
            m_revindex.push_back(reactionNumber());
            m_nrev++;
        }
        else {
            m_irrevProductStoich.add(reactionNumber(), pk);
            //m_dn.push_back(pk.size() - rk.size());            
            m_irrev.push_back( reactionNumber() );
            m_nirrev++;
        }        
    }


    void EdgeKinetics::installGroups(int irxn, 
        const vector<grouplist_t>& r, const vector<grouplist_t>& p) {
        if (!r.empty()) {
            m_rgroups[reactionNumber()] = r;
            m_pgroups[reactionNumber()] = p;
        }
    }

    /**
     * Prepare the class for the addition of reactions. This function
     * must be called after instantiation of the class, but before
     * any reactions are actually added to the mechanism.
     * This function calculates m_kk the number of species in all
     * phases participating in the reaction mechanism. We don't know
     * m_kk previously, before all phases have been added. 
     */
    void EdgeKinetics::init() {
        int n;
        m_kk = 0;
        int np = nPhases();
        for (n = 0; n < np; n++) {
            m_kk += thermo(n).nSpecies();
        }
        m_rrxn.resize(m_kk);
        m_prxn.resize(m_kk);
        m_conc.resize(m_kk);
        m_mu0.resize(m_kk);
        m_pot.resize(m_kk, 0.0);
        m_phi.resize(np, 0.0);
    }

    /**
     * Finish adding reactions and prepare for use. This function
     * must be called after all reactions are entered into the mechanism
     * and before the mechanism is used to calculate reaction rates.
     *
     * Here, we resize work arrays based on the number of reactions,
     * since we don't know this number up to now.
     */
    void EdgeKinetics::finalize() {
        m_rwork.resize(nReactions());
        m_finalized = true;
    }


    bool EdgeKinetics::ready() const {
        return (m_finalized);
    }

}








