/**
 *  @file InterfaceKinetics.cpp 
 *
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "InterfaceKinetics.h"
#include "SurfPhase.h"

#include "ReactionData.h"
#include "StoichManager.h"
#include "RateCoeffMgr.h"

#include "ImplicitSurfChem.h"

#include <iostream>
using namespace std;


namespace Cantera {

    //////////////////////////////////////////////////////////////////

    /**
     * Construct an empty InterfaceKinetics reaction mechanism.
     *  @param thermo This is an optional parameter that may be
     *         used to initialize the inherited Kinetics class with
     *         one ThermoPhase class object -> in other words it's
     *         useful for initialization of homogeneous kinetics 
     *         mechanisms.
     */    
    InterfaceKinetics::
    InterfaceKinetics(thermo_t* thermo) :
        Kinetics(thermo),
        m_kk(0), 
        m_redo_rates(false),
        m_nirrev(0), 
        m_nrev(0),
        m_surf(0),
        m_integrator(0),
        m_finalized(false),
        m_has_coverage_dependence(false)
    {
        m_kdata = new InterfaceKineticsData;
        m_kdata->m_temp = 0.0;
    }

    /**
     * Destructor
     */
    InterfaceKinetics::
    ~InterfaceKinetics(){
        delete m_kdata;
        delete m_integrator;
    }


    /**
     * Update properties that depend on temperature
     *
     */ 
    void InterfaceKinetics::
    _update_rates_T() {
        _update_rates_phi();
        if (m_has_coverage_dependence) {
            m_surf->getCoverages(m_conc.begin());
            m_rates.update_C(m_conc.begin());
            m_redo_rates = true;
        }
        doublereal T = thermo(surfacePhaseIndex()).temperature();
        if (T != m_kdata->m_temp || m_redo_rates) {
            m_kdata->m_logtemp = log(T);
            m_rates.update(T, m_kdata->m_logtemp, m_kdata->m_rfn.begin());
            applyButlerVolmerCorrection(m_kdata->m_rfn.begin());
            m_kdata->m_temp = T;
            updateKc();
            m_kdata->m_ROP_ok = false;
            m_redo_rates = false;
        }
    }

    void InterfaceKinetics::
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
    void InterfaceKinetics::
    _update_rates_C() {
        int n;

        int np = nPhases();
        for (n = 0; n < np; n++) {
	  /*
	   * We call the getActivityConcentrations function of each
	   * ThermoPhase class that makes up this kinetics object to 
	   * obtain the generalized concentrations for species within that 
	   * class. This is collected in the vector m_conc. m_start[]
	   * are integer indecises for that vector denoting the start of the
	   * species for each phase.
	   */
	  thermo(n).getActivityConcentrations(m_conc.begin() + m_start[n]);
        }
        m_kdata->m_ROP_ok = false;
    }


    /**
     * Update the equilibrium constants in molar units for all
     * reversible reactions. Irreversible reactions have their 
     * equilibrium constant set to zero.
     */
    void InterfaceKinetics::updateKc() {
        int i, irxn;

        vector_fp& m_rkc = m_kdata->m_rkcn;
        fill(m_rkc.begin(), m_rkc.end(), 0.0);

        static vector_fp mu(nTotalSpecies());
        if (m_nrev > 0) {

            int n, nsp, k, ik=0;
            doublereal rt = GasConstant*thermo(0).temperature();
            doublereal rrt = 1.0/rt;
            int np = nPhases();
            for (n = 0; n < np; n++) {
                thermo(n).getStandardChemPotentials(m_mu0.begin() + m_start[n]);
                nsp = thermo(n).nSpecies();
                for (k = 0; k < nsp; k++) {
                    m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
                    m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                    ik++;
                }
            }

            // compute Delta mu^0 for all reversible reactions
            //m_reactantStoich.decrementReactions(m_mu0.begin(), m_rkc.begin()); 
            //m_revProductStoich.incrementReactions(m_mu0.begin(), m_rkc.begin());
            m_rxnstoich.getRevReactionDelta(m_ii, m_mu0.begin(), m_rkc.begin());

            for (i = 0; i < m_nrev; i++) {
                irxn = m_revindex[i];
                if (irxn < 0 || irxn >= nReactions()) {
                    throw CanteraError("InterfaceKinetics","illegal value: irxn = "+int2str(irxn));
                }
                m_rkc[irxn] = exp(m_rkc[irxn]*rrt);
            }
            for (i = 0; i != m_nirrev; ++i) {
                m_rkc[ m_irrev[i] ] = 0.0;
            }
        }
    }


    void InterfaceKinetics::checkPartialEquil() {
        int i, irxn;
        vector_fp dmu(nTotalSpecies(), 0.0);
        vector_fp rmu(nReactions(), 0.0);
        if (m_nrev > 0) {

            int n, nsp, k, ik=0;
            //doublereal rt = GasConstant*thermo(0).temperature();
            //            doublereal rrt = 1.0/rt;
            int np = nPhases();
            for (n = 0; n < np; n++) {
                thermo(n).getChemPotentials(dmu.begin() + m_start[n]);
                nsp = thermo(n).nSpecies();
                for (k = 0; k < nsp; k++) {
                    dmu[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                    cout << thermo(n).speciesName(k) << "   " << dmu[ik] << endl;
                    ik++;
                }
            }

            // compute Delta mu^ for all reversible reactions
            //m_reactantStoich.decrementReactions(dmu.begin(), rmu.begin()); 
            //m_revProductStoich.incrementReactions(dmu.begin(), rmu.begin());
            m_rxnstoich.getRevReactionDelta(m_ii, dmu.begin(), rmu.begin());

            for (i = 0; i < m_nrev; i++) {
                irxn = m_revindex[i];
                cout << "Reaction " << reactionString(irxn) << "  " << rmu[irxn] << endl;
            }
        }
    }


    /**
     * Get the equilibrium constants of all reactions, whether
     * reversible or not.
     */
    void InterfaceKinetics::getEquilibriumConstants(doublereal* kc) {
        int i;

        int n, nsp, k, ik=0;
        doublereal rt = GasConstant*thermo(0).temperature();
        doublereal rrt = 1.0/rt;
        int np = nPhases();
        for (n = 0; n < np; n++) {
            thermo(n).getStandardChemPotentials(m_mu0.begin() + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (k = 0; k < nsp; k++) {
                m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
                m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                ik++;
            }
        }

        fill(kc, kc + m_ii, 0.0);

        //m_reactantStoich.decrementReactions(m_mu0.begin(), kc); 
        //m_revProductStoich.incrementReactions(m_mu0.begin(), kc);
        //m_irrevProductStoich.incrementReactions(m_mu0.begin(), kc);
        m_rxnstoich.getReactionDelta(m_ii, m_mu0.begin(), kc);

        for (i = 0; i < m_ii; i++) {
            kc[i] = exp(-kc[i]*rrt);
        }
    }


    /**
     * For reactions that transfer charge across a potential difference,
     * the activation energies are modified by the potential difference.
     * (see, for example, ...). This method applies this correction.
     */
    void InterfaceKinetics::applyButlerVolmerCorrection(doublereal* kf) {
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
        //fill(m_rwork.begin(), m_rwork.begin() + m_ii, 0.0);
        //m_reactantStoich.decrementReactions(m_pot.begin(), m_rwork.begin()); 
        //m_revProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());
        //m_irrevProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());
        m_rxnstoich.getReactionDelta(m_ii, m_pot.begin(), m_rwork.begin());

        // modify the reaction rates. Only modify those with a
        // non-zero activation energy, and do not decrease the
        // activation energy below zero.
        doublereal ea, eamod;

        for (i = 0; i < m_ii; i++) {
            eamod = 0.5*m_rwork[i];
            if (eamod != 0.0 && m_E[i] != 0.0) {
                ea = GasConstant * m_E[i];
                if (eamod + ea < 0.0) eamod = -ea;
                kf[i] *= exp(-eamod*rrt);
            }
        }
    }




    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::getFwdRateConstants(doublereal* kfwd) {

        _update_rates_T();
        _update_rates_C();

        const vector_fp& rf = m_kdata->m_rfn;

        // copy rate coefficients into kfwd
        copy(rf.begin(), rf.end(), kfwd);

        // multiply by perturbation factor
        multiply_each(kfwd, kfwd + nReactions(), m_perturb.begin());
           
    }


    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::getRevRateConstants(doublereal* krev) {
        getFwdRateConstants(krev);
        const vector_fp& rkc = m_kdata->m_rkcn;
        multiply_each(krev, krev + nReactions(), rkc.begin());
    }

    void InterfaceKinetics::getActivationEnergies(doublereal *E) {
        copy(m_E.begin(), m_E.end(), E);
    }

    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::updateROP() {

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
        m_rxnstoich.multiplyReactants(m_conc.begin(), ropf.begin()); 
        //m_reactantStoich.multiply(m_conc.begin(), ropf.begin()); 

        // for reversible reactions, multiply ropr by concentration
        // products
        m_rxnstoich.multiplyRevProducts(m_conc.begin(), ropr.begin()); 
        //m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

        // do global reactions
        //m_globalReactantStoich.power(m_conc.begin(), ropf.begin());

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
    void InterfaceKinetics::
    addReaction(const ReactionData& r) {

        //        int nr = r.reactants.size();

        // a global reaction is identified as one with 
        // a reactant stoichiometric coefficient not equal 
        // to the molecularity for some reactant
//         bool isglobal = false;
//         for (int n = 0; n < nr; n++) {
//            if (r.rstoich[n] != int(r.order[n])) {
//                isglobal = true; break;
//            }
//         }
        // if (isglobal)
        //     addGlobalReaction(r);
        //else

        //        if (r.global) 
        //    cout << r.equation << " is global " << endl;

        addElementaryReaction(r);

        //if (r.global) {
        //    int nr = r.order.size();
        //    vector_fp ordr(nr);
        //    for (int n = 0; n < nr; n++) {
        //        ordr[n] = r.order[n] - r.rstoich[n];
                //      cout << r.reactants[n] << "  " << r.order[n] << "  " << ordr[n] << endl;
        //  }
            //m_globalReactantStoich.add( reactionNumber(),
            //    r.reactants, ordr);
            //}


//         if (r.reactionType == ELEMENTARY_RXN)      
//             addElementaryReaction(r);
//         if (r.reactionType == SURFACE_RXN)      
//             addElementaryReaction(r);
//         else if (r.reactionType == GLOBAL_RXN)     
//             addGlobalReaction(r);

        // operations common to all reaction types
        installReagents( r );
        installGroups(reactionNumber(), r.rgroups, r.pgroups);
        incrementRxnCount();
        m_rxneqn.push_back(r.equation);
    }


    void InterfaceKinetics::
    addElementaryReaction(const ReactionData& r) {
        int iloc;
        // install rate coeff calculator
        vector_fp rp = r.rateCoeffParameters;
        int ncov = r.cov.size();
        if (ncov > 3) {
            m_has_coverage_dependence = true;
        }
        for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, rp.size(), 
            rp.begin() );
        // store activation energy
        m_E.push_back(r.rateCoeffParameters[2]);
        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                
        registerReaction( reactionNumber(), ELEMENTARY_RXN, iloc);
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
//         m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);

//         int nr = r.order.size();
//         vector_fp ordr(nr);
//         for (int n = 0; n < nr; n++) {
//             ordr[n] = r.order[n] - r.rstoich[n];
//         }
//         m_globalReactantStoich.add( reactionNumber(),
//             r.reactants, ordr);

//         registerReaction( reactionNumber(), GLOBAL_RXN, iloc);
//     }

        
    void InterfaceKinetics::installReagents(const ReactionData& r) {
            
        m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn
        m_kdata->m_ropr.push_back(0.0);
        m_kdata->m_ropnet.push_back(0.0);
        int n, ns, m;

        int rnum = reactionNumber();

        // vectors rk and pk are lists of species numbers, with
        // repeated entries for species with stoichiometric
        // coefficients > 1. This allows the reaction to be defined
        // with unity reaction order for each reactant, and so the
        // faster method 'multiply' can be used to compute the rate of
        // progress instead of 'power'.

        // Note that this procedure is used for global reactions also. 
        // The 

        vector_int rk;
        int nr = r.reactants.size();
        for (n = 0; n < nr; n++) {
            ns = r.rstoich[n];
            m_rrxn[r.reactants[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                rk.push_back(r.reactants[n]);
            }
        }
        m_reactants.push_back(rk);
        
        vector_int pk;
        int np = r.products.size();
        for (n = 0; n < np; n++) {
            ns = r.pstoich[n];
            m_prxn[r.products[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                pk.push_back(r.products[n]);
            }
        }
        m_products.push_back(pk);

        m_kdata->m_rkcn.push_back(0.0);

        m_rxnstoich.add( reactionNumber(), r);

        //m_reactantStoich.add( reactionNumber(), rk);

        if (r.reversible) {
            //  m_revProductStoich.add(reactionNumber(), pk);
            m_revindex.push_back(reactionNumber());
            m_nrev++;
        }
        else {
            //m_irrevProductStoich.add(reactionNumber(), pk);
            m_irrev.push_back( reactionNumber() );
            m_nirrev++;
        }        
    }


    void InterfaceKinetics::installGroups(int irxn, 
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
    void InterfaceKinetics::init() {
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
    void InterfaceKinetics::finalize() {
        m_rwork.resize(nReactions());
        int ks = surfacePhaseIndex();
        if (ks < 0) throw CanteraError("InterfaceKinetics::finalize",
             "no surface phase is present.");
        m_surf = (SurfPhase*)&thermo(ks);
        m_finalized = true;
    }


    bool InterfaceKinetics::ready() const {
        return (m_finalized);
    }

    void InterfaceKinetics::
    advanceCoverages(doublereal tstep) {
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

}








