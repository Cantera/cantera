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


    ///////////////////////////////////////////////////////////
    //
    //    class SurfPhase methods
    //
    ///////////////////////////////////////////////////////////

    SurfPhase::
    SurfPhase(doublereal n0): m_n0(n0), m_tlast(0.0) {
        setNDim(2);
    }

    doublereal SurfPhase::
    enthalpy_mole() const {
        if (m_n0 <= 0.0) return 0.0; 
        _updateThermo();
        return mean_X(m_h0.begin());
    }

    SurfPhase::
    ~SurfPhase() { }

    /**
     * For a surface phase, the pressure is not a relevant
     * thermodynamic variable, and so the enthalpy is equal to the
     * internal energy.
     */
    doublereal SurfPhase::
    intEnergy_mole() const { return enthalpy_mole(); }

    void SurfPhase::
    getStandardChemPotentials(doublereal* mu0) const {
        _updateThermo();
        copy(m_mu0.begin(), m_mu0.end(), mu0);
    }

    void SurfPhase::
    getActivityConcentrations(doublereal* c) const { 
        getConcentrations(c); 
    }

    doublereal SurfPhase::
    standardConcentration(int k) const {
        return m_n0/size(k); 
    }

    doublereal SurfPhase::
    logStandardConc(int k) const {
        return m_logn0 - m_logsize[k];
    }

    void SurfPhase::
    setParameters(int n, doublereal* c) {
        m_n0 = c[0];
        m_logn0 = log(m_n0);
    }

    void SurfPhase::
    initThermo() {
        m_h0.resize(m_kk);
        m_s0.resize(m_kk);
        m_cp0.resize(m_kk);
        m_mu0.resize(m_kk);
        m_work.resize(m_kk);
        m_pe.resize(m_kk, 0.0);
        vector_fp cov(m_kk, 0.0);
        cov[0] = 1.0;
        setCoverages(cov.begin());
        m_logsize.resize(m_kk);
        for (int k = 0; k < m_kk; k++) 
            m_logsize[k] = log(size(k));
    }

    void SurfPhase::
    setPotentialEnergy(int k, doublereal pe) {
        m_pe[k] = pe;
        _updateThermo(true);
    }

    void SurfPhase::
    setSiteDensity(doublereal n0) {
        doublereal x = n0;
        setParameters(1, &x);
    }


    //void SurfPhase::
    //setElectricPotential(doublereal V) {
    //    for (int k = 0; k < m_kk; k++) {
    //        m_pe[k] = charge(k)*Faraday*V;
    //    }
    //    _updateThermo(true);
    //}

    void SurfPhase::
    setCoverages(const doublereal* theta) {
        for (int k = 0; k < m_kk; k++) {
            m_work[k] = m_n0*theta[k]/size(k);
        }
        setConcentrations(m_work.begin());
    }

    void SurfPhase::
    getCoverages(doublereal* theta) const {
        getConcentrations(theta);
        for (int k = 0; k < m_kk; k++) {
            theta[k] *= size(k)/m_n0; 
        }
    }

    void SurfPhase::
    _updateThermo(bool force) const {
        doublereal tnow = temperature();
        if (m_tlast != tnow || force) {
            m_spthermo->update(tnow, m_cp0.begin(), m_h0.begin(), 
                m_s0.begin());
            m_tlast = tnow;
            doublereal rt = GasConstant * tnow;
            int k;
            //doublereal deltaE;
            for (k = 0; k < m_kk; k++) {
                m_h0[k] *= rt;
                m_s0[k] *= GasConstant;
                m_cp0[k] *= GasConstant;
                //deltaE = m_pe[k];
                //m_h0[k] += deltaE;
                m_mu0[k] = m_h0[k] - tnow*m_s0[k];
            }
            m_tlast = tnow;
        }
    }

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
        m_integrator(0),
        m_finalized(false)
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
        doublereal T = thermo().temperature();
        if (T != m_kdata->m_temp || m_redo_rates) {
            doublereal logT = log(T);
            m_rates.update(T, logT, m_kdata->m_rfn.begin());
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
     * reversible reactions.
     */
    void InterfaceKinetics::updateKc() {
        int i, irxn;
        vector_fp& m_rkc = m_kdata->m_rkcn;

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

        fill(m_rkc.begin(), m_rkc.end(), 0.0);

        // compute Delta mu^0 for all reversible reactions
        m_reactantStoich.decrementReactions(m_mu0.begin(), m_rkc.begin()); 
        m_revProductStoich.incrementReactions(m_mu0.begin(), m_rkc.begin());

        for (i = 0; i < m_nrev; i++) {
            irxn = m_revindex[i];
            m_rkc[irxn] = exp(m_rkc[irxn]*rrt);
        }

        for (i = 0; i != m_nirrev; ++i) {
            m_rkc[ m_irrev[i] ] = 0.0;
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

        m_reactantStoich.decrementReactions(m_mu0.begin(), kc); 
        m_revProductStoich.incrementReactions(m_mu0.begin(), kc);
        m_irrevProductStoich.incrementReactions(m_mu0.begin(), kc);

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
        fill(m_rwork.begin(), m_rwork.begin() + m_ii, 0.0);
        m_reactantStoich.decrementReactions(m_pot.begin(), m_rwork.begin()); 
        m_revProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());
        m_irrevProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());

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
        m_reactantStoich.multiply(m_conc.begin(), ropf.begin()); 

        // for reversible reactions, multiply ropr by concentration
        // products
        m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

        // do global reactions
        m_globalReactantStoich.power(m_conc.begin(), ropf.begin());

        for (int j = 0; j != m_ii; ++j) {
            ropnet[j] = ropf[j] - ropr[j];
        }        

        m_kdata->m_ROP_ok = true;
    }

    void InterfaceKinetics::
    addReaction(const ReactionData& r) {

        if (r.reactionType == ELEMENTARY_RXN)      
            addElementaryReaction(r);
        else if (r.reactionType == GLOBAL_RXN)     
            addGlobalReaction(r);

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
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(), 
            r.rateCoeffParameters.begin() );
        // store activation energy
        m_E.push_back(r.rateCoeffParameters[2]);
        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                
        registerReaction( reactionNumber(), ELEMENTARY_RXN, iloc);
    }


    void InterfaceKinetics::
    addGlobalReaction(const ReactionData& r) {
            
        int iloc;
        // install rate coeff calculator
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );

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

        
    void InterfaceKinetics::installReagents(const ReactionData& r) {
            
        m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn
        m_kdata->m_ropr.push_back(0.0);
        m_kdata->m_ropnet.push_back(0.0);
        int n, ns, m;

        int rnum = reactionNumber();

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


    void InterfaceKinetics::installGroups(int irxn, 
        const vector<grouplist_t>& r, const vector<grouplist_t>& p) {
        if (!r.empty()) {
            m_rgroups[reactionNumber()] = r;
            m_pgroups[reactionNumber()] = p;
        }
    }


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
        m_phi.resize(np,0.0);
    }

    void InterfaceKinetics::finalize() {
        m_rwork.resize(nReactions());
        m_finalized = true;
    }


    bool InterfaceKinetics::ready() const {
        return (m_finalized);
    }

    void InterfaceKinetics::
    advanceCoverages(doublereal tstep) {
        if (m_integrator == 0) {
            m_integrator = new ImplicitSurfChem(*this);
            m_integrator->initialize();
        }
        m_integrator->integrate(0.0, tstep);
    }

}








