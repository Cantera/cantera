/**
 *  @file GasKinetics.cpp 
 *
 * Homogeneous kinetics in ideal gases
 * 
 */

// Copyright 2001  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "GasKinetics.h"

#include "ReactionData.h"
#include "StoichManager.h"
#include "Enhanced3BConc.h"
#include "ThirdBodyMgr.h"
#include "RateCoeffMgr.h"

#include <iostream>
using namespace std;

#ifdef HAVE_INTEL_MKL
#include "mkl_vml.h"
#endif

void update_kc(const double* grt, double c0, double* rkc);
void update_rates(double t, double tlog, double* rf);
void mult_by_conc(const double* c, double* ropf, double* ropr);
void eval_ropnet(const double* c, const double* rf, const double* rkc, double* r);

namespace Cantera {


    /**
     * Construct an empty reaction mechanism.
     */    
    GasKinetics::
    GasKinetics(thermo_t* thermo) :
        Kinetics(thermo),
        m_kk(0), 
        m_nfall(0), 
        m_dt_threshold(0.0), // 1.e-6),
        m_nirrev(0), 
        m_nrev(0),
        m_finalized(false)
    {
        m_kdata = new GasKineticsData;
        m_kdata->m_temp = 0.0;
    }

    /**
     * Update temperature-dependent portions of reaction rates and
     * falloff functions.
     */
    void GasKinetics::
    update_T() {}

    void GasKinetics::
    update_C() {}

    void GasKinetics::
    _update_rates_T() {
        doublereal T = thermo().temperature();
        m_kdata->m_logc0 = log(thermo().standardConcentration()); 
        if (fabs(T - m_kdata->m_temp) > m_dt_threshold) {
            doublereal logT = log(T);
            //m_kdata->m_logp0 - logT;
            m_rates.update(T, logT, m_kdata->m_rfn.begin());
            m_falloff_low_rates.update(T, logT, m_kdata->m_rfn_low.begin()); 
            m_falloff_high_rates.update(T, logT, m_kdata->m_rfn_high.begin());
            m_falloffn.updateTemp(T, m_kdata->falloff_work.begin());
            m_kdata->m_temp = T;
            updateKc();
            m_kdata->m_ROP_ok = false;
        }
        else {
            doublereal logT = log(T);
            doublereal dT = T - m_kdata->m_temp;
            //m_kdata->m_logc0 = m_kdata->m_logp0 - logT;
            m_rates.update_dT(T, logT, dT, m_kdata->m_rfn.begin());
            m_falloff_low_rates.update_dT(T, logT, dT, 
                m_kdata->m_rfn_low.begin()); 
            m_falloff_high_rates.update_dT(T, logT, dT, 
                m_kdata->m_rfn_high.begin());
            m_falloffn.updateTemp(T, m_kdata->falloff_work.begin());
            m_kdata->m_temp = T;
            updateKc();
            m_kdata->m_ROP_ok = false;
        }
    };


    /**
     * Update properties that depend on concentrations. Currently only
     * the enhanced collision partner concentrations are updated here.
     */         
    void GasKinetics::
    _update_rates_C() {
        thermo().getActivityConcentrations(m_conc.begin());
        doublereal ctot = thermo().molarDensity();
        m_3b_concm.update(m_conc, ctot, m_kdata->concm_3b_values.begin());
        m_falloff_concm.update(m_conc, ctot, 
            m_kdata->concm_falloff_values.begin());
        m_kdata->m_ROP_ok = false;
    }

    /**
     * Update the equilibrium constants in molar units.
     */
    void GasKinetics::updateKc() {
        int i, irxn;
        vector_fp& m_rkc = m_kdata->m_rkcn;
        
#ifdef HWMECH
        const vector_fp& expg0_RT = thermo().expGibbs_RT();
        doublereal exp_c0 = exp(m_kdata->m_logc0);
        update_kc(expg0_RT.begin(), exp_c0, m_rkc.begin());
#else

        //thermo().getGibbs_RT(m_grt.begin());
        thermo().getStandardChemPotentials(m_grt.begin());
        fill(m_rkc.begin(), m_rkc.end(), 0.0);

        // compute Delta G^0 for all reversible reactions
        m_reactantStoich.decrementReactions(m_grt.begin(), m_rkc.begin()); 
        m_revProductStoich.incrementReactions(m_grt.begin(), m_rkc.begin());
 
        doublereal logc0 = m_kdata->m_logc0;
        doublereal rrt = 1.0/(GasConstant * thermo().temperature());
        for (i = 0; i < m_nrev; i++) {
            irxn = m_revindex[i];
            m_rkc[irxn] = exp(m_rkc[irxn]*rrt - m_dn[irxn]*logc0);
        }

        for(i = 0; i != m_nirrev; ++i) {
            m_rkc[ m_irrev[i] ] = 0.0;
        }
#endif
    }

    /**
     * Get the equilibrium constants of all reactions, whether
     * reversible or not.
     */
    void GasKinetics::getEquilibriumConstants(doublereal* kc) {
        int i;
        _update_rates_T();
        vector_fp& rkc = m_kdata->m_rkcn;
        //thermo().getGibbs_RT(m_grt.begin());
        thermo().getStandardChemPotentials(m_grt.begin());
        fill(rkc.begin(), rkc.end(), 0.0);
        
        // compute Delta G^0 for all reactions
        m_reactantStoich.decrementReactions(m_grt.begin(), rkc.begin()); 
        m_revProductStoich.incrementReactions(m_grt.begin(), 
            rkc.begin());
        m_irrevProductStoich.incrementReactions(m_grt.begin(), 
            rkc.begin());
 
        doublereal logc0 = m_kdata->m_logc0;
        doublereal rrt = 1.0/(GasConstant * thermo().temperature());
        for (i = 0; i < m_ii; i++) {
            kc[i] = exp(-rkc[i]*rrt + m_dn[i]*logc0);
        }
    }


    void GasKinetics::processFalloffReactions() {

        int i;
        const vector_fp& fc = m_kdata->concm_falloff_values;
        const array_fp& m_rf_low = m_kdata->m_rfn_low;
        const array_fp& m_rf_high = m_kdata->m_rfn_high;

        // use m_ropr for temporary storage of reduced pressure
        array_fp& pr = m_kdata->m_ropr;

        array_fp& ropf = m_kdata->m_ropf;

        for (i = 0; i < m_nfall; i++) {
            pr[i] = fc[i] * m_rf_low[i] / m_rf_high[i];
        }

        m_falloffn.pr_to_falloff( pr.begin(), m_kdata->falloff_work.begin() );
        
        for (i = 0; i < m_nfall; i++) {
            pr[i] *= m_rf_high[i]; 
        }

        _scatter_copy(pr.begin(), pr.begin() + m_nfall, 
            ropf.begin(), m_fallindx.begin());
    }


    void GasKinetics::updateROP() {

        _update_rates_T();
        _update_rates_C();

        if (m_kdata->m_ROP_ok) return;

        const vector_fp& rf = m_kdata->m_rfn;
        const vector_fp& m_rkc = m_kdata->m_rkcn;
        array_fp& ropf = m_kdata->m_ropf;
        array_fp& ropr = m_kdata->m_ropr;
        array_fp& ropnet = m_kdata->m_ropnet;

#ifdef HWMECH
        copy(rf.begin(), rf.end(), ropf.begin());
        m_3b_concm.multiply( ropf, m_kdata->concm_3b_values.begin() );
        processFalloffReactions();
        multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());
        eval_ropnet(m_conc.begin(), ropf.begin(), m_rkc.begin(), ropnet.begin());
#else

        // copy rate coefficients into ropf
        copy(rf.begin(), rf.end(), ropf.begin());

        // multiply ropf by enhanced 3b conc for all 3b rxns
        m_3b_concm.multiply( ropf.begin(), m_kdata->concm_3b_values.begin() );

        processFalloffReactions();

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

        for (int j = 0; j != m_ii; ++j) {
            ropnet[j] = ropf[j] - ropr[j];
        }

#endif
        m_kdata->m_ROP_ok = true;
    }


    void GasKinetics::
    addReaction(const ReactionData& r) {

        if (r.reactionType == ELEMENTARY_RXN)      addElementaryReaction(r);
        else if (r.reactionType == THREE_BODY_RXN) addThreeBodyReaction(r);
        else if (r.reactionType == FALLOFF_RXN)    addFalloffReaction(r);

        // operations common to all reaction types
        //installReagents( r.reactants, r.products, r.reversible );
        installReagents( r );
        installGroups(reactionNumber(), r.rgroups, r.pgroups);
        incrementRxnCount();
        m_rxneqn.push_back(r.equation);
    }


    void GasKinetics::
    addFalloffReaction(const ReactionData& r) {

        // install high and low rate coeff calculators

        int iloc = m_falloff_high_rates.install( m_nfall,
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );     
    
        m_falloff_low_rates.install( m_nfall, 
            r.rateCoeffType, r.auxRateCoeffParameters.size(), 
            r.auxRateCoeffParameters.begin() );
                
        // add constant terms to high and low rate
        // coeff value vectors
        m_kdata->m_rfn_high.push_back(r.rateCoeffParameters[0]);
        m_kdata->m_rfn_low.push_back(r.auxRateCoeffParameters[0]);
                
        // add a dummy entry in m_rf, where computed falloff
        // rate coeff will be put
        m_kdata->m_rfn.push_back(0.0);
                
        // add this reaction number to the list of 
        // falloff reactions
        m_fallindx.push_back( reactionNumber() );
                
        // install the enhanced third-body concentration
        // calculator for this reaction
        m_falloff_concm.install( m_nfall, r.thirdBodyEfficiencies, 
            r.default_3b_eff);
                
        // install the falloff function calculator for
        // this reaction
        m_falloffn.install( m_nfall, r.falloffType, r.falloffParameters );
                
        // forward rxn order equals number of reactants, since rate
        // coeff is defined in terms of the high-pressure limit
        m_fwdOrder.push_back(r.reactants.size());

        // increment the falloff reaction counter
        ++m_nfall;
        registerReaction( reactionNumber(), FALLOFF_RXN, iloc);
    }


    void GasKinetics::
    addElementaryReaction(const ReactionData& r) {
        int iloc;

        // install rate coeff calculator
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(), 
            r.rateCoeffParameters.begin() );

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                

        // forward rxn order equals number of reactants
        m_fwdOrder.push_back(r.reactants.size());
        registerReaction( reactionNumber(), ELEMENTARY_RXN, iloc);
    }


    void GasKinetics::
    addThreeBodyReaction(const ReactionData& r) {
            
        int iloc;
        // install rate coeff calculator
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                

        // forward rxn order equals number of reactants + 1
        m_fwdOrder.push_back(r.reactants.size() + 1);

        m_3b_concm.install( reactionNumber(), r.thirdBodyEfficiencies, 
            r.default_3b_eff );
        registerReaction( reactionNumber(), THREE_BODY_RXN, iloc);
    }

        
    void GasKinetics::installReagents(const ReactionData& r) {
        //const vector_int& r,
        //const vector_int& p, bool reversible) {
            
        m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn
        m_kdata->m_ropr.push_back(0.0);
        m_kdata->m_ropnet.push_back(0.0);
        int n, ns, m;

        int rnum = reactionNumber();

        vector_int rk;
        int nr = r.reactants.size();
        for (n = 0; n < nr; n++) {
            ns = r.rstoich[n];
            if (ns != 0) m_rrxn[r.reactants[n]][rnum] += ns;
            for (m = 0; m < ns; m++) {
                rk.push_back(r.reactants[n]);
            }
        }
        m_reactants.push_back(rk);

        vector_int pk;
        int np = r.products.size();
        for (n = 0; n < np; n++) {
            ns = r.pstoich[n];
            if (ns != 0) m_prxn[r.products[n]][rnum] += ns;
            for (m = 0; m < ns; m++) {
                pk.push_back(r.products[n]);
            }
        }
        m_products.push_back(pk);

        m_kdata->m_rkcn.push_back(0.0);
        //        int nr = r.size();
        
        //int i;
        //for (i = 0; i < nr; i++) {
        //     m_rrxn[r[i]][rnum] += 1.0;
        //}
        m_reactantStoich.add( reactionNumber(), rk);

        //int np = p.size();
        //
        //for (i = 0; i < np; i++) {
        //    m_prxn[p[i]][rnum] += 1.0;
        // }

        if (r.reversible) {
            m_revProductStoich.add(reactionNumber(), pk);
            m_dn.push_back(pk.size() - rk.size());
            m_revindex.push_back(reactionNumber());
            m_nrev++;
        }
        else {
            m_irrevProductStoich.add(reactionNumber(), pk);
            m_dn.push_back(pk.size() - rk.size());            
            m_irrev.push_back( reactionNumber() );
            m_nirrev++;
        }        
    }


    void GasKinetics::installGroups(int irxn, 
        const vector<grouplist_t>& r, const vector<grouplist_t>& p) {
        if (!r.empty()) {
            m_rgroups[reactionNumber()] = r;
            m_pgroups[reactionNumber()] = p;
        }
    }


    void GasKinetics::init() { 
        m_kk = thermo().nSpecies();
        m_rrxn.resize(m_kk);
        m_prxn.resize(m_kk);
        m_conc.resize(m_kk);
        m_grt.resize(m_kk);
        m_kdata->m_logp0 = log(thermo().refPressure()) - log(GasConstant);
    }

    void GasKinetics::finalize() {
        if (!m_finalized) {
            int i, j, nr, np;
            m_kdata->falloff_work.resize(m_falloffn.workSize());
            m_kdata->concm_3b_values.resize(m_3b_concm.workSize());
            m_kdata->concm_falloff_values.resize(m_falloff_concm.workSize());

            for (i = 0; i < m_ii; i++) {
                nr = m_reactants[i].size();
                for (j = 0; j < nr; j++) {
                    m_rstoich[i][m_reactants[i][j]]++;
                }
                np = m_products[i].size();
                for (j = 0; j < np; j++) {
                    m_pstoich[i][m_products[i][j]]++;
                }
            }
            m_finalized = true;
        }
    }

    bool GasKinetics::ready() const {
        return (m_finalized);
    }

}








