/**
 *  @file SurfKinetics.cpp 
 *
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


#include "surfKinetics.h"
#include "ReactionData.h"
#include "RateCoeffMgr.h"
#include "ImplicitSurfChem.h"
#include <iostream>
using namespace std;

#include "ctml.h"
using namespace ctml;

#include <time.h>

namespace Cantera {

    void importInterfaceData(SurfacePhase* ph, SurfKinetics* kin, string fname, string id);

    /**
     * Construct an empty surface reaction mechanism.
     */    
    SurfKinetics::
    SurfKinetics(SurfacePhase* surfphase, 
        thermo_t* th1, 
        thermo_t* th2, string fname, string id) :
        Kinetics(),
        m_surfphase(surfphase),
        m_kk(0),
        m_kk1(0),
        m_kk2(0),
        m_ktot(0),
        m_nirrev(0), 
        m_integrator(0),
        m_finalized(false),
        m_twobulk(false),
        m_xml(new XML_Node("interface_reactions"))
    {

        // add the two bulk phases
        addPhase(*th1);
        if (th2) {
            m_twobulk = true;
            addPhase(*th2);
        }

        m_kk1  = phase(0).nSpecies();

        if (th2) {
            m_kk2 = phase(1).nSpecies();
        }
        m_kk = m_surfphase->nSpecies();
        m_kdata = new SurfKineticsData;
        m_kdata->m_temp = 0.0;

        if (fname != "") importInterfaceData(surfphase, this, fname, id);
    }

    void SurfKinetics::
    _update_rates_T() {
        doublereal T = m_surfphase->temperature();
        if (T != m_kdata->m_temp) {
            doublereal logT = log(T);
            m_rates.update(T, logT, m_kdata->m_rfn.begin());
            m_kdata->m_temp = T;
            m_kdata->m_ROP_ok = false;
        }
    };

    void SurfKinetics::
    _update_rates_C() {
        phase(0).getConcentrations(m_conc.begin());
        if (m_twobulk) {
            phase(1).getConcentrations(m_conc.begin() + m_kk1);
        }
        m_surfphase->getConcentrations(m_conc.begin() + m_kk1 + m_kk2);
        m_rates.update_C(m_conc.begin());
        m_kdata->m_ROP_ok = false;
    }

    void SurfKinetics::updateROP() {

        _update_rates_C();
        _update_rates_T();
        
        if (m_kdata->m_ROP_ok) return;

        const vector_fp& rf = m_kdata->m_rfn;
        vector_fp& ropf = m_kdata->m_ropf;

        // copy rate coefficients into ropf
        copy(rf.begin(), rf.end(), ropf.begin());

        // multiply by perturbation factor
        multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());
           
        // multiply ropf by concentration products
        int i, j, k, o;
        for (i = 0; i < m_ii; i++) {
            for (j = 0; j < m_nr[i]; j++) {
                k = m_reactants[i][j];
                o = m_order[i][j];
                ropf[i] *= pow(m_conc[k],m_order[i][j]);
            }
        }
        m_kdata->m_ROP_ok = true;
    }


    void SurfKinetics::
    getNetProductionRates(doublereal* net) {
        updateROP();
        int i, n, k;
        doublereal q;
        for (k = 0; k < m_ktot; k++) net[k] = 0.0;
        
        for (i = 0; i < m_ii; i++) {
            q = m_kdata->m_ropf[i];
            for (n = 0; n < m_nr[i]; n++) {
                k = m_reactants[i][n];
                net[k] -= q*m_rst[i][n];
            }
            for (n = 0; n < m_np[i]; n++) {
                k = m_products[i][n];
                net[k] += q*m_pst[i][n];
            }
        }
    }

    void SurfKinetics::
    getCreationRates(doublereal* cdot) {
        updateROP();
        int i, n, k;
        doublereal q;
        fill(cdot, cdot + m_ktot, 0.0);
        for (i = 0; i < m_ii; i++) {
            q = m_kdata->m_ropf[i];
            for (n = 0; n < m_np[i]; n++) {
                k = m_products[i][n];
                cdot[k] += q*m_pst[i][n];
            }
        }    
    }

    void SurfKinetics::
    getDestructionRates(doublereal* ddot) {
        updateROP();
        int i, n, k;
        doublereal q;
        fill(ddot, ddot + m_ktot, 0.0);
        for (i = 0; i < m_ii; i++) {
            q = m_kdata->m_ropf[i];
            for (n = 0; n < m_nr[i]; n++) {
                k = m_reactants[i][n];
                ddot[k] += q*m_rst[i][n];
            }
        }
    }

    void SurfKinetics::
    getChemRates(doublereal* rtau) {
        updateROP();
        int i, n, k;
        doublereal q;
        fill(rtau, rtau + m_ktot, 0.0);
        for (i = 0; i < m_ii; i++) {
            q = m_kdata->m_ropf[i];
            for (n = 0; n < m_nr[i]; n++) {
                k = m_reactants[i][n];
                rtau[k] += q*m_rst[i][n];
            }
        }
        for (k = 0;k < m_ktot; k++) {
            if (m_conc[k] != 0.0) 
                rtau[k] = fabs(rtau[k]/m_conc[k]);
            else
                rtau[k] = 0.0;
        }
    }

    void SurfKinetics::
    saveReactionData(
        const vector_int& r, 
        const vector_int& rstoich, 
        const vector_int& order, 
        const vector_int& p, 
        const vector_int& pstoich, 
        const vector_fp& rateParams) {

        if (nReactions() == 0) 
            m_xml->addChild("ReactionArray");

        XML_Node& rxndata = *new XML_Node("reaction");
        int n, k;
        string nm, ph, ustr, comment;
        for (n = 0; n < r.size(); n++) {
            XML_Node& reac = rxndata.addChild("reactant");
            if (r[n] < m_kk1) {
                k = r[n]; 
                nm = phase(0).speciesName(k);
                ph = phase(0).id();
                ustr = "kmol/m^3";
                m_bsp1[nm] = 1;
            }
            else if (r[n] < m_kk1 + m_kk2) {
                k = r[n] - m_kk1; 
                nm = phase(1).speciesName(k);
                ph = phase(1).id();
                ustr = "kmol/m^3";
                m_bsp2[nm] = 1;
            }
            else {
                k = r[n] - m_kk1 - m_kk2; nm = m_surfphase->speciesName(k);
                ph = ""; // m_surfphase->id();
                ustr = "kmol/m^2";
            }
            if (ph != "") reac.addAttribute("phase",ph);
            reac.addAttribute("name",nm);
            reac.addAttribute("stoich",rstoich[n]);
            reac.addAttribute("order",order[n]);
            //            reac.addAttribute("units",ustr);
            comment += nm+" + ";
        }
        comment = comment.substr(0, comment.size() - 2) + " => ";

        for (n = 0; n < p.size(); n++) {
            XML_Node& prod = rxndata.addChild("product");
            if (p[n] < m_kk1) {
                k = p[n]; nm = phase(0).speciesName(k);
                ph = phase(0).id();
                ustr = "kmol/m^3";
            }
            else if (p[n] < m_kk1 + m_kk2) {
                k = p[n] - m_kk1; 
                nm = phase(1).speciesName(k);
                ph = phase(1).id();
                ustr = "kmol/m^3";
            }
            else {
                k = p[n] - m_kk1 - m_kk2; 
                nm = m_surfphase->speciesName(k);
                ph = "";
                ustr = "kmol/m^2";
            }
            if (ph != "") prod.addAttribute("phase",ph);
            prod.addAttribute("name",nm);
            prod.addAttribute("stoich",pstoich[n]);
            comment += nm+" + ";
        }
        comment = "   "+comment.substr(0, comment.size() - 2)+"   ";
        
        XML_Node& rate = rxndata.addChild("rate");
        rate.addAttribute("type","Arrhenius");
        rate.addAttribute("units","kmol/m^2/s");
        addFloat(rate, "A", rateParams[0]);
        addFloat(rate, "n", rateParams[1]);
        addFloat(rate, "E", rateParams[2], "K");

        XML_Node& rxns = m_xml->child("ReactionArray");
        rxns.addComment(comment);
        rxns.addChild(rxndata);
    };


    void SurfKinetics::
    addReaction(const vector_int& r, 
        const vector_int& rstoich, 
        const vector_int& order, 
        const vector_int& p, 
        const vector_int& pstoich, 
        const vector_fp& rateParams) {

        // record reaction parameters
        saveReactionData(r, rstoich, order, p, pstoich, rateParams);

        // prohibit adding more species
        if (!m_surfphase->speciesFrozen())
            m_surfphase->freezeSpecies();

        // if init() hasn't been called yet, call it
        if (m_kk == 0) init();

        int iloc;
        // install rate coeff calculator
        iloc = m_rates.install( m_ii,
            ARRHENIUS, rateParams.size(), rateParams.begin());

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(rateParams[0]);            

        // forward rxn order
        m_order.push_back(order);

        m_kdata->m_ropf.push_back(0.0);     // extend by one for new rxn

        m_reactants.push_back(r);
        m_rst.push_back(rstoich);
        m_products.push_back(p);
        m_pst.push_back(pstoich);

        m_nr.push_back(r.size());
        m_np.push_back(p.size());

        incrementRxnCount();
    }


    void SurfKinetics::init() { 
        m_kk = m_surfphase->nSpecies();
        m_ktot = m_kk + m_kk1 + m_kk2;
        m_conc.resize(m_ktot);
	Kinetics::init();
    }

    void SurfKinetics::save(string fname, string idtag, string comment) {
        struct tm *newtime;
        time_t aclock;
        ::time( &aclock );                  /* Get time in seconds */
        newtime = localtime( &aclock );     /* Convert time to struct tm form */

        ofstream fout(fname.c_str());
        XML_Node root("doc");
        XML_Node& ct = root.addChild("ctml");
        ct.addComment(comment);

        XML_Node& iface = ct.addChild("interface");
        addString(iface,"timestamp",asctime(newtime));
        iface.addAttribute("id",idtag);
        addFloat(iface, "site_density", m_surfphase->siteDensity());
        XML_Node& bp1 = iface.addChild("phase");
        bp1.addAttribute("id",phase(0).id());
        map<string,int>::const_iterator b = m_bsp1.begin(), e = m_bsp1.end();
        for (; b != e; ++b) {
            bp1.addChild("species").addAttribute("name",b->first);
        }
        bp1.addChild(thermo(0).xml());
        if (m_twobulk) {
            XML_Node& bp2 = iface.addChild("phase");
            bp2.addAttribute("id",phase(1).id());
            map<string,int>::const_iterator b = m_bsp2.begin(), e = m_bsp2.end();
            for (; b != e; ++b) {
                bp2.addChild("species").addAttribute("name",b->first);
            }
            bp2.addChild(thermo(1).xml());
        }
        iface.addChild(m_surfphase->xml().child("SpeciesArray"));
        iface.addChild(m_xml->child("ReactionArray"));
        ct.writeHeader(fout);
        ct.write(fout);
        fout.close();
    }

    void SurfKinetics::finalize() {
        if (!m_finalized) {
            m_finalized = true;
        }
    }

    bool SurfKinetics::ready() const {
        return (m_finalized);
    }

    void SurfKinetics::integrate(doublereal dt) {
        finalize();
        if (m_integrator == 0) {
            m_integrator = new ImplicitSurfChem(*this);
            m_integrator->initialize(0.0);
        }
        m_integrator->integrate(0.0, dt);
    }
}








