/**
 * @file ctonedim.cpp
 */
/*
 *      $Id: ctonedim.cpp,v 1.25 2009/07/11 17:16:09 hkmoffa Exp $
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#pragma warning(disable:4800)
#endif

#define CANTERA_USE_INTERNAL
#include "ctonedim.h"


// Cantera includes
#include "config.h"
#include "Sim1D.h"
#include "StFlow.h"
#include "Inlet1D.h"
#include "DenseMatrix.h"



// local includes
#include "Cabinet.h"
#include "Storage.h"

using namespace std;
using namespace Cantera;

template<> Cabinet<Sim1D>*   Cabinet<Sim1D>::__storage = 0;
template<> Cabinet<Domain1D>*   Cabinet<Domain1D>::__storage = 0;


inline Sim1D* _sim1D(int i) {
    return Cabinet<Sim1D>::cabinet()->item(i);
}

inline Domain1D* _domain(int i) {
    return Cabinet<Domain1D>::cabinet()->item(i);
}

static StFlow* _stflow(int i) {
    Domain1D* d = _domain(i);
    if (d->domainType() == cFlowType) return (StFlow*)d;
    else
        throw CanteraError("_stflow","wrong domain type");
}

static Bdry1D* _bdry(int i) {
    Domain1D* d = _domain(i);
    if (d->isConnector()) return (Bdry1D*)d;
    else
        throw CanteraError("_bdry","wrong domain type: "
            +int2str(d->domainType()));
}

inline ThermoPhase* _phase(int n) {
    return Storage::__storage->__thtable[n];
}

inline Kinetics* _kinetics(int n) {
    return Storage::__storage->__ktable[n];
}

inline ThermoPhase* _thermo(int n) {
    return Storage::__storage->__thtable[n];
}

inline Transport* _transport(int n) {
    return Storage::__storage->__trtable[n];
}


extern "C" {

    int DLL_EXPORT domain_clear() {
        try {
            Cabinet<Domain1D>::cabinet()->clear();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT domain_del(int i) {
        Cabinet<Domain1D>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT domain_type(int i) {
        return _domain(i)->domainType();
    }

    int DLL_EXPORT domain_index(int i) {
        return _domain(i)->domainIndex();
    }

    int DLL_EXPORT domain_nComponents(int i) {
        return _domain(i)->nComponents();
    }

    int DLL_EXPORT domain_nPoints(int i) {
        return _domain(i)->nPoints();
    }


    int DLL_EXPORT domain_componentName(int i, int n, int sz, char* buf) {
        try {
            string nm = _domain(i)->componentName(n);
            int lout = min(sz, nm.size());
            copy(nm.c_str(), nm.c_str() + lout, buf);
            buf[lout] = '\0';
            return static_cast<int>(nm.size());
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT domain_componentIndex(int i, char* name) {
        try {
            int n = _domain(i)->componentIndex(string(name));
            return n;
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT domain_grid(int i, int n) {
        try {
            return _domain(i)->grid(n);
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT domain_setBounds(int i, int n, double lower, double upper) {
        try {
            _domain(i)->setBounds(n, lower, upper);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT domain_upperBound(int i, int n) {
        try {
            return _domain(i)->upperBound(n);
        }
        catch (CanteraError) { return DERR; }
    }

    double DLL_EXPORT domain_lowerBound(int i, int n) {
        try {
            return _domain(i)->lowerBound(n);
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT domain_setTolerances(int i, int n, double rtol, 
        double atol, int itime) {
        try {
            _domain(i)->setTolerances(n, rtol, atol, itime);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT domain_rtol(int i, int n) {
        try {
            return _domain(i)->rtol(n);
        }
        catch (CanteraError) { return DERR; }
    }

    double DLL_EXPORT domain_atol(int i, int n) {
        try {
            return _domain(i)->atol(n);
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT domain_setupGrid(int i, int npts, double* grid) {
        try {
            _domain(i)->setupGrid(npts, grid);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT domain_setID(int i, char* id) {
        try {
            string s = string(id);
            _domain(i)->setID(s);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT domain_setDesc(int i, char* desc) {
        try {
            string s = string(desc);
            _domain(i)->setDesc(s);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT inlet_new() {
        try {
            Inlet1D* i = new Inlet1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT surf_new() {
        try {
            Surf1D* i = new Surf1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT reactingsurf_new() {
        try {
	  //writelog("in reactingsurf_new\n");
            Domain1D* i = new ReactingSurf1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { writelog("error"); return -1; }
    }

    int DLL_EXPORT symm_new() {
        try {
            Symm1D* i = new Symm1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT outlet_new() {
        try {
            Outlet1D* i = new Outlet1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT outletres_new() {
        try {
            OutletRes1D* i = new OutletRes1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT bdry_setMdot(int i, double mdot) {
        try {
            _bdry(i)->setMdot(mdot);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT bdry_setTemperature(int i, double t) {
        try {
            _bdry(i)->setTemperature(t);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT bdry_setMoleFractions(int i, char* x) {
        try {
            _bdry(i)->setMoleFractions(string(x));
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT bdry_temperature(int i) {
        try {
            return _bdry(i)->temperature();
        }
        catch (CanteraError) { return DERR; }
    }

    double DLL_EXPORT bdry_massFraction(int i, int k) {
        try {
            return _bdry(i)->massFraction(k);
        }
        catch (CanteraError) { return DERR; }
    }

    double DLL_EXPORT bdry_mdot(int i) {
        try {
            return _bdry(i)->mdot();
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT reactingsurf_setkineticsmgr(int i, int j) {
        try {
            ReactingSurf1D* srf = (ReactingSurf1D*)_bdry(i);
            InterfaceKinetics* k = (InterfaceKinetics*)_kinetics(j);
            srf->setKineticsMgr(k);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT reactingsurf_enableCoverageEqs(int i, int onoff) {
        try {
            ReactingSurf1D* srf = (ReactingSurf1D*)_bdry(i);
            srf->enableCoverageEquations(bool(onoff));
            return 0;
        }

        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT inlet_setSpreadRate(int i, double v) {
        try {
            Inlet1D* inlt = (Inlet1D*)_bdry(i);
            inlt->setSpreadRate(v);
            return 0;
        }

        catch (CanteraError) { return -1; }
    }

    //------------------ stagnation flow domains --------------------

    int DLL_EXPORT stflow_new(int iph, int ikin, int itr, int itype) {
        try {
            IdealGasPhase* ph = (IdealGasPhase*)_thermo(iph);
            if (itype == 1) {
                AxiStagnFlow* x = new AxiStagnFlow(ph, ph->nSpecies(), 2);
                x->setKinetics(*_kinetics(ikin));
                x->setTransport(*_transport(itr));
                return Cabinet<Domain1D>::cabinet()->add(x);
            }
            else if (itype == 2) {
                FreeFlame* x = new FreeFlame(ph, ph->nSpecies(), 2);
                x->setKinetics(*_kinetics(ikin));
                x->setTransport(*_transport(itr));
                return Cabinet<Domain1D>::cabinet()->add(x);
            }
            else {
                return -2;
            }
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT stflow_setTransport(int i, int itr, int iSoret) {
        bool withSoret = false;
        if (iSoret > 0) withSoret = true;
        try {
            _stflow(i)->setTransport(*_transport(itr), withSoret);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT stflow_enableSoret(int i, int iSoret) {
        bool withSoret = false;
        if (iSoret > 0) withSoret = true;
        try {
            _stflow(i)->enableSoret(withSoret);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT stflow_setPressure(int i, double p) {
        try {
            _stflow(i)->setPressure(p);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT stflow_setFixedTempProfile(int i, int n, double* pos, 
        int m, double* temp) {
        try {
            int j;
            vector_fp vpos(n), vtemp(n);
            for (j = 0; j < n; j++) {
                vpos[j] = pos[j]; 
                vtemp[j] = temp[j];
            }
            _stflow(i)->setFixedTempProfile(vpos, vtemp);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT stflow_solveSpeciesEqs(int i, int flag) {
        try {
            if (flag > 0) 
                _stflow(i)->solveSpecies(-1);
            else
                _stflow(i)->fixSpecies(-1);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT stflow_solveEnergyEqn(int i, int flag) {
        try {
            if (flag > 0)
                _stflow(i)->solveEnergyEqn(-1);
            else
                _stflow(i)->fixTemperature(-1);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    //------------------- Sim1D --------------------------------------

    int DLL_EXPORT sim1D_new(int nd, int* domains) {
        vector<Domain1D*> d;
        try {
	  //  cout << "nd = " << nd << endl;
            for (int n = 0; n < nd; n++) {
	      //writelog("n = "+int2str(n)+"\n");
	      //writelog("dom = "+int2str(domains[n])+"\n");
                d.push_back(_domain(domains[n]));
            }
            //writelog("in sim1D_new, calling new Sim1D\n");
            Sim1D* s = new Sim1D(d);
            //writelog("in sim1D_new, ret Sim1D\n");
            return Cabinet<Sim1D>::cabinet()->add(s);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_clear() {
        try {
            Cabinet<Sim1D>::cabinet()->clear();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_del(int i) {
        Cabinet<Sim1D>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT sim1D_setValue(int i, int dom, int comp, 
        int localPoint, double value) {
        try {
            _sim1D(i)->setValue(dom, comp, localPoint, value);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setProfile(int i, int dom, int comp, 
        int np, double* pos, int nv, double* v) {
        try {
            vector_fp vv, pv;
            for (int n = 0; n < np; n++) {
                vv.push_back(v[n]);
                pv.push_back(pos[n]);
            }
            _sim1D(i)->setProfile(dom, comp, pv, vv);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setFlatProfile(int i, int dom, int comp, double v) {
        try {
            _sim1D(i)->setFlatProfile(dom, comp, v);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_showSolution(int i, char* fname) {
        string fn = string(fname);
        if (fn == "-")
            _sim1D(i)->showSolution();
        else {
            ofstream fout(fname);
            _sim1D(i)->showSolution(fout);
            fout.close();
        }
        return 0;
    }

    int DLL_EXPORT sim1D_setTimeStep(int i, double stepsize, int ns, integer* nsteps) {
        try {
            _sim1D(i)->setTimeStep(stepsize, ns, nsteps);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_getInitialSoln(int i) {
        try {
            _sim1D(i)->getInitialSoln();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_solve(int i, int loglevel, int refine_grid) {
        try {
            bool r = (refine_grid == 0 ? false : true);
            _sim1D(i)->solve(loglevel, r);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_refine(int i, int loglevel) {
        try {
            _sim1D(i)->refine(loglevel);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setRefineCriteria(int i, int dom, double ratio,
        double slope, double curve, double prune) {
        try {
            _sim1D(i)->setRefineCriteria(dom, ratio, slope, curve, prune);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_save(int i, char* fname, char* id, 
        char* desc) {
        try {
            string sname = string(fname);
            string sid = string(id);
            string sdesc = string(desc);
            _sim1D(i)->save(sname, sid, sdesc);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_restore(int i, char* fname, char* id) {
        try {
            string sname = string(fname);
            string sid = string(id);
            _sim1D(i)->restore(sname, sid);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_writeStats(int i) {
        try {
            _sim1D(i)->writeStats();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_domainIndex(int i, char* name) {
        try {
            return _sim1D(i)->domainIndex(string(name));
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT sim1D_value(int i, int idom, int icomp, int localPoint) {
        try {
            return _sim1D(i)->value(idom, icomp, localPoint);
        }
        catch (CanteraError) { return DERR; }
    }

    double DLL_EXPORT sim1D_workValue(int i, int idom, int icomp, int localPoint) {
        try {
            return _sim1D(i)->workValue(idom, icomp, localPoint);
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT sim1D_eval(int i, double rdt, int count) {
        try {
            _sim1D(i)->eval(rdt, count);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setMaxJacAge(int i, int ss_age, int ts_age) {
        try {
            _sim1D(i)->setJacAge(ss_age, ts_age);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_timeStepFactor(int i, double tfactor) {
        try {
            _sim1D(i)->setTimeStepFactor(tfactor);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setTimeStepLimits(int i, double tsmin, double tsmax) {
        try {
            if (tsmin > 0.0)
                _sim1D(i)->setMinTimeStep(tsmin);
            if (tsmax > 0.0)
                _sim1D(i)->setMaxTimeStep(tsmax);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_setFixedTemperature(int i, double temp) {
        try {
            _sim1D(i)->setFixedTemperature(temp);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT sim1D_evalSSJacobian(int i) {
        try {
            _sim1D(i)->evalSSJacobian();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT sim1D_jacobian(int i, int m, int n) {
        try {
            return _sim1D(i)->jacobian(m,n);
        }
        catch (CanteraError) { return DERR; }
    }

    int DLL_EXPORT sim1D_size(int i) {
        try {
            return _sim1D(i)->size();
        }
        catch (CanteraError) { return -1; }
    }

}
