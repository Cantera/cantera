
// Cantera includes
#include "oneD/Sim1D.h"
#include "oneD/StFlow.h"
#include "oneD/Inlet1D.h"
#include "DenseMatrix.h"
#include "Cabinet.h"
#include "Storage.h"

// Build as a DLL under Windows
#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#else
#define DLL_EXPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999


Cabinet<Sim1D>*   Cabinet<Sim1D>::__storage = 0;
Cabinet<Domain1D>*   Cabinet<Domain1D>::__storage = 0;


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
            return nm.size();
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
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT domain_setBounds(int i, int nl, double* lower, 
        int nu, double* upper) {
        try {
            _domain(i)->setBounds(nl, lower, nu, upper);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT domain_setTolerances(int i, int nr, double* rtol, 
        int na, double* atol, int itime) {
        try {
            _domain(i)->setTolerances(nr, rtol, na, atol, itime);
            return 0;
        }
        catch (CanteraError) { return -1; }
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


    //------------------ inlet domains ------------------------------

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
            ReactingSurf1D* i = new ReactingSurf1D();
            return Cabinet<Domain1D>::cabinet()->add(i);
        }
        catch (CanteraError) { return -1; }
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
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT bdry_massFraction(int i, int k) {
        try {
            return _bdry(i)->massFraction(k);
        }
        catch (CanteraError) { return -1; }
    }

    double DLL_EXPORT bdry_mdot(int i) {
        try {
            return _bdry(i)->mdot();
        }
        catch (CanteraError) { return -1; }
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

    //------------------ stagnation flow domains --------------------

    int DLL_EXPORT stflow_new(int iph, int ikin, int itr) {
        try {
            IdealGasPhase* ph = (IdealGasPhase*)_thermo(iph);
            AxiStagnFlow* x = new AxiStagnFlow(ph, ph->nSpecies(), 2);
            x->setKinetics(*_kinetics(ikin));
            x->setTransport(*_transport(ikin));

            return Cabinet<Domain1D>::cabinet()->add(x);
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
        double* temp) {
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
            for (int n = 0; n < nd; n++) {
                d.push_back(_domain(domains[n]));
            }
            Sim1D* s = new Sim1D(d);
            return Cabinet<Sim1D>::cabinet()->add(s);
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
        int np, double* pos, double* v) {
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

    int DLL_EXPORT sim1D_setTimeStep(int i, double stepsize, int ns, int* nsteps) {
        try {
            _sim1D(i)->setTimeStep(stepsize, ns, nsteps);
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
        double slope, double curve) {
        try {
            _sim1D(i)->setRefineCriteria(dom, ratio, slope, curve);
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
        catch (CanteraError) { return -1.0; }
    }

    double DLL_EXPORT sim1D_workValue(int i, int idom, int icomp, int localPoint) {
        try {
            return _sim1D(i)->workValue(idom, icomp, localPoint);
        }
        catch (CanteraError) { return -1.0; }
    }

    int DLL_EXPORT sim1D_eval(int i, double rdt, int count) {
        try {
            _sim1D(i)->eval(rdt, count);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

}
