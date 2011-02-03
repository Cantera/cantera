
// Cantera includes
#include "oneD/OneDim.h"
#include "oneD/StFlow.h"
#include "oneD/Inlet1D.h"
#include "oneD/MultiNewton.h"
#include "DenseMatrix.h"
#include "Cabinet.h"
#include "Storage.h"

// Build as a DLL under Windows
#ifdef WIN32
#ifdef NO_DLL_BUILD
#define DLL_EXPORT
#else
#define DLL_EXPORT __declspec(dllexport)
#endif
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#else
#define DLL_EXPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999

using namespace FlowBdry;

Cabinet<OneDim>*   Cabinet<OneDim>::__storage = 0;
Cabinet<StFlow>*   Cabinet<StFlow>::__storage = 0;
Cabinet<Boundary>* Cabinet<Boundary>::__storage = 0;
//Cabinet<Surf1D>* Cabinet<Surf1D>::__storage = 0;

inline OneDim* _onedim(int i) {
    return Cabinet<OneDim>::cabinet()->item(i);
}

inline StFlow* _flow(int i) {
    return Cabinet<StFlow>::cabinet()->item(i);
}

inline Boundary* _boundary(int i) {
    return Cabinet<Boundary>::cabinet()->item(i);
}

inline Bdry1D* _bndry(int i) {
    return Cabinet<Bdry1D>::cabinet()->item(i);
}

//inline SurfKinetics* _surfkin(int i) {
//    return Cabinet<SurfKinetics>::cabinet()->item(i);
//}

//inline Surf1D* _surface(int i) {
//    return Cabinet<Surf1D>::cabinet()->item(i);
//}

inline DenseMatrix* _matrix(int i) {
    return Cabinet<DenseMatrix>::cabinet()->item(i);
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

    int DLL_EXPORT flow_new(int type, int iph, int np) {
        IdealGasPhase* ph = (IdealGasPhase*)_thermo(iph);
        StFlow* x;
        try {
            switch (type) {
            case 0:
                x = new AxiStagnFlow(ph, ph->nSpecies(), np); break;
            case 1:
                x = new OneDFlow(ph, ph->nSpecies(), np); break;
            default:
                return -2;
            }
            return Cabinet<StFlow>::cabinet()->add(x);
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT flow_del(int i) {
        Cabinet<StFlow>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT flow_copy(int i) {
        return Cabinet<StFlow>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT flow_assign(int i, int j) {
        return Cabinet<StFlow>::cabinet()->assign(i,j);
    }

//     int DLL_EXPORT flow_readinputs(int i, char* infile) {
//         try {
//             ifstream f(infile);
//             if (!f) throw CanteraError("flow_readinputs",
//                 "error opening input file");
//             //            _flow(i)->readInputs(f);
//             f.close();
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//         catch (...) { return ERR; }
//     }

    int DLL_EXPORT flow_setupgrid(int i, int npts, double* grid) {
        try {
            _flow(i)->setupGrid(npts, grid);
            return 0;
        }
        catch (CanteraError) { return -1; }
        //catch (...) { return ERR; }
    }

    int DLL_EXPORT flow_setthermo(int i, int k) {
        IdealGasPhase* th = (IdealGasPhase*)_thermo(k);
        _flow(i)->setThermo(*th);
        return 0;
    }

    int DLL_EXPORT flow_setkinetics(int i, int k) {
        Kinetics* kin = _kinetics(k);
        _flow(i)->setKinetics(*kin);
        return 0;
    }

    int DLL_EXPORT flow_settransport(int i, int k, int soret) {
        try {
            Transport* tr = _transport(k);
            bool withSoret = (soret == 1); 
            _flow(i)->setTransport(*tr, withSoret);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT flow_settemperature(int i, int j, double t) {
        _flow(i)->setTemperature(j, t);
        return 0;
    }

    int DLL_EXPORT flow_setmassfraction(int i, int j, int k, double t) {
        _flow(i)->setMassFraction(j, k, t);
        return 0;
    }

    int DLL_EXPORT flow_setpressure(int i, double p) {
        _flow(i)->setPressure(p);
        return 0;
    }

    int DLL_EXPORT flow_showsolution(int i, char* fname, double* soln) {
        string fn = string(fname);
        if (fn == "-")
            _flow(i)->showSolution(cout, soln);
        else {
            ofstream fout(fname);
            _flow(i)->showSolution(fout, soln);
            fout.close();
        }
        return 0;
    }

    int DLL_EXPORT flow_outputtec(int i, doublereal* x, 
        char* fname, char* title, int zone) {
        ofstream f(fname);
        //DenseMatrix* mat = _matrix(m);
        _flow(i)->outputTEC(f, x, string(title), zone);
        return 0;
    }


    // solve / fix

    int DLL_EXPORT flow_solveenergyeqn(int i, int j) {
        _flow(i)->solveEnergyEqn(j);
        return 0;
    }

    int DLL_EXPORT flow_fixtemperature(int i, int j) {
        _flow(i)->fixTemperature(j);
        return 0;
    }

    int DLL_EXPORT flow_setenergyfactor(int i, double e) {
        _flow(i)->setEnergyFactor(e);
        return 0;
    }

    int DLL_EXPORT flow_fixspecies(int i, int j) {
        _flow(i)->fixSpecies(j);
        return 0;
    }

    int DLL_EXPORT flow_solvespecies(int i, int j) {
        _flow(i)->solveSpecies(j);
        return 0;
    }

    int DLL_EXPORT flow_resize(int i, int points) {
        _flow(i)->resize(points);
        return 0;
    }

//     int DLL_EXPORT flow_integratechem(int i, doublereal* x, double dt) {
//         try{
//             _flow(i)->integrateChem(x, dt);
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//     }

    int DLL_EXPORT flow_settolerances(int i, int nr, 
        doublereal* rtol, int na, doublereal* atol) {
        try {
            _flow(i)->setTolerances(nr, rtol, na, atol);
            return 0;
        }
        catch (CanteraError) { return -1; }
        //catch (...) { return ERR; }
    }

    int DLL_EXPORT flow_eval(int i, int j, doublereal* x, doublereal* r, integer* m) {
        try {
            _flow(i)->eval(j, x, r, m);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT flow_restore(int i, int job, char* fname, char* id, 
        int& size_z, doublereal* z, int& size_soln, doublereal* soln) {
        try {
            _flow(i)->restore(job, fname, string(id), size_z, z, 
                size_soln, soln);
            return 0;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }

    int DLL_EXPORT flow_setfixedpoint(int i, int j0, doublereal t0) {
        _flow(i)->setFixedPoint(j0, t0);
        return 0;
    }


    int DLL_EXPORT flow_setboundaries(int i, int nleft, int nright) {
        Boundary *left=0, *right=0;
        if (nleft > 0) left = _boundary(nleft);
        if (nright > 0) right = _boundary(nright);
        _flow(i)->setBoundaries(left, right);
        return 0;
    }


    //==========================================================

    int DLL_EXPORT bdry_new(int type, int iph, int kin) {
        Boundary* x=0;
        //const doublereal* wt = _phase(iph)->molecularWeights().begin();
        int nsp = _phase(iph)->nSpecies();
        switch (type) {
        case 0:
            x = new Inlet(nsp); break;
        case 1:
            x = new Outlet(nsp); break;
            //case 2:
            //if (kin > 0)
            //    x = new Surface(nsp, _surfkin(kin));
            //else
            //    x = new Surface(nsp, 0);
            //break;
        case 3:
            x = new SymmPlane(nsp); break;
        default:
            return -2;
        }
        return Cabinet<Boundary>::cabinet()->add(x);
    }


    int DLL_EXPORT bdry_del(int i) {
        Cabinet<Boundary>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT bdry_copy(int i) {
        return Cabinet<Boundary>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT bdry_assign(int i, int j) {
        return Cabinet<Boundary>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT bdry_set(int i, int n, doublereal* v) {
        switch (n) {
        case 1:
            _boundary(i)->set_mdot(*v); break;
        case 2:
            _boundary(i)->set_V(*v); break;
        case 3:
            _boundary(i)->set_T(*v); break;
        case 4:
            _boundary(i)->set_Y(v); break;
        default:
            throw CanteraError("bdry_set","unknown option");
        }
        return 0;
    }

    //=========================================================


    int DLL_EXPORT onedim_new(int nd, int* domains, int* types) {
        int i;
        vector<Domain1D*> doms;
        for (i = 0; i < nd; i++) {
            switch (types[i]) {
            case 0:
                doms.push_back(_flow(domains[i])); break;
                //case 1:
                //doms.push_back(_surface(domains[i])); break;
            case 2:
                doms.push_back(_bndry(domains[i])); break;
            default:
                throw CanteraError("onedim_new", "unknown domain type");
            }
        }
        try {
            OneDim* x = new OneDim(doms);
            return Cabinet<OneDim>::cabinet()->add(x);
        }
        catch (CanteraError) { return -1; }
    }


    int DLL_EXPORT onedim_del(int i) {
        Cabinet<OneDim>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT onedim_addFlow(int i, int n) {
        try {
            _onedim(i)->addDomain(_flow(n));
            return 0;
        }
        catch (CanteraError) { return -1; }
        //        catch (...) { return ERR; }
    }

//     int DLL_EXPORT onedim_addSurf(int i, int n) {
//         try {
//             _onedim(i)->addDomain(_surface(n));
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//     }

    int DLL_EXPORT onedim_eval(int i, doublereal* x0, doublereal* r) {
        try {
            _onedim(i)->eval(-1, x0, r, 0.0);
            return 0;
        }
        catch (CanteraError) { return -1; }
        //    catch (...) { return ERR; }
    }

    int DLL_EXPORT onedim_solve(int i, doublereal* x0, doublereal* x1, 
        int loglevel) {
        try {
            int m = _onedim(i)->solve(x0, x1, loglevel);
            return m;
        }
        catch (CanteraError) { return -1; }
        //catch (...) { return ERR; }
    }

    double DLL_EXPORT onedim_ssnorm(int i, doublereal* x0, doublereal* x1) {
        return _onedim(i)->ssnorm(x0, x1);
    }

    int DLL_EXPORT onedim_setsteadymode(int i) {
        if (_onedim(i)->transient()) {
            _onedim(i)->setSteadyMode();
            //_onedim(i)->jacobian().setAge(10000);
            return 1;
        }
        return 0;
    }

    int DLL_EXPORT onedim_settransientmode(int i, doublereal dt, doublereal* x) {
        _onedim(i)->initTimeInteg(dt, x);
        double rr = fabs(_onedim(i)->rdt()*dt - 1.0); 
        if ((rr > 1.e-5) || _onedim(i)->steady()) {
            //_onedim(i)->jacobian().setAge(10000);
            return 1;
        }
        return 0;
    }

    int DLL_EXPORT onedim_setnewtonoptions(int i, int maxage) {
        _onedim(i)->newton().setOptions(maxage);
        return 0;
    }

    int DLL_EXPORT onedim_resize(int i) {
        _onedim(i)->resize();
        return 0;
    }

    int DLL_EXPORT onedim_writeStats(int i) {
        _onedim(i)->writeStats();
        return 0;
    }

    double DLL_EXPORT onedim_timestep(int i, int nsteps, doublereal dt,
        doublereal* x, doublereal* xnew, int loglevel) {
        try {
            return _onedim(i)->timeStep(nsteps, dt, x, xnew, loglevel);
        }
        catch (CanteraError) { return -1.0; }
    }

    int DLL_EXPORT onedim_save(int i, char* fname, char* id, 
        char* desc, doublereal* soln) {
        try {
            _onedim(i)->save(string(fname), string(id), string(desc), soln);
            return 0;
        }
        catch (CanteraError) { return -1; }
        //catch (...) { return ERR; }
    }


}
