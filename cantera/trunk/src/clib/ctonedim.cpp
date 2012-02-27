/**
 * @file ctonedim.cpp
 */
#define CANTERA_USE_INTERNAL

#include "ctonedim.h"

// Cantera includes
#include "cantera/base/config.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/numerics/DenseMatrix.h"

// local includes
#include "Cabinet.h"
#include "Storage.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<Sim1D> SimCabinet;
typedef Cabinet<Domain1D> DomainCabinet;
template<> SimCabinet* SimCabinet::__storage = 0;
template<> DomainCabinet* DomainCabinet::__storage = 0;

static StFlow* _stflow(int i)
{
    Domain1D* d = &DomainCabinet::item(i);
    if (d->domainType() == cFlowType) {
        return dynamic_cast<StFlow*>(d);
    } else {
        throw CanteraError("_stflow","wrong domain type");
    }
    return 0;
}

static Bdry1D* _bdry(int i)
{
    Domain1D* d = &DomainCabinet::item(i);
    if (! d->isConnector()) {
        throw CanteraError("_bdry","wrong domain type: " +int2str(d->domainType()));
    }
    return dynamic_cast<Bdry1D*>(d);
}

inline ThermoPhase* _phase(int n)
{
    return Storage::__storage->__thtable[n];
}

inline Kinetics* _kinetics(int n)
{
    return Storage::__storage->__ktable[n];
}

inline ThermoPhase* _thermo(int n)
{
    return Storage::__storage->__thtable[n];
}

inline Transport* _transport(int n)
{
    return Storage::__storage->__trtable[n];
}


extern "C" {

    int DLL_EXPORT domain_clear()
    {
        try {
            DomainCabinet::clear();
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT domain_del(int i)
    {
        DomainCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT domain_type(int i)
    {
        return DomainCabinet::item(i).domainType();
    }

    size_t DLL_EXPORT domain_index(int i)
    {
        return DomainCabinet::item(i).domainIndex();
    }

    size_t DLL_EXPORT domain_nComponents(int i)
    {
        return DomainCabinet::item(i).nComponents();
    }

    size_t DLL_EXPORT domain_nPoints(int i)
    {
        return DomainCabinet::item(i).nPoints();
    }

    int DLL_EXPORT domain_componentName(int i, int n, int sz, char* buf)
    {
        try {
            string nm = DomainCabinet::item(i).componentName(n);
            size_t lout = std::min<size_t>(sz, nm.size());
            copy(nm.c_str(), nm.c_str() + lout, buf);
            buf[lout] = '\0';
            return static_cast<int>(nm.size());
        } catch (CanteraError) {
            return -1;
        }
    }

    size_t DLL_EXPORT domain_componentIndex(int i, char* name)
    {
        try {
            size_t n = DomainCabinet::item(i).componentIndex(string(name));
            return n;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT domain_grid(int i, int n)
    {
        try {
            return DomainCabinet::item(i).grid(n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    int DLL_EXPORT domain_setBounds(int i, int n, double lower, double upper)
    {
        try {
            DomainCabinet::item(i).setBounds(n, lower, upper);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT domain_upperBound(int i, int n)
    {
        try {
            return DomainCabinet::item(i).upperBound(n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    double DLL_EXPORT domain_lowerBound(int i, int n)
    {
        try {
            return DomainCabinet::item(i).lowerBound(n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    int DLL_EXPORT domain_setTolerances(int i, int n, double rtol,
                                        double atol, int itime)
    {
        try {
            DomainCabinet::item(i).setTolerances(n, rtol, atol, itime);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT domain_rtol(int i, int n)
    {
        try {
            return DomainCabinet::item(i).rtol(n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    double DLL_EXPORT domain_atol(int i, int n)
    {
        try {
            return DomainCabinet::item(i).atol(n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    int DLL_EXPORT domain_setupGrid(int i, size_t npts, double* grid)
    {
        try {
            DomainCabinet::item(i).setupGrid(npts, grid);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT domain_setID(int i, char* id)
    {
        try {
            string s = string(id);
            DomainCabinet::item(i).setID(s);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }


    int DLL_EXPORT domain_setDesc(int i, char* desc)
    {
        try {
            string s = string(desc);
            DomainCabinet::item(i).setDesc(s);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }


    int DLL_EXPORT inlet_new()
    {
        try {
            Inlet1D* i = new Inlet1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT surf_new()
    {
        try {
            Surf1D* i = new Surf1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT reactingsurf_new()
    {
        try {
            //writelog("in reactingsurf_new\n");
            Domain1D* i = new ReactingSurf1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            writelog("error");
            return -1;
        }
    }

    int DLL_EXPORT symm_new()
    {
        try {
            Symm1D* i = new Symm1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT outlet_new()
    {
        try {
            Outlet1D* i = new Outlet1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT outletres_new()
    {
        try {
            OutletRes1D* i = new OutletRes1D();
            return DomainCabinet::add(i);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT bdry_setMdot(int i, double mdot)
    {
        try {
            _bdry(i)->setMdot(mdot);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT bdry_setTemperature(int i, double t)
    {
        try {
            _bdry(i)->setTemperature(t);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT bdry_setMoleFractions(int i, char* x)
    {
        try {
            _bdry(i)->setMoleFractions(string(x));
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT bdry_temperature(int i)
    {
        try {
            return _bdry(i)->temperature();
        } catch (CanteraError) {
            return DERR;
        }
    }

    double DLL_EXPORT bdry_massFraction(int i, int k)
    {
        try {
            return _bdry(i)->massFraction(k);
        } catch (CanteraError) {
            return DERR;
        }
    }

    double DLL_EXPORT bdry_mdot(int i)
    {
        try {
            return _bdry(i)->mdot();
        } catch (CanteraError) {
            return DERR;
        }
    }

    int DLL_EXPORT reactingsurf_setkineticsmgr(int i, int j)
    {
        try {
            ReactingSurf1D* srf = (ReactingSurf1D*)_bdry(i);
            InterfaceKinetics* k = (InterfaceKinetics*)_kinetics(j);
            srf->setKineticsMgr(k);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT reactingsurf_enableCoverageEqs(int i, int onoff)
    {
        try {
            ReactingSurf1D* srf = (ReactingSurf1D*)_bdry(i);
            srf->enableCoverageEquations(onoff != 0);
            return 0;
        }

        catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT inlet_setSpreadRate(int i, double v)
    {
        try {
            Inlet1D* inlt = (Inlet1D*)_bdry(i);
            inlt->setSpreadRate(v);
            return 0;
        }

        catch (CanteraError) {
            return -1;
        }
    }

    //------------------ stagnation flow domains --------------------

    int DLL_EXPORT stflow_new(int iph, int ikin, int itr, int itype)
    {
        try {
            IdealGasPhase* ph = (IdealGasPhase*)_thermo(iph);
            if (itype == 1) {
                AxiStagnFlow* x = new AxiStagnFlow(ph, ph->nSpecies(), 2);
                x->setKinetics(*_kinetics(ikin));
                x->setTransport(*_transport(itr));
                return DomainCabinet::add(x);
            } else if (itype == 2) {
                FreeFlame* x = new FreeFlame(ph, ph->nSpecies(), 2);
                x->setKinetics(*_kinetics(ikin));
                x->setTransport(*_transport(itr));
                return DomainCabinet::add(x);
            } else {
                return -2;
            }
        } catch (CanteraError) {
            return -1;
        }
    }


    int DLL_EXPORT stflow_setTransport(int i, int itr, int iSoret)
    {
        bool withSoret = false;
        if (iSoret > 0) {
            withSoret = true;
        }
        try {
            _stflow(i)->setTransport(*_transport(itr), withSoret);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT stflow_enableSoret(int i, int iSoret)
    {
        bool withSoret = false;
        if (iSoret > 0) {
            withSoret = true;
        }
        try {
            _stflow(i)->enableSoret(withSoret);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT stflow_setPressure(int i, double p)
    {
        try {
            _stflow(i)->setPressure(p);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT stflow_setFixedTempProfile(int i, size_t n, double* pos,
            size_t m, double* temp)
    {
        try {
            vector_fp vpos(n), vtemp(n);
            for (size_t j = 0; j < n; j++) {
                vpos[j] = pos[j];
                vtemp[j] = temp[j];
            }
            _stflow(i)->setFixedTempProfile(vpos, vtemp);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }


    int DLL_EXPORT stflow_solveSpeciesEqs(int i, int flag)
    {
        try {
            if (flag > 0) {
                _stflow(i)->solveSpecies(-1);
            } else {
                _stflow(i)->fixSpecies(-1);
            }
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }


    int DLL_EXPORT stflow_solveEnergyEqn(int i, int flag)
    {
        try {
            if (flag > 0) {
                _stflow(i)->solveEnergyEqn(-1);
            } else {
                _stflow(i)->fixTemperature(-1);
            }
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }


    //------------------- Sim1D --------------------------------------

    int DLL_EXPORT sim1D_new(size_t nd, int* domains)
    {
        vector<Domain1D*> d;
        try {
            //  cout << "nd = " << nd << endl;
            for (size_t n = 0; n < nd; n++) {
                //writelog("n = "+int2str(n)+"\n");
                //writelog("dom = "+int2str(domains[n])+"\n");
                d.push_back(&DomainCabinet::item(domains[n]));
            }
            //writelog("in sim1D_new, calling new Sim1D\n");
            Sim1D* s = new Sim1D(d);
            //writelog("in sim1D_new, ret Sim1D\n");
            return SimCabinet::add(s);
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_clear()
    {
        try {
            SimCabinet::clear();
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_del(int i)
    {
        SimCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT sim1D_setValue(int i, int dom, int comp,
                                  int localPoint, double value)
    {
        try {
            SimCabinet::item(i).setValue(dom, comp, localPoint, value);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setProfile(int i, int dom, int comp,
                                    size_t np, double* pos, size_t nv, double* v)
    {
        try {
            vector_fp vv, pv;
            for (size_t n = 0; n < np; n++) {
                vv.push_back(v[n]);
                pv.push_back(pos[n]);
            }
            SimCabinet::item(i).setProfile(dom, comp, pv, vv);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setFlatProfile(int i, int dom, int comp, double v)
    {
        try {
            SimCabinet::item(i).setFlatProfile(dom, comp, v);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_showSolution(int i, char* fname)
    {
        string fn = string(fname);
        if (fn == "-") {
            SimCabinet::item(i).showSolution();
        } else {
            ofstream fout(fname);
            SimCabinet::item(i).showSolution(fout);
            fout.close();
        }
        return 0;
    }

    int DLL_EXPORT sim1D_setTimeStep(int i, double stepsize, size_t ns, integer* nsteps)
    {
        try {
            SimCabinet::item(i).setTimeStep(stepsize, ns, nsteps);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_getInitialSoln(int i)
    {
        try {
            SimCabinet::item(i).getInitialSoln();
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_solve(int i, int loglevel, int refine_grid)
    {
        try {
            bool r = (refine_grid == 0 ? false : true);
            SimCabinet::item(i).solve(loglevel, r);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_refine(int i, int loglevel)
    {
        try {
            SimCabinet::item(i).refine(loglevel);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setRefineCriteria(int i, int dom, double ratio,
                                           double slope, double curve, double prune)
    {
        try {
            SimCabinet::item(i).setRefineCriteria(dom, ratio, slope, curve, prune);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_save(int i, char* fname, char* id,
                              char* desc)
    {
        try {
            string sname = string(fname);
            string sid = string(id);
            string sdesc = string(desc);
            SimCabinet::item(i).save(sname, sid, sdesc);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_restore(int i, char* fname, char* id)
    {
        try {
            string sname = string(fname);
            string sid = string(id);
            SimCabinet::item(i).restore(sname, sid);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_writeStats(int i, int printTime)
    {
        try {
            SimCabinet::item(i).writeStats(printTime);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_domainIndex(int i, char* name)
    {
        try {
            return (int) SimCabinet::item(i).domainIndex(string(name));
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT sim1D_value(int i, int idom, int icomp, int localPoint)
    {
        try {
            return SimCabinet::item(i).value(idom, icomp, localPoint);
        } catch (CanteraError) {
            return DERR;
        }
    }

    double DLL_EXPORT sim1D_workValue(int i, int idom, int icomp, int localPoint)
    {
        try {
            return SimCabinet::item(i).workValue(idom, icomp, localPoint);
        } catch (CanteraError) {
            return DERR;
        }
    }

    int DLL_EXPORT sim1D_eval(int i, double rdt, int count)
    {
        try {
            SimCabinet::item(i).eval(rdt, count);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setMaxJacAge(int i, int ss_age, int ts_age)
    {
        try {
            SimCabinet::item(i).setJacAge(ss_age, ts_age);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_timeStepFactor(int i, double tfactor)
    {
        try {
            SimCabinet::item(i).setTimeStepFactor(tfactor);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setTimeStepLimits(int i, double tsmin, double tsmax)
    {
        try {
            if (tsmin > 0.0) {
                SimCabinet::item(i).setMinTimeStep(tsmin);
            }
            if (tsmax > 0.0) {
                SimCabinet::item(i).setMaxTimeStep(tsmax);
            }
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_setFixedTemperature(int i, double temp)
    {
        try {
            SimCabinet::item(i).setFixedTemperature(temp);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT sim1D_evalSSJacobian(int i)
    {
        try {
            SimCabinet::item(i).evalSSJacobian();
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT sim1D_jacobian(int i, int m, int n)
    {
        try {
            return SimCabinet::item(i).jacobian(m,n);
        } catch (CanteraError) {
            return DERR;
        }
    }

    size_t DLL_EXPORT sim1D_size(int i)
    {
        try {
            return SimCabinet::item(i).size();
        } catch (CanteraError) {
            return -1;
        }
    }
}
