/**
 * @file ctonedim.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CANTERA_USE_INTERNAL

#include "cantera/clib/ctonedim.h"

// Cantera includes
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/transport/TransportBase.h"
#include "Cabinet.h"

#include <fstream>

using namespace std;
using namespace Cantera;

typedef Cabinet<Sim1D> SimCabinet;
typedef Cabinet<Domain1D> DomainCabinet;
template<> SimCabinet* SimCabinet::s_storage = 0;
template<> DomainCabinet* DomainCabinet::s_storage = 0;

typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Transport> TransportCabinet;

extern "C" {

    int ct_clearOneDim()
    {
        try {
            DomainCabinet::clear();
            SimCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int domain_del(int i)
    {
        try {
            DomainCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1 , ERR);
        }
    }

    int domain_type(int i)
    {
        try {
            return DomainCabinet::item(i).domainType();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t domain_index(int i)
    {
        try {
            return DomainCabinet::item(i).domainIndex();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t domain_nComponents(int i)
    {
        try {
            return DomainCabinet::item(i).nComponents();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t domain_nPoints(int i)
    {
        try {
            return DomainCabinet::item(i).nPoints();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int domain_componentName(int i, int n, int sz, char* buf)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            return static_cast<int>(copyString(dom.componentName(n), buf, sz));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t domain_componentIndex(int i, const char* name)
    {
        try {
            return DomainCabinet::item(i).componentIndex(name);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double domain_grid(int i, int n)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkPointIndex(n);
            return dom.grid(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int domain_setBounds(int i, int n, double lower, double upper)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            dom.setBounds(n, lower, upper);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double domain_upperBound(int i, int n)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            return dom.upperBound(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double domain_lowerBound(int i, int n)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            return dom.lowerBound(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int domain_setSteadyTolerances(int i, int n, double rtol,
                                   double atol)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            dom.setSteadyTolerances(rtol, atol, n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int domain_setTransientTolerances(int i, int n, double rtol,
                                      double atol)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            dom.setTransientTolerances(rtol, atol, n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double domain_rtol(int i, int n)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            return dom.rtol(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double domain_atol(int i, int n)
    {
        try {
            Domain1D& dom = DomainCabinet::item(i);
            dom.checkComponentIndex(n);
            return dom.atol(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int domain_setupGrid(int i, size_t npts, const double* grid)
    {
        try {
            DomainCabinet::item(i).setupGrid(npts, grid);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int domain_setID(int i, const char* id)
    {
        try {
            DomainCabinet::item(i).setID(id);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int inlet_new()
    {
        try {
            return DomainCabinet::add(new Inlet1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_new()
    {
        try {
            return DomainCabinet::add(new Surf1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactingsurf_new()
    {
        try {
            return DomainCabinet::add(new ReactingSurf1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int symm_new()
    {
        try {
            return DomainCabinet::add(new Symm1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int outlet_new()
    {
        try {
            return DomainCabinet::add(new Outlet1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int outletres_new()
    {
        try {
            return DomainCabinet::add(new OutletRes1D());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int bdry_setMdot(int i, double mdot)
    {
        try {
            DomainCabinet::get<Boundary1D>(i).setMdot(mdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int bdry_setTemperature(int i, double t)
    {
        try {
            DomainCabinet::get<Boundary1D>(i).setTemperature(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int bdry_setMoleFractions(int i, const char* x)
    {
        try {
            DomainCabinet::get<Boundary1D>(i).setMoleFractions(x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double bdry_temperature(int i)
    {
        try {
            return DomainCabinet::get<Boundary1D>(i).temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double bdry_massFraction(int i, int k)
    {
        try {
            return DomainCabinet::get<Boundary1D>(i).massFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double bdry_mdot(int i)
    {
        try {
            return DomainCabinet::get<Boundary1D>(i).mdot();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int reactingsurf_setkineticsmgr(int i, int j)
    {
        try {
            InterfaceKinetics& k = Cabinet<Kinetics>::get<InterfaceKinetics>(j);
            DomainCabinet::get<ReactingSurf1D>(i).setKineticsMgr(&k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactingsurf_enableCoverageEqs(int i, int onoff)
    {
        try {
            DomainCabinet::get<ReactingSurf1D>(i).enableCoverageEquations(onoff != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int inlet_setSpreadRate(int i, double v)
    {
        try {
            DomainCabinet::get<Inlet1D>(i).setSpreadRate(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //------------------ stagnation flow domains --------------------

    int stflow_new(int iph, int ikin, int itr, int itype)
    {
        try {
            IdealGasPhase& ph = ThermoCabinet::get<IdealGasPhase>(iph);
            StFlow* x = new StFlow(&ph, ph.nSpecies(), 2);
            if (itype == 1) {
                x->setAxisymmetricFlow();
            } else if (itype == 2) {
                x->setFreeFlow();
            } else {
                delete x;
                return -2;
            }
            x->setKinetics(KineticsCabinet::item(ikin));
            x->setTransport(TransportCabinet::item(itr));
            return DomainCabinet::add(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    int stflow_setTransport(int i, int itr)
    {
        try {
            DomainCabinet::get<StFlow>(i).setTransport(TransportCabinet::item(itr));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int stflow_enableSoret(int i, int iSoret)
    {
        try {
            bool withSoret = (iSoret > 0);
            DomainCabinet::get<StFlow>(i).enableSoret(withSoret);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int stflow_setPressure(int i, double p)
    {
        try {
            DomainCabinet::get<StFlow>(i).setPressure(p);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double stflow_pressure(int i)
    {
        try {
            return DomainCabinet::get<StFlow>(i).pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int stflow_setFixedTempProfile(int i, size_t n, const double* pos,
                                   size_t m, const double* temp)
    {
        try {
            vector_fp vpos(n), vtemp(n);
            for (size_t j = 0; j < n; j++) {
                vpos[j] = pos[j];
                vtemp[j] = temp[j];
            }
            DomainCabinet::get<StFlow>(i).setFixedTempProfile(vpos, vtemp);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int stflow_solveEnergyEqn(int i, int flag)
    {
        try {
            if (flag > 0) {
                DomainCabinet::get<StFlow>(i).solveEnergyEqn(npos);
            } else {
                DomainCabinet::get<StFlow>(i).fixTemperature(npos);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //------------------- Sim1D --------------------------------------

    int sim1D_new(size_t nd, const int* domains)
    {
        try {
            vector<Domain1D*> d;
            for (size_t n = 0; n < nd; n++) {
                d.push_back(&DomainCabinet::item(domains[n]));
            }
            return SimCabinet::add(new Sim1D(d));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_del(int i)
    {
        try {
            SimCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setValue(int i, int dom, int comp, int localPoint, double value)
    {
        try {
            SimCabinet::item(i).setValue(dom, comp, localPoint, value);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setProfile(int i, int dom, int comp, size_t np, const double* pos,
                         size_t nv, const double* v)
    {
        try {
            Sim1D& sim = SimCabinet::item(i);
            sim.checkDomainIndex(dom);
            sim.domain(dom).checkComponentIndex(comp);
            vector_fp vv, pv;
            for (size_t n = 0; n < np; n++) {
                vv.push_back(v[n]);
                pv.push_back(pos[n]);
            }
            sim.setProfile(dom, comp, pv, vv);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setFlatProfile(int i, int dom, int comp, double v)
    {
        try {
            Sim1D& sim = SimCabinet::item(i);
            sim.checkDomainIndex(dom);
            sim.domain(dom).checkComponentIndex(comp);
            sim.setFlatProfile(dom, comp, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_showSolution(int i, const char* fname)
    {
        try {
            string fn = string(fname);
            if (fn == "-") {
                SimCabinet::item(i).showSolution();
            } else {
                ofstream fout(fname);
                SimCabinet::item(i).showSolution(fout);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setTimeStep(int i, double stepsize, size_t ns, const int* nsteps)
    {
        try {
            SimCabinet::item(i).setTimeStep(stepsize, ns, nsteps);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_getInitialSoln(int i)
    {
        try {
            SimCabinet::item(i).getInitialSoln();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_solve(int i, int loglevel, int refine_grid)
    {
        try {
            bool r = (refine_grid == 0 ? false : true);
            SimCabinet::item(i).solve(loglevel, r);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_refine(int i, int loglevel)
    {
        try {
            SimCabinet::item(i).refine(loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setRefineCriteria(int i, int dom, double ratio,
                                double slope, double curve, double prune)
    {
        try {
            SimCabinet::item(i).setRefineCriteria(dom, ratio, slope, curve, prune);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setGridMin(int i, int dom, double gridmin)
    {
        try {
            SimCabinet::item(i).setGridMin(dom, gridmin);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_save(int i, const char* fname, const char* id, const char* desc)
    {
        try {
            SimCabinet::item(i).save(fname, id, desc);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_restore(int i, const char* fname, const char* id)
    {
        try {
            SimCabinet::item(i).restore(fname, id);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_writeStats(int i, int printTime)
    {
        try {
            SimCabinet::item(i).writeStats(printTime);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_domainIndex(int i, const char* name)
    {
        try {
            return (int) SimCabinet::item(i).domainIndex(name);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double sim1D_value(int i, int idom, int icomp, int localPoint)
    {
        try {
            Sim1D& sim = SimCabinet::item(i);
            sim.checkDomainIndex(idom);
            sim.domain(idom).checkComponentIndex(icomp);
            return sim.value(idom, icomp, localPoint);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double sim1D_workValue(int i, int idom, int icomp, int localPoint)
    {
        try {
            Sim1D& sim = SimCabinet::item(i);
            sim.checkDomainIndex(idom);
            sim.domain(idom).checkComponentIndex(icomp);
            return sim.workValue(idom, icomp, localPoint);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int sim1D_eval(int i, double rdt, int count)
    {
        try {
            SimCabinet::item(i).eval(rdt, count);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setMaxJacAge(int i, int ss_age, int ts_age)
    {
        try {
            SimCabinet::item(i).setJacAge(ss_age, ts_age);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int sim1D_setFixedTemperature(int i, double temp)
    {
        try {
            SimCabinet::item(i).setFixedTemperature(temp);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
