// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <iostream>
#include <string>
#include <vector>

#include "ctmatutils.h"
#include "cantera/clib/ctonedim.h"

using namespace std;

void onedimmethods(int nlhs, mxArray* plhs[],
                   int nrhs, const mxArray* prhs[])
{
    double vv = 0.0;
    int job = getInt(prhs[2]);
    size_t n, m;
    double* dom_ids, *h;
    int indx = 0;
    char* nm;
    int dom;
    dom = getInt(prhs[1]);
    int idom, icomp, localPoint;
    if (job < 10) {
        int ph, kin, tr, itype;
        size_t sz, nd;

        switch (job) {
        case 1:
            // construct a new stagnation flow instance
            checkNArgs(7, nrhs);
            ph = getInt(prhs[3]);
            kin = getInt(prhs[4]);
            tr = getInt(prhs[5]);
            itype = getInt(prhs[6]);
            indx = stflow_new(ph, kin, tr, itype);
            break;
        case 2:
            // construct a new Inlet1D instance
            checkNArgs(3, nrhs);
            indx = inlet_new();
            break;
        case 3:
            // construct a new Surf1D instance
            checkNArgs(3, nrhs);
            indx = surf_new();
            break;
        case 4:
            // construct a new Symm1D instance
            checkNArgs(3, nrhs);
            indx = symm_new();
            break;
        case 5:
            // construct a new Outlet1D instance
            checkNArgs(3, nrhs);
            indx = outlet_new();
            break;
        case 6:
            // construct a new ReactingSurf1D instance
            checkNArgs(4, nrhs);
            indx = reactingsurf_new();
            reactingsurf_setkineticsmgr(indx, getInt(prhs[3]));
            break;
        case 8: {
            // construct a new Sim1D instance
            checkNArgs(5, nrhs);
            nd = getInt(prhs[3]);
            dom_ids = mxGetPr(prhs[4]);
            m = mxGetM(prhs[4]);
            n = mxGetN(prhs[4]);
            sz = (m == 1) ? n : m;
            if (sz != nd) {
                mexErrMsgTxt("wrong size for domain array");
            }

            std::vector<int> ptrs(sz);
            for (size_t k = 0; k < sz; k++) {
                ptrs[k] = int(dom_ids[k]);
            }
            indx = sim1D_new(sz, &ptrs[0]);
            break;
        }

        // construct a new OutletRes1D instance
        case -2:
            checkNArgs(3,nrhs);
            indx = outletres_new();
            break;

        default:
            mexErrMsgTxt("onedimmethods: unknown object type");
        }

        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        h = mxGetPr(plhs[0]);
        *h = double(indx);
        if (indx < 0) {
            reportError();
        }
        return;
    } else if (job < 40) {
        // methods
        int k;

        switch (job) {
        case 10:
            checkNArgs(3, nrhs);
            vv = domain_del(dom);
            break;
        case 11:
            checkNArgs(3, nrhs);
            vv = (double) domain_nComponents(dom);
            break;
        case 12:
            checkNArgs(3, nrhs);
            vv = domain_type(dom);
            break;
        case 13:
            checkNArgs(3, nrhs);
            vv = (double) domain_index(dom);
            if (vv >= 0.0) {
                vv += 1.0;
            }
            break;
        case 14:
            checkNArgs(3, nrhs);
            vv = (double) domain_nPoints(dom);
            break;
        case 15:
            checkNArgs(3, nrhs);
            vv = bdry_temperature(dom);
            break;
        case 16:
            checkNArgs(4, nrhs);
            k = getInt(prhs[3]);
            vv = bdry_massFraction(dom, k);
            break;
        case 17:
            checkNArgs(3, nrhs);
            vv = bdry_mdot(dom);
            break;
        case 18:
            checkNArgs(4, nrhs);
            nm = getString(prhs[3]);
            vv = (double) domain_componentIndex(dom, nm);
            if (vv >= 0.0) {
                vv += 1.0;
            }
            break;
        case 19:
            checkNArgs(4, nrhs);
            localPoint = getInt(prhs[3]) - 1;
            vv = domain_grid(dom, localPoint);
            break;
        case 30:
            checkNArgs(6, nrhs);
            idom = getInt(prhs[3]) - 1;
            icomp = getInt(prhs[4]) - 1;
            localPoint = getInt(prhs[5]) - 1;
            vv = sim1D_value(dom, idom, icomp, localPoint);
            break;
        case 31:
            checkNArgs(6, nrhs);
            idom = getInt(prhs[3]) - 1;
            icomp = getInt(prhs[4]) - 1;
            localPoint = getInt(prhs[5]) - 1;
            vv = sim1D_workValue(dom, idom, icomp, localPoint);
            break;
        default:
            mexErrMsgTxt("unknown job");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = vv;
        if ((job != 30) && (vv == -1.0)) {
            reportError();
        }
        return;
    } else if (job < 50) {
        int iok = -1;
        int buflen, icomp;
        char* output_buf;
        switch (job) {
        case 40:
            icomp = getInt(prhs[3]) - 1;
            buflen = domain_componentName(dom, icomp, 0, 0);
            if (buflen > 0) {
                output_buf = (char*) mxCalloc(buflen, sizeof(char));
                iok = domain_componentName(dom, icomp, buflen, output_buf);
            }
            break;
        default:
            iok = -1;
        }
        if (iok >= 0) {
            plhs[0] = mxCreateString(output_buf);
            return;
        } else {
            mexErrMsgTxt("error or unknown method.");
            return;
        }
    } else {
    // set parameters
        int iok = -1;
        double lower, upper, rtol, atol, *grid, *pos, *values,
               mdot, t, p, val, *temp, ratio, slope, curve, tstep, *dts,
               rdt, prune;
        size_t npts, np, nv;
        int comp, localPoint, idom,
            loglevel, refine_grid, n, flag, itime, ns, icount,
            onoff, ss_age, ts_age;
        char* xstr, *fname, *id, *desc, *name;
        switch (job) {
        case 51:
            checkNArgs(6, nrhs);
            n = getInt(prhs[3]) - 1;
            lower = getDouble(prhs[4]);
            upper = getDouble(prhs[5]);
            iok = domain_setBounds(dom, n, lower, upper);
            break;
        case 53:
            checkNArgs(4, nrhs);
            grid = mxGetPr(prhs[3]);
            npts = mxGetM(prhs[3]) * mxGetN(prhs[3]);
            iok = domain_setupGrid(dom, npts, grid);
            break;
        case 54:
            id = getString(prhs[3]);
            iok = domain_setID(dom, id);
            break;
        case 55:
        case 56:
            checkNArgs(6, nrhs);
            n = getInt(prhs[3]) - 1;
            rtol = getDouble(prhs[4]);
            atol = getDouble(prhs[5]);
            if (job == 55) {
                iok = domain_setSteadyTolerances(dom, n, rtol, atol);
            } else {
                iok = domain_setTransientTolerances(dom, n, rtol, atol);
            }
            break;
        case 60:
            checkNArgs(4, nrhs);
            mdot = getDouble(prhs[3]);
            iok = bdry_setMdot(dom, mdot);
            break;
        case 61:
            checkNArgs(4, nrhs);
            t = getDouble(prhs[3]);
            iok = bdry_setTemperature(dom, t);
            break;
        case 62:
            checkNArgs(4, nrhs);
            xstr = getString(prhs[3]);
            iok = bdry_setMoleFractions(dom, xstr);
            break;
        case 63:
            checkNArgs(4, nrhs);
            p = getDouble(prhs[3]);
            iok = stflow_setPressure(dom, p);
            break;
        case 64:
            checkNArgs(5, nrhs);
            pos = mxGetPr(prhs[3]);
            temp = mxGetPr(prhs[4]);
            np = mxGetM(prhs[3])*mxGetN(prhs[3]);
            nv = mxGetM(prhs[4])*mxGetN(prhs[4]);
            iok = stflow_setFixedTempProfile(dom, np, pos, nv, temp);
            break;
        case 66:
            checkNArgs(4, nrhs);
            flag = getInt(prhs[3]);
            iok = stflow_solveEnergyEqn(dom, flag);
            break;
        case 100:
            checkNArgs(7, nrhs);
            idom = getInt(prhs[3]) - 1;
            comp = getInt(prhs[4]) - 1;
            localPoint = getInt(prhs[5]) -1;
            val = getDouble(prhs[6]);
            iok = sim1D_setValue(dom, idom, comp, localPoint, val);
            break;
        case 101:
            checkNArgs(7, nrhs);
            idom = getInt(prhs[3]) - 1;
            comp = getInt(prhs[4]) - 1;
            pos = mxGetPr(prhs[5]);
            values = mxGetPr(prhs[6]);
            np = mxGetM(prhs[5])*mxGetN(prhs[5]);
            nv = mxGetM(prhs[6])*mxGetN(prhs[6]);
            iok = sim1D_setProfile(dom, idom, comp, np, pos, nv, values);
            break;
        case 102:
            checkNArgs(6, nrhs);
            idom = getInt(prhs[3]) - 1;
            comp = getInt(prhs[4]) - 1;
            val = getDouble(prhs[5]);
            iok = sim1D_setFlatProfile(dom, idom, comp, val);
            break;
        case 103:
            checkNArgs(4, nrhs);
            fname = getString(prhs[3]);
            iok = sim1D_showSolution(dom, fname);
            break;
        case 104:
            checkNArgs(5, nrhs);
            loglevel = getInt(prhs[3]);
            refine_grid = getInt(prhs[4]);
            iok = sim1D_solve(dom, loglevel, refine_grid);
            break;
        case 105:
            checkNArgs(4, nrhs);
            loglevel = getInt(prhs[3]);
            iok = sim1D_refine(dom, loglevel);
            break;
        case 106:
            checkNArgs(8, nrhs);
            idom = getInt(prhs[3]) - 1;
            ratio = getDouble(prhs[4]);
            slope = getDouble(prhs[5]);
            curve = getDouble(prhs[6]);
            prune = getDouble(prhs[7]);
            iok = sim1D_setRefineCriteria(dom, idom, ratio, slope, curve, prune);
            break;
        case 107:
            iok = 0;
            checkNArgs(6, nrhs);
            fname = getString(prhs[3]);
            id = getString(prhs[4]);
            desc = getString(prhs[5]);
            iok = sim1D_save(dom, fname, id, desc);
            break;
        case 108:
            checkNArgs(3, nrhs);
            iok = sim1D_writeStats(dom, 1);
            break;
        case 109:
            checkNArgs(4, nrhs);
            name = getString(prhs[3]);
            iok = sim1D_domainIndex(dom, name);
            if (iok >= 0) {
                iok++;
            }
            break;
        case 110:
            checkNArgs(3, nrhs);
            iok = sim1D_del(dom);
            break;
        case 111:
            iok = 0;
            checkNArgs(5, nrhs);
            fname = getString(prhs[3]);
            id = getString(prhs[4]);
            iok = sim1D_restore(dom, fname, id);
            break;
        case 112: {
            tstep = getDouble(prhs[3]);
            ns = getInt(prhs[4]);
            dts = mxGetPr(prhs[5]);
            std::vector<int> nsteps(ns);
            for (n = 0; n < ns; n++) {
                nsteps[n] = int(dts[n]);
            }
            iok = sim1D_setTimeStep(dom, tstep, ns, &nsteps[0]);
            break;
        }
        case 113:
            checkNArgs(5, nrhs);
            rdt = getDouble(prhs[3]);
            icount = getInt(prhs[4]);
            iok = sim1D_eval(dom, rdt, icount);
            break;
        case 114:
            checkNArgs(5, nrhs);
            ss_age = getInt(prhs[3]);
            ts_age = getInt(prhs[4]);
            iok = sim1D_setMaxJacAge(dom, ss_age, ts_age);
            break;
        case 120:
            checkNArgs(4, nrhs);
            onoff = getInt(prhs[3]);
            iok = reactingsurf_enableCoverageEqs(dom, onoff);
            break;
        default:
            mexPrintf(" job = %d ",job);
            mexErrMsgTxt("unknown parameter");
        }
        if (iok < 0) {
            reportError();
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(iok);
        return;
    }
}
