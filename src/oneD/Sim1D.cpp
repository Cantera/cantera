/**
 * @file Sim1D.cpp
 */

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/MultiNewton.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/xml.h"

#include <fstream>

using namespace std;

namespace Cantera
{

Sim1D::Sim1D(vector<Domain1D*>& domains) :
    OneDim(domains)
{
    // resize the internal solution vector and the work array, and perform
    // domain-specific initialization of the solution vector.
    resize();
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._getInitialSoln(&m_x[start(n)]);
    }

    // set some defaults
    m_tstep = 1.0e-5;
    m_steps = { 10 };
}

void Sim1D::setInitialGuess(const std::string& component, vector_fp& locs, vector_fp& vals)
{
    for (size_t dom=0; dom<nDomains(); dom++) {
        Domain1D& d = domain(dom);
        size_t ncomp = d.nComponents();
        for (size_t comp=0; comp<ncomp; comp++) {
            if (d.componentName(comp)==component) {
                setProfile(dom,comp,locs,vals);
            }
        }
    }
}

void Sim1D::setValue(size_t dom, size_t comp, size_t localPoint, doublereal value)
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::setValue",
                   "Index out of bounds: {} > {}", iloc, m_x.size());
    m_x[iloc] = value;
}

doublereal Sim1D::value(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::value",
                   "Index out of bounds: {} > {}", iloc, m_x.size());
    return m_x[iloc];
}

doublereal Sim1D::workValue(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::workValue",
                   "Index out of bounds: {} > {}", iloc, m_x.size());
    return m_xnew[iloc];
}

void Sim1D::setProfile(size_t dom, size_t comp,
                       const vector_fp& pos, const vector_fp& values)
{
    Domain1D& d = domain(dom);
    doublereal z0 = d.zmin();
    doublereal z1 = d.zmax();
    doublereal zpt, frac, v;
    for (size_t n = 0; n < d.nPoints(); n++) {
        zpt = d.z(n);
        frac = (zpt - z0)/(z1 - z0);
        v = linearInterp(frac, pos, values);
        setValue(dom, comp, n, v);
    }
}

void Sim1D::save(const std::string& fname, const std::string& id,
                 const std::string& desc, int loglevel)
{
    OneDim::save(fname, id, desc, m_x.data(), loglevel);
}

void Sim1D::saveResidual(const std::string& fname, const std::string& id,
                         const std::string& desc, int loglevel)
{
    vector_fp res(m_x.size(), -999);
    OneDim::eval(npos, &m_x[0], &res[0], 0.0);
    OneDim::save(fname, id, desc, &res[0], loglevel);
}

void Sim1D::restore(const std::string& fname, const std::string& id,
                    int loglevel)
{
    ifstream s(fname);
    if (!s) {
        throw CanteraError("Sim1D::restore",
                           "could not open input file "+fname);
    }

    XML_Node root;
    root.build(s);
    s.close();

    XML_Node* f = root.findID(id);
    if (!f) {
        throw CanteraError("Sim1D::restore","No solution with id = "+id);
    }

    vector<XML_Node*> xd = f->getChildren("domain");
    if (xd.size() != nDomains()) {
        throw CanteraError("Sim1D::restore", "Solution does not contain the "
            " correct number of domains. Found {} expected {}.\n",
            xd.size(), nDomains());
    }
    for (size_t m = 0; m < nDomains(); m++) {
        Domain1D& dom = domain(m);
        if (loglevel > 0 && xd[m]->attrib("id") != dom.id()) {
            writelog("Warning: domain names do not match: '" +
                     (*xd[m])["id"] + + "' and '" + dom.id() + "'\n");
        }
        dom.resize(domain(m).nComponents(), intValue((*xd[m])["points"]));
    }
    resize();
    m_xlast_ts.clear();
    for (size_t m = 0; m < nDomains(); m++) {
        domain(m).restore(*xd[m], &m_x[domain(m).loc()], loglevel);
    }
    finalize();
}

void Sim1D::setFlatProfile(size_t dom, size_t comp, doublereal v)
{
    size_t np = domain(dom).nPoints();
    size_t n;
    for (n = 0; n < np; n++) {
        setValue(dom, comp, n, v);
    }
}

void Sim1D::showSolution(ostream& s)
{
    for (size_t n = 0; n < nDomains(); n++) {
        if (domain(n).domainType() != cEmptyType) {
            domain(n).showSolution_s(s, &m_x[start(n)]);
        }
    }
}

void Sim1D::showSolution()
{
    for (size_t n = 0; n < nDomains(); n++) {
        if (domain(n).domainType() != cEmptyType) {
            writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+domain(n).id()
                     +" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
            domain(n).showSolution(&m_x[start(n)]);
        }
    }
}

void Sim1D::restoreTimeSteppingSolution()
{
    if (m_xlast_ts.empty()) {
        throw CanteraError("Sim1D::restoreTimeSteppingSolution",
                           "No successful time steps taken on this grid.");
    }
    m_x = m_xlast_ts;
}

void Sim1D::restoreSteadySolution()
{
    if (m_xlast_ss.empty()) {
        throw CanteraError("Sim1D::restoreSteadySolution",
                           "No successful steady state solution");
    }
    m_x = m_xlast_ss;
    for (size_t n = 0; n < nDomains(); n++) {
        vector_fp& z = m_grid_last_ss[n];
        domain(n).setupGrid(z.size(), z.data());
    }
}

void Sim1D::getInitialSoln()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._getInitialSoln(&m_x[start(n)]);
    }
}

void Sim1D::finalize()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._finalize(&m_x[start(n)]);
    }
}

void Sim1D::setTimeStep(doublereal stepsize, size_t n, integer* tsteps)
{
    m_tstep = stepsize;
    m_steps.resize(n);
    for (size_t i = 0; i < n; i++) {
        m_steps[i] = tsteps[i];
    }
}

int Sim1D::newtonSolve(int loglevel)
{
    int m = OneDim::solve(m_x.data(), m_xnew.data(), loglevel);
    if (m >= 0) {
        m_x = m_xnew;
        return 0;
    } else if (m > -10) {
        return -1;
    } else {
        throw CanteraError("Sim1D::newtonSolve",
                           "ERROR: OneDim::solve returned m = {}", m);
    }
}

void Sim1D::solve(int loglevel, bool refine_grid)
{
    int new_points = 1;
    int nsteps;
    doublereal dt = m_tstep;
    m_nsteps = 0;
    int soln_number = -1;
    finalize();

    while (new_points > 0) {
        size_t istep = 0;
        nsteps = m_steps[istep];

        bool ok = false;
        if (loglevel > 0) {
            writeline('.', 78, true, true);
        }
        while (!ok) {
            // Attempt to solve the steady problem
            setSteadyMode();
            newton().setOptions(m_ss_jac_age);
            debuglog("Attempt Newton solution of steady-state problem...", loglevel);
            int status = newtonSolve(loglevel-1);

            if (status == 0) {
                if (loglevel > 0) {
                    writelog("    success.\n\n");
                    writelog("Problem solved on [");
                    for (size_t mm = 1; mm < nDomains(); mm+=2) {
                        writelog("{}", domain(mm).nPoints());
                        if (mm + 2 < nDomains()) {
                            writelog(", ");
                        }
                    }
                    writelog("] point grid(s).\n");
                }
                if (loglevel > 6) {
                    save("debug_sim1d.xml", "debug",
                         "After successful Newton solve");
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.xml", "residual",
                                 "After successful Newton solve");
                }
                ok = true;
                soln_number++;
            } else {
                debuglog("    failure. \n", loglevel);
                if (loglevel > 6) {
                    save("debug_sim1d.xml", "debug",
                         "After unsuccessful Newton solve");
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.xml", "residual",
                                 "After unsuccessful Newton solve");
                }
                if (loglevel > 0) {
                    writelog("Take {} timesteps   ", nsteps);
                }
                dt = timeStep(nsteps, dt, m_x.data(), m_xnew.data(),
                              loglevel-1);
                m_xlast_ts = m_x;
                if (loglevel > 6) {
                    save("debug_sim1d.xml", "debug", "After timestepping");
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.xml", "residual",
                                 "After timestepping");
                }

                if (loglevel == 1) {
                    writelog(" {:10.4g} {:10.4g}\n", dt,
                             log10(ssnorm(m_x.data(), m_xnew.data())));
                }
                istep++;
                if (istep >= m_steps.size()) {
                    nsteps = m_steps.back();
                } else {
                    nsteps = m_steps[istep];
                }
                dt = std::min(dt, m_tmax);
            }
        }
        if (loglevel > 0) {
            writeline('.', 78, true, true);
        }
        if (loglevel > 2) {
            showSolution();
        }

        if (refine_grid) {
            new_points = refine(loglevel);
            if (new_points) {
                // If the grid has changed, preemptively reduce the timestep
                // to avoid multiple successive failed time steps.
                dt = m_tstep;
            }
            if (new_points && loglevel > 6) {
                save("debug_sim1d.xml", "debug", "After regridding");
            }
            if (new_points && loglevel > 7) {
                saveResidual("debug_sim1d.xml", "residual",
                             "After regridding");
            }
            if (new_points < 0) {
                writelog("Maximum number of grid points reached.");
                new_points = 0;
            }
        } else {
            debuglog("grid refinement disabled.\n", loglevel);
            new_points = 0;
        }
    }
}

int Sim1D::refine(int loglevel)
{
    int ianalyze, np = 0;
    vector_fp znew, xnew;
    doublereal xmid, zmid;
    std::vector<size_t> dsize;

    m_xlast_ss = m_x;
    m_grid_last_ss.clear();

    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        Refiner& r = d.refiner();

        // Save the old grid corresponding to the converged solution
        m_grid_last_ss.push_back(d.grid());

        // determine where new points are needed
        ianalyze = r.analyze(d.grid().size(), d.grid().data(), &m_x[start(n)]);
        if (ianalyze < 0) {
            return ianalyze;
        }

        if (loglevel > 0) {
            r.show();
        }

        np += r.nNewPoints();
        size_t comp = d.nComponents();

        // loop over points in the current grid
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        for (size_t m = 0; m < npnow; m++) {
            if (r.keepPoint(m)) {
                // add the current grid point to the new grid
                znew.push_back(d.grid(m));

                // do the same for the solution at this point
                for (size_t i = 0; i < comp; i++) {
                    xnew.push_back(value(n, i, m));
                }

                // now check whether a new point is needed in the interval to
                // the right of point m, and if so, add entries to znew and xnew
                // for this new point
                if (r.newPointNeeded(m) && m + 1 < npnow) {
                    // add new point at midpoint
                    zmid = 0.5*(d.grid(m) + d.grid(m+1));
                    znew.push_back(zmid);
                    np++;

                    // for each component, linearly interpolate
                    // the solution to this point
                    for (size_t i = 0; i < comp; i++) {
                        xmid = 0.5*(value(n, i, m) + value(n, i, m+1));
                        xnew.push_back(xmid);
                    }
                }
            } else {
                if (loglevel > 0) {
                    writelog("refine: discarding point at {}\n", d.grid(m));
                }
            }
        }
        dsize.push_back(znew.size() - nstart);
    }

    // At this point, the new grid znew and the new solution vector xnew have
    // been constructed, but the domains themselves have not yet been modified.
    // Now update each domain with the new grid.

    size_t gridstart = 0, gridsize;
    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        gridsize = dsize[n];
        d.setupGrid(gridsize, &znew[gridstart]);
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    m_x = xnew;
    resize();
    finalize();
    return np;
}

int Sim1D::setFixedTemperature(doublereal t)
{
    int np = 0;
    vector_fp znew, xnew;
    doublereal xmid;
    doublereal zfixed,interp_factor;
    doublereal z1 = 0.0, z2 = 0.0, t1,t2;
    size_t n, m, i;
    size_t m1 = 0;
    std::vector<size_t> dsize;

    for (n = 0; n < nDomains(); n++) {
        bool addnewpt=false;
        Domain1D& d = domain(n);
        size_t comp = d.nComponents();

        // loop over points in the current grid to determine where new point is
        // needed.
        FreeFlame* d_free = dynamic_cast<FreeFlame*>(&domain(n));
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        if (d_free) {
            for (m = 0; m < npnow-1; m++) {
                if (value(n,2,m) == t) {
                    zfixed = d.grid(m);
                    d_free->m_zfixed = zfixed;
                    d_free->m_tfixed = t;
                    addnewpt = false;
                    break;
                } else if ((value(n,2,m)<t) && (value(n,2,m+1)>t)) {
                    z1 = d.grid(m);
                    m1 = m;
                    z2 = d.grid(m+1);
                    t1 = value(n,2,m);
                    t2 = value(n,2,m+1);

                    zfixed = (z1-z2)/(t1-t2)*(t-t2)+z2;
                    d_free->m_zfixed = zfixed;
                    d_free->m_tfixed = t;
                    addnewpt = true;
                    break;
                    //copy solution domain and push back values
                }
            }
        }

        for (m = 0; m < npnow; m++) {
            // add the current grid point to the new grid
            znew.push_back(d.grid(m));

            // do the same for the solution at this point
            for (i = 0; i < comp; i++) {
                xnew.push_back(value(n, i, m));
            }
            if (m==m1 && addnewpt) {
                //add new point at zfixed
                znew.push_back(zfixed);
                np++;
                interp_factor = (zfixed-z2) / (z1-z2);
                // for each component, linearly interpolate
                // the solution to this point
                for (i = 0; i < comp; i++) {
                    xmid = interp_factor*(value(n, i, m) - value(n, i, m+1)) + value(n,i,m+1);
                    xnew.push_back(xmid);
                }
            }
        }
        dsize.push_back(znew.size() - nstart);
    }

    // At this point, the new grid znew and the new solution vector xnew have
    // been constructed, but the domains themselves have not yet been modified.
    // Now update each domain with the new grid.
    size_t gridstart = 0, gridsize;
    for (n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        gridsize = dsize[n];
        d.setupGrid(gridsize, &znew[gridstart]);
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    m_x = xnew;

    resize();
    finalize();
    return np;
}

void Sim1D::setRefineCriteria(int dom, doublereal ratio,
                              doublereal slope, doublereal curve, doublereal prune)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setCriteria(ratio, slope, curve, prune);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setCriteria(ratio, slope, curve, prune);
        }
    }
}

void Sim1D::setGridMin(int dom, double gridmin)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setGridMin(gridmin);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setGridMin(gridmin);
        }
    }
}

void Sim1D::setMaxGridPoints(int dom, int npoints)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setMaxPoints(npoints);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setMaxPoints(npoints);
        }
    }
}

doublereal Sim1D::jacobian(int i, int j)
{
    return OneDim::jacobian().value(i,j);
}

void Sim1D::evalSSJacobian()
{
    OneDim::evalSSJacobian(m_x.data(), m_xnew.data());
}

void Sim1D::resize()
{
    OneDim::resize();
    m_x.resize(size(), 0.0);
    m_xnew.resize(size(), 0.0);
}

}
