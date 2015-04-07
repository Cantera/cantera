/**
 * @file Sim1D.cpp
 */

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/MultiJac.h"

#include <fstream>

using namespace std;

namespace Cantera
{

static void sim1D_drawline()
{
    string s(78,'.');
    s += '\n';
    writelog(s.c_str());
}

Sim1D::Sim1D() :
    OneDim()
{
    //writelog("Sim1D default constructor\n");
}

Sim1D::Sim1D(vector<Domain1D*>& domains) :
    OneDim(domains)
{
    // resize the internal solution vector and the work array, and perform
    // domain-specific initialization of the solution vector.

    m_x.resize(size(), 0.0);
    m_xnew.resize(size(), 0.0);
    for (size_t n = 0; n < m_nd; n++) {
        domain(n)._getInitialSoln(DATA_PTR(m_x) + start(n));
        domain(n).m_adiabatic=false;
    }

    // set some defaults
    m_tstep = 1.0e-5;
    //m_maxtimestep = 10.0;
    m_steps.push_back(1);
    m_steps.push_back(2);
    m_steps.push_back(5);
    m_steps.push_back(10);
}

void Sim1D::setInitialGuess(const std::string& component, vector_fp& locs, vector_fp& vals)
{
    for (size_t dom=0; dom<m_nd; dom++) {
        Domain1D& d = domain(dom);
        size_t ncomp = d.nComponents();
        for (size_t comp=0; comp<ncomp; comp++) {
            if (d.componentName(comp)==component) {
                setProfile(dom,comp,locs,vals);
            }
        }
    }
}

void Sim1D::setValue(size_t dom, size_t comp, size_t localPoint,  doublereal value)
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::setValue",
                   "Index out of bounds:" + int2str(iloc) + " > " +
                   int2str(m_x.size()));
    m_x[iloc] = value;
}

doublereal Sim1D::value(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::value",
                   "Index out of bounds:" + int2str(iloc) + " > " +
                   int2str(m_x.size()));
    return m_x[iloc];
}

doublereal Sim1D::workValue(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_x.size(), "Sim1D::workValue",
                   "Index out of bounds:" + int2str(iloc) + " > " +
                   int2str(m_x.size()));
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
    OneDim::save(fname, id, desc, DATA_PTR(m_x), loglevel);
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
    ifstream s(fname.c_str());
    if (!s)
        throw CanteraError("Sim1D::restore",
                           "could not open input file "+fname);

    XML_Node root;
    root.build(s);
    s.close();

    XML_Node* f = root.findID(id);
    if (!f) {
        throw CanteraError("Sim1D::restore","No solution with id = "+id);
    }

    vector<XML_Node*> xd;
    f->getChildren("domain", xd);
    if (xd.size() != m_nd) {
        throw CanteraError("Sim1D::restore", "Solution does not contain the "
                           " correct number of domains. Found " +
                           int2str(xd.size()) + "expected " +
                           int2str(m_nd) + ".\n");
    }
    size_t sz = 0;
    for (size_t m = 0; m < m_nd; m++) {
        if (loglevel > 0 && xd[m]->attrib("id") != domain(m).id()) {
            writelog("Warning: domain names do not match: '" +
                     (*xd[m])["id"] + + "' and '" + domain(m).id() + "'\n");
        }
        sz += domain(m).nComponents() * intValue((*xd[m])["points"]);
    }
    m_x.resize(sz);
    m_xnew.resize(sz);
    for (size_t m = 0; m < m_nd; m++) {
        domain(m).restore(*xd[m], DATA_PTR(m_x) + domain(m).loc(), loglevel);
    }
    resize();
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
    for (size_t n = 0; n < m_nd; n++) {
        if (domain(n).domainType() != cEmptyType) {
            domain(n).showSolution_s(s, DATA_PTR(m_x) + start(n));
        }
    }
}

void Sim1D::showSolution()
{
    for (size_t n = 0; n < m_nd; n++) {
        if (domain(n).domainType() != cEmptyType) {
            writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+domain(n).id()
                     +" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
            domain(n).showSolution(DATA_PTR(m_x) + start(n));
        }
    }
}

void Sim1D::getInitialSoln()
{
    for (size_t n = 0; n < m_nd; n++) {
        domain(n)._getInitialSoln(DATA_PTR(m_x) + start(n));
    }
}

void Sim1D::finalize()
{
    for (size_t n = 0; n < m_nd; n++) {
        domain(n)._finalize(DATA_PTR(m_x) + start(n));
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
    int m = OneDim::solve(DATA_PTR(m_x), DATA_PTR(m_xnew), loglevel);
    if (m >= 0) {
        copy(m_xnew.begin(), m_xnew.end(), m_x.begin());
        return 0;
    } else if (m > -10) {
        return -1;
    } else {
        throw CanteraError("Sim1D::newtonSolve",
                           "ERROR: OneDim::solve returned m = " + int2str(m) + "\n");
    }
}

void Sim1D::solve(int loglevel, bool refine_grid)
{
    int new_points = 1;
    int nsteps;
    doublereal dt = m_tstep;
    int soln_number = -1;
    finalize();

    while (new_points > 0) {
        size_t istep = 0;
        nsteps = m_steps[istep];

        bool ok = false;
        if (loglevel > 0) {
            writelog("\n");
            sim1D_drawline();
        }
        while (!ok) {
            writelog("Attempt Newton solution of steady-state problem...", loglevel);
            int status = newtonSolve(loglevel-1);

            if (status == 0) {
                if (loglevel > 0) {
                    writelog("    success.\n\n");
                    writelog("Problem solved on [");
                    for (size_t mm = 1; mm < nDomains(); mm+=2) {
                        writelog(int2str(domain(mm).nPoints()));
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
                char buf[100];
                writelog("    failure. \n", loglevel);
                if (loglevel > 6) {
                    save("debug_sim1d.xml", "debug",
                         "After unsuccessful Newton solve");
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.xml", "residual",
                                 "After unsuccessful Newton solve");
                }
                writelog("Take "+int2str(nsteps)+" timesteps   ", loglevel);
                dt = timeStep(nsteps, dt, DATA_PTR(m_x), DATA_PTR(m_xnew),
                              loglevel-1);
                if (loglevel > 6) {
                    save("debug_sim1d.xml", "debug", "After timestepping");
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.xml", "residual",
                                 "After timestepping");
                }

                if (loglevel == 1) {
                    sprintf(buf, " %10.4g %10.4g \n", dt,
                            log10(ssnorm(DATA_PTR(m_x), DATA_PTR(m_xnew))));
                    writelog(buf);
                }
                istep++;
                if (istep >= m_steps.size()) {
                    nsteps = m_steps.back();
                } else {
                    nsteps = m_steps[istep];
                }
                if (dt > m_tmax) {
                    dt = m_tmax;
                }
            }
        }
        if (loglevel > 0) {
            sim1D_drawline();
            writelog("\n");
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
            writelog("grid refinement disabled.\n", loglevel);
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

    for (size_t n = 0; n < m_nd; n++) {
        Domain1D& d = domain(n);
        Refiner& r = d.refiner();

        // determine where new points are needed
        ianalyze = r.analyze(d.grid().size(),
                             DATA_PTR(d.grid()), DATA_PTR(m_x) + start(n));
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

                // now check whether a new point is needed in the
                // interval to the right of point m, and if so, add
                // entries to znew and xnew for this new point

                if (r.newPointNeeded(m) && m + 1 < npnow) {

                    // add new point at midpoint
                    zmid = 0.5*(d.grid(m) + d.grid(m+1));
                    znew.push_back(zmid);
                    np++;
                    //writelog(string("refine: adding point at ")+fp2str(zmid)+"\n");

                    // for each component, linearly interpolate
                    // the solution to this point
                    for (size_t i = 0; i < comp; i++) {
                        xmid = 0.5*(value(n, i, m) + value(n, i, m+1));
                        xnew.push_back(xmid);
                    }
                }
            } else {
                writelog("refine: discarding point at "+fp2str(d.grid(m))+"\n", loglevel);
            }
        }
        dsize.push_back(znew.size() - nstart);
    }

    // At this point, the new grid znew and the new solution
    // vector xnew have been constructed, but the domains
    // themselves have not yet been modified.  Now update each
    // domain with the new grid.

    size_t gridstart = 0, gridsize;
    for (size_t n = 0; n < m_nd; n++) {
        Domain1D& d = domain(n);
        //            Refiner& r = d.refiner();
        gridsize = dsize[n]; // d.nPoints() + r.nNewPoints();
        d.setupGrid(gridsize, DATA_PTR(znew) + gridstart);
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    m_x.resize(xnew.size());
    copy(xnew.begin(), xnew.end(), m_x.begin());

    // resize the work array
    m_xnew.resize(xnew.size());

    //        copy(xnew.begin(), xnew.end(), m_xnew.begin());

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


    for (n = 0; n < m_nd; n++) {
        bool addnewpt=false;
        Domain1D& d = domain(n);

        size_t comp = d.nComponents();

        // loop over points in the current grid to determine where new point is needed.
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        for (m = 0; m < npnow-1; m++) {
            if (value(n,2,m) == t) {
                zfixed = d.grid(m);
                //set d.zfixed, d.ztemp
                d.m_zfixed = zfixed;
                d.m_tfixed = t;
                addnewpt = false;
                break;
            } else if ((value(n,2,m)<t) && (value(n,2,m+1)>t)) {
                z1 = d.grid(m);
                m1 = m;
                z2 = d.grid(m+1);
                t1 = value(n,2,m);
                t2 = value(n,2,m+1);

                zfixed = (z1-z2)/(t1-t2)*(t-t2)+z2;
                //set d.zfixed, d.ztemp;
                d.m_zfixed = zfixed;
                d.m_tfixed = t;
                addnewpt = true;
                break;
                //copy solution domain and push back values
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

    // At this point, the new grid znew and the new solution
    // vector xnew have been constructed, but the domains
    // themselves have not yet been modified.  Now update each
    // domain with the new grid.

    size_t gridstart = 0, gridsize;
    for (n = 0; n < m_nd; n++) {
        Domain1D& d = domain(n);
        //            Refiner& r = d.refiner();
        gridsize = dsize[n]; // d.nPoints() + r.nNewPoints();
        d.setupGrid(gridsize, DATA_PTR(znew) + gridstart);
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    m_x.resize(xnew.size());
    copy(xnew.begin(), xnew.end(), m_x.begin());

    // resize the work array
    m_xnew.resize(xnew.size());

    copy(xnew.begin(), xnew.end(), m_xnew.begin());

    resize();
    finalize();
    return np;
}

void Sim1D::setAdiabaticFlame(void)
{
    for (size_t n = 0; n < m_nd; n++) {
        Domain1D& d = domain(n);
        d.m_adiabatic=true;
    }
}

void Sim1D::setRefineCriteria(int dom, doublereal ratio,
                              doublereal slope, doublereal curve, doublereal prune)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setCriteria(ratio, slope, curve, prune);
    } else {
        for (size_t n = 0; n < m_nd; n++) {
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
        for (size_t n = 0; n < m_nd; n++) {
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
        for (size_t n = 0; n < m_nd; n++) {
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
    OneDim::evalSSJacobian(DATA_PTR(m_x), DATA_PTR(m_xnew));
}
}
