/**
 * @file Sim1D.cpp
 */

#include "Sim1D.h"

namespace Cantera {


    Sim1D::Sim1D(vector<Resid1D*>& domains) : OneDim(domains) {

        // resize the internal solution vector and the wprk array,
        // and perform domain-specific initialization of the
        // solution vector.
        m_x.resize(size(), 0.0);
        m_xnew.resize(size(), 0.0);
        for (int n = 0; n < m_nd; n++) {
            domain(n)._getInitialSoln(m_x.begin() + start(n));
        }

        // set some defaults
        m_tstep = 1.0e-5;
        //m_maxtimestep = 10.0;
        m_steps.push_back(1);
        m_steps.push_back(2);
        m_steps.push_back(5);
        m_steps.push_back(10);

    }


    /**
     * Set a single value in the solution vector.
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param localPoint grid point within the domain, beginning with 0 for
     * the leftmost grid point in the domain.
     * @param value the value.
     */
    void Sim1D::setValue(int dom, int comp, int localPoint,  doublereal value) {
        int iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
        m_x[iloc] = value;
    }


    /**
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param localPoint grid point within the domain, beginning with 0 for
     * the leftmost grid point in the domain.
     */
    doublereal Sim1D::value(int dom, int comp, int localPoint) const {
        int iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
        return m_x[iloc];
    }


    /**
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param pos A vector of relative positions, beginning with 0.0 at the
     * left of the domain, and ending with 1.0 at the right of the domain.
     * @param values A vector of values corresponding to the relative position
     * locations. 
     *
     * Note that the vector pos and values can have lengths
     * different than the number of grid points, but their lengths
     * must be equal. The values at the grid points will be
     * linearly interpolated based on the (pos, values)
     * specification.
     */
    void Sim1D::setProfile(int dom, int comp, 
        const vector_fp& pos, const vector_fp& values) {

        Resid1D& d = domain(dom);
        int np = d.nPoints();
        int n;
        doublereal z0 = d.zmin();
        doublereal z1 = d.zmax();
        doublereal zpt, frac, v;
        for (n = 0; n < np; n++) {
            zpt = d.z(n);
            frac = (zpt - z0)/(z1 - z0);
            v = linearInterp(frac, pos, values);
            setValue(dom, comp, n, v);
        }
    }



    void Sim1D::setFlatProfile(int dom, int comp, doublereal v) {
        int np = domain(dom).nPoints();
        int n;
        for (n = 0; n < np; n++) { setValue(dom, comp, n, v); }
    }


    void Sim1D::showSolution(ostream& s) {
        for (int n = 0; n < m_nd; n++) {
            domain(n).showSolution(s, m_x.begin() + start(n));
        }
    }


    void Sim1D::finalize() {
        for (int n = 0; n < m_nd; n++) {
            domain(n)._finalize(m_x.begin() + start(n));
        }
    }


    void Sim1D::setTimeStep(doublereal stepsize, int n, integer* tsteps) {
        m_tstep = stepsize;
        m_steps.resize(n);
        for (int i = 0; i < n; i++) m_steps[i] = tsteps[i];
    }


    void Sim1D::newtonSolve(int loglevel) {
        int m = OneDim::solve(m_x.begin(), m_xnew.begin(), loglevel);
        if (m >= 0) 
            copy(m_xnew.begin(), m_xnew.end(), m_x.begin());
        else if (m > -10)
            throw CanteraError("Sim1D::newtonSolve","no solution found");
        else {
            cout << "ERROR: solve returned m = " << m << endl;
            exit(-1);
        }
    }


    void Sim1D::solve(int loglevel, bool refine_grid) {

        int new_points = 1;
        int istep, nsteps;
        doublereal dt = m_tstep;
        int soln_number = -1;

        finalize();

        while (new_points > 0) {

            istep = 0;
            nsteps = m_steps[istep];

            bool ok = false;
            while (!ok) {
                    
                try {
                    if (loglevel > 0) {
                        writelog("Attempt Newton solution of steady-state problem...");
                    }
                    newtonSolve(loglevel-1);
                    
                    if (loglevel > 0) {
                        writelog("success.\n\n");
                        //writelog("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n"); 
                        writelog("Problem solved on [");
                        for (int mm = 1; mm < nDomains(); mm+=2) {
                            writelog(int2str(domain(mm).nPoints()));
                            if (mm < nDomains() - 2) writelog(", ");
                        }
                        writelog("]");
                        writelog(" point grid(s).\n\n");
                        //writelog("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"); 
                    }
                    ok = true;
                    soln_number++;
                    
                }

                catch (CanteraError) {

                    char buf[100];
                    if (loglevel > 0) writelog("failure. \n\n");
                    if (loglevel == 1) writelog("Take "+int2str(nsteps)+" timesteps   ");
                    dt = timeStep(nsteps, dt, m_x.begin(), m_xnew.begin(), loglevel-1);
                    if (loglevel == 1) {
                        sprintf(buf, " %10.4g %10.4g \n", dt, 
                            log10(ssnorm(m_x.begin(), m_xnew.begin())));
                        writelog(buf);
                    }
                    istep++;
                    if (istep >= int(m_steps.size())) {
                        nsteps = m_steps.back();
                        dt *= 2.0;
                        cout << " doubled dt = " << dt << endl;
                    }
                    else {
                        nsteps = m_steps[istep];
                    }
                    if (dt > m_tmax) dt = m_tmax;
                }
            }
            if (loglevel > 2) showSolution(cout);
            
            if (refine_grid) {
                new_points = refine(loglevel);
            }
            else {
                new_points = 0;
            }
        }
    }

    /**
     * Refine the grid in all domains.
     */
    int Sim1D::refine(int loglevel) {
        int np = 0;
        vector_fp znew, xnew;
        doublereal xmid, zmid;
        int strt, n, m, i;

        for (n = 0; n < m_nd; n++) {
            strt = znew.size();
            Resid1D& d = domain(n);
            Refiner& r = d.refiner();

            // determine where new points are needed
            r.analyze(d.grid().size(), d.grid().begin(), m_x.begin() + start(n));
            if (loglevel > 0) { r.show(); }

            np += r.nNewPoints();
            int comp = d.nComponents();

            // loop over points in the current grid
            int npnow = d.nPoints();
            for (m = 0; m < npnow; m++) {

                // add the current grid point to the new grid
                znew.push_back(d.grid(m));

                // do the same for the solution at this point
                for (i = 0; i < comp; i++) {
                    xnew.push_back(value(n, i, m));
                }

                // now check whether a new point is needed in the interval to the
                // right of point m, and if so, add entries to znew and xnew for
                // this new point

                if (r.newPointNeeded(m)) {

                    // add new point at midpoint
                    zmid = 0.5*(d.grid(m) + d.grid(m+1));
                    znew.push_back(zmid);

                    // for each component, linearly interpolate the solution to 
                    // this point
                    for (i = 0; i < comp; i++) {
                        xmid = 0.5*(value(n, i, m) + value(n, i, m+1));
                        xnew.push_back(xmid);
                    }
                }
            }
        }

        // At this point, the new grid znew and the new solution vector xnew have
        // been constructed, but the domains themselves have not yet been modified.
        // Now update each domain with the new grid.

        int gridstart = 0, gridsize;
        for (n = 0; n < m_nd; n++) {
            Resid1D& d = domain(n);
            Refiner& r = d.refiner();
            gridsize = d.nPoints() + r.nNewPoints();
            d.setupGrid(gridsize, znew.begin() + gridstart);
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


    /**
     * Set grid refinement criteria. If dom >= 0, then the settings
     * apply only to the specified domain.  If dom < 0, the settings
     * are applied to each domain.  @see Refiner::setCriteria.
     */
    void Sim1D::setRefineCriteria(int dom, doublereal ratio,
        doublereal slope, doublereal curve) {
        if (dom >= 0) {
            Refiner& r = domain(dom).refiner();
            r.setCriteria(ratio, slope, curve);
        }
        else {
            for (int n = 0; n < m_nd; n++) {
                Refiner& r = domain(n).refiner();
                r.setCriteria(ratio, slope, curve);
            }                    
        }
    }
            
}
