/**
 * @file Sim1D.cpp
 */

#include "Sim1D.h"

namespace Cantera {


    static void drawline() {
        string s(78,'.');
        s += '\n';
        writelog(s.c_str());
    }

    Sim1D::Sim1D() : OneDim() { 
      //writelog("Sim1D default constructor\n"); 
    }

    Sim1D::Sim1D(vector<Domain1D*>& domains) : OneDim(domains) {

        // resize the internal solution vector and the wprk array,
        // and perform domain-specific initialization of the
        // solution vector.
      //writelog("size = "+int2str(size())+"\n");
        m_x.resize(size(), 0.0);
        m_xnew.resize(size(), 0.0);
        for (int n = 0; n < m_nd; n++) {
	  //    writelog("calling domain "+int2str(n)+" _getInitialSoln\n");
            domain(n)._getInitialSoln(m_x.begin() + start(n));
	    // writelog("ret\n");
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

    doublereal Sim1D::workValue(int dom, int comp, int localPoint) const {
        int iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
        return m_xnew[iloc];
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

        Domain1D& d = domain(dom);
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


    void Sim1D::save(string fname, string id, string desc) {
        OneDim::save(fname, id, desc, m_x.begin());
    }

    /**
     * Initialize the solution with a previously-saved solution.
     */
    void Sim1D::restore(string fname, string id) {
        ifstream s(fname.c_str());
        //char buf[100];
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
        int sz = 0, np, nv, m;
        for (m = 0; m < m_nd; m++) {
            XML_Node* d = f->findID(domain(m).id());
            if (!d) {
                writelog("No data for domain "+domain(m).id());
                xd.push_back(0);
                sz += domain(m).nComponents();
            }
            else {
                const XML_Node& node = *d;
                xd.push_back(d);
                np = intValue(node["points"]);
                nv = intValue(node["components"]);
                sz += np*domain(m).nComponents();
            }
        } 
        m_x.resize(sz);
        m_xnew.resize(sz);
        for (m = 0; m < m_nd; m++) {
            if (xd[m]) {
                domain(m).restore(*xd[m], m_x.begin() + domain(m).loc());
            }
        }
        resize();
        finalize();
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

    void Sim1D::showSolution() {
        for (int n = 0; n < m_nd; n++) {
            writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+domain(n).id()
                       +" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
            domain(n).showSolution(m_x.begin() + start(n));
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
            writelog(string("ERROR: solve returned m = ") + int2str(m) + "\n");
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
                        drawline();
                        writelog("\nAttempt Newton solution of steady-state problem...");
                    }
                    newtonSolve(loglevel-1);
                    
                    if (loglevel > 0) {
                        writelog("    success.\n\n");
                        writelog("Problem solved on [");
                        for (int mm = 1; mm < nDomains(); mm+=2) {
                            writelog(int2str(domain(mm).nPoints()));
                            if (mm < nDomains() - 2) writelog(", ");
                        }
                        writelog("]");
                        writelog(" point grid(s).\n\n");
                    }
                    ok = true;
                    soln_number++;
                    
                }

                catch (CanteraError) {

                    char buf[100];
                    if (loglevel > 0) {
                        writelog("    failure. \n\n");
                        drawline();
                        //                    }
                        //if (loglevel == 1) 
                        writelog("Take "+int2str(nsteps)+
                            " timesteps   ");
                    }
                    dt = timeStep(nsteps, dt, m_x.begin(), m_xnew.begin(), 
                        loglevel-1);
                    if (loglevel == 1) {
                        sprintf(buf, " %10.4g %10.4g \n", dt, 
                            log10(ssnorm(m_x.begin(), m_xnew.begin())));
                        writelog(buf);
                    }
                    istep++;
                    if (istep >= int(m_steps.size())) {
                        nsteps = m_steps.back();
                        //           dt *= 2.0;
                    }
                    else {
                        nsteps = m_steps[istep];
                    }
                    if (dt > m_tmax) dt = m_tmax;
                }
            }
            if (loglevel > 2) showSolution();
            
            if (refine_grid) {
                //writelog("calling refine.\n");
                new_points = refine(loglevel);
            }
            else {
                writelog("grid refinement disabled.\n");
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
        vector_int dsize;

        for (n = 0; n < m_nd; n++) {
            strt = znew.size();
            Domain1D& d = domain(n);
            Refiner& r = d.refiner();

            // determine where new points are needed
            r.analyze(d.grid().size(), d.grid().begin(), m_x.begin() + start(n));
            if (loglevel > 0) { r.show(); }

            np += r.nNewPoints();
            int comp = d.nComponents();

            // loop over points in the current grid
            int npnow = d.nPoints();
            int nstart = znew.size();
            for (m = 0; m < npnow; m++) {

                if (r.keepPoint(m)) {
                    // add the current grid point to the new grid
                    znew.push_back(d.grid(m));
                    
                    // do the same for the solution at this point
                    for (i = 0; i < comp; i++) {
                        xnew.push_back(value(n, i, m));
                    }

                    // now check whether a new point is needed in the
                    // interval to the right of point m, and if so, add
                    // entries to znew and xnew for this new point

                    if (r.newPointNeeded(m) && m < npnow - 1) {

                        // add new point at midpoint
                        zmid = 0.5*(d.grid(m) + d.grid(m+1));
                        znew.push_back(zmid);
                        np++;
                        //writelog(string("refine: adding point at ")+fp2str(zmid)+"\n");

                        // for each component, linearly interpolate
                        // the solution to this point
                        for (i = 0; i < comp; i++) {
                            xmid = 0.5*(value(n, i, m) + value(n, i, m+1));
                            xnew.push_back(xmid);
                        }
                    }
                }
                else {
                    writelog(string("refine: discarding point at ")+fp2str(d.grid(m))+"\n");
                    ; // throw CanteraError("refine","keepPoint is false at m = "+int2str(m));
                }
            }
            dsize.push_back(znew.size() - nstart);
        }

        // At this point, the new grid znew and the new solution
        // vector xnew have been constructed, but the domains
        // themselves have not yet been modified.  Now update each
        // domain with the new grid.

        int gridstart = 0, gridsize;
        for (n = 0; n < m_nd; n++) {
            Domain1D& d = domain(n);
            //            Refiner& r = d.refiner();
            gridsize = dsize[n]; // d.nPoints() + r.nNewPoints();
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
        doublereal slope, doublereal curve, doublereal prune) {
        if (dom >= 0) {
            Refiner& r = domain(dom).refiner();
            r.setCriteria(ratio, slope, curve, prune);
        }
        else {
            for (int n = 0; n < m_nd; n++) {
                Refiner& r = domain(n).refiner();
                r.setCriteria(ratio, slope, curve, prune);
            }                    
        }
    }
            
}
