
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "MultiJac.h"
#include "MultiNewton.h"
#include "OneDim.h"

#include "ctml.h"
using namespace ctml;

namespace Cantera {

    /**
     * Default constructor. Create an empty object.
     */
    OneDim::OneDim() 
      : m_jac(0), m_newt(0), 
	m_rdt(0.0), m_jac_ok(false),
	m_nd(0), m_bw(0), m_size(0),
	m_nflow(0), m_nbdry(0), m_init(false),
	m_ss_jac_age(10), m_ts_jac_age(20),
	m_nevals(0), m_evaltime(0.0)
    {
        m_newt = new MultiNewton(1);
        m_solve_time = 0.0;
    }


    /**
     * Construct a OneDim container for the domains pointed at by the
     * input vector of pointers.
    */ 
    OneDim::OneDim(vector<Resid1D*> domains) : 
        m_jac(0), m_newt(0), 
	m_rdt(0.0), m_jac_ok(false),
	m_nd(0), m_bw(0), m_size(0),
	m_nflow(0), m_nbdry(0), m_init(false),
	m_ss_jac_age(10), m_ts_jac_age(20),
	m_nevals(0), m_evaltime(0.0)
    {

        // create a Newton iterator, and add each domain.
        m_newt = new MultiNewton(1);
        int nd = domains.size();
        int i;
        for (i = 0; i < nd; i++) addDomain(domains[i]);
        init();
        resize();
    }

    void OneDim::addDomain(Resid1D* d) {

        // if not the first domain, link it to the last (rightmost)
        // current domain
        int n = m_dom.size();
        if (m_dom.size() > 0) m_dom.back()->append(d);
        if (2*(n/2) == n)
            m_connect.push_back(d);
        else
            m_bulk.push_back(d);

        // add it also to the global domain list, and set its 
        // container and position
        m_dom.push_back(d);
        d->setContainer(this, m_nd);
        m_nd++;
        resize();
    }


    OneDim::~OneDim() {
        delete m_jac;
        delete m_newt;
    }

    MultiJac& OneDim::jacobian() { return *m_jac; }
    MultiNewton& OneDim::newton() { return *m_newt; }

    void OneDim::writeStats() {
        saveStats();
        char buf[100];
        sprintf(buf,"\nStatistics:\n\n Grid   Functions   Time      Jacobians   Time \n");
        writelog(buf);
        int n = m_gridpts.size();
        for (int i = 0; i < n; i++) {
            sprintf(buf,"%5i   %5i    %9.4f    %5i    %9.4f \n", 
                m_gridpts[i], m_funcEvals[i], m_funcElapsed[i], 
                m_jacEvals[i], m_jacElapsed[i]);
            writelog(buf);
        }
    }

    void OneDim::saveStats() {
        if (m_jac) {
            int nev = m_jac->nEvals();
            if (nev > 0) {
                m_gridpts.push_back(m_pts);
                m_jacEvals.push_back(m_jac->nEvals());
                m_jacElapsed.push_back(m_jac->elapsedTime());
                m_funcEvals.push_back(m_nevals);
                m_nevals = 0;
                m_funcElapsed.push_back(m_evaltime);
                m_evaltime = 0.0;
            }
        }
    }

    void OneDim::resize() {
        int i;
        m_bw = 0;
        m_nvars.clear();
        m_loc.clear();

        int lc = 0;

        saveStats();

        m_pts = 0;
        for (i = 0; i < m_nd; i++) {
            Resid1D* d = m_dom[i];

            int np = d->nPoints();
            int nv = d->nComponents();
            for (int n = 0; n < np; n++) {
                m_nvars.push_back(nv);
                m_loc.push_back(lc);
                lc += nv;
                m_pts++;
            }

            // update the Jacobian bandwidth
            int bw1, bw2 = 0;

            // bandwidth of the local block
            bw1 = 2*d->nComponents() - 1;

            // bandwidth of the block coupling the first point of this
            // domain to the last point of the previous domain
            if (i > 0) {
                bw2 = d->nComponents() + m_dom[i-1]->nComponents() - 1;
            }
            if (bw1 > m_bw) m_bw = bw1;
            if (bw2 > m_bw) m_bw = bw2;

            m_size = d->loc() + d->size();
        }
        m_newt->resize(size());
        m_mask.resize(size());
        
        delete m_jac;
        m_jac = new MultiJac(*this);
        m_jac_ok = false;

        for (i = 0; i < m_nd; i++)
            m_dom[i]->setJac(m_jac);
    }


    int OneDim::solve(doublereal* x, doublereal* xnew, int loglevel) {
        clock_t t0 = clock();
        static int iok = 0, inotok = 0;
        if (!m_jac_ok) {
            eval(-1, x, xnew, 0.0, 0); // m_rdt);
            m_jac->eval(x, xnew, 0.0); // m_rdt);
            m_jac->updateTransient(m_rdt, m_mask.begin());
            m_jac_ok = true;
            inotok++;
        }
        else
            iok++;
        int m = m_newt->solve(x, xnew, *this, *m_jac, loglevel);
        clock_t t1 = clock();
        m_solve_time += (t1 - t0)/(1.0*CLOCKS_PER_SEC);
        return m;
    }

    /**
     * Return a pointer to the domain that contains component i of the
     * global solution vector. The domains are scanned right-to-left,
     * and the first one with starting location less or equal to i is
     * returned.
     *
     * 8/26/02 changed '<' to '<='  DGG
     *
     */
    Resid1D* OneDim::pointDomain(int i) {
        Resid1D* d = right();
        while (d) {
            if (d->loc() <= i) return d;
            d = d->left();
        }
        return 0;
    }
 

    /**
     * Evaluate the multi-domain residual function, and return the
     * result in array r.  
     */
    void OneDim::eval(int j, double* x, double* r, doublereal rdt, int count) {
        clock_t t0 = clock();
        fill(r, r+m_size, 0.0);
        fill(m_mask.begin(), m_mask.end(), 0);
        if (rdt < 0.0) rdt = m_rdt;

        vector<Resid1D*>::iterator d; 
        for (d = m_bulk.begin(); d != m_bulk.end(); ++d)
            (*d)->eval(j, x, r, m_mask.begin(), rdt);
        for (d = m_connect.begin(); d != m_connect.end(); ++d)
            (*d)->eval(j, x, r, m_mask.begin(), rdt);

        if (count) {
            clock_t t1 = clock();
            m_evaltime += double(t1 - t0)/CLOCKS_PER_SEC;
            m_nevals++;
        }
    }


    /** 
     * The 'infinity' (maximum magnitude) norm of the steady-state
     * residual. Used only for diagnostic output.
     */
    doublereal OneDim::ssnorm(doublereal* x, doublereal* r) {
        eval(-1, x, r, 0.0);
        doublereal ss = 0.0;
        for (int i = 0; i < m_size; i++) { 
            ss = fmaxx(fabs(r[i]),ss);
        }
        return ss;
    }


    /**
     * Prepare for time stepping with timestep dt. 
     */
    void OneDim::initTimeInteg(doublereal dt, doublereal* x) {
        doublereal rdt_old = m_rdt;
        m_rdt = 1.0/dt;
        if (fabs(rdt_old - m_rdt) > Tiny) {
            m_jac->updateTransient(m_rdt, m_mask.begin());
        }
        Resid1D* d = left();
        while (d) {
            d->initTimeInteg(dt, x);
            d = d->right();
        }
    }
    
    /** 
     * Prepare to solve the steady-state problem.  Set the reciprocal
     * of the time step to zero, and, if it was previously non-zero,
     * signal that a new Jacobian will be needed.
     */
    void OneDim::setSteadyMode() {
        m_rdt = 0.0;
        m_jac->updateTransient(m_rdt, m_mask.begin());
    }

    /**
     * Initialize all domains. On the first call, this methods calls
     * the init method of each domain, proceeding from left to right.
     * Subsequent calls do nothing.
     */
    void OneDim::init() {
        if (!m_init) {
            Resid1D* d = left();
            while (d) {
                d->init();
                d = d->right();
            }
        }
        m_init = true;
    }

    /**
     * Signal that the current Jacobian is no longer valid.
     */
    void Resid1D::needJacUpdate() { 
        if (m_container)
            m_container->jacobian().setAge(10000);
    }

    /**
     * Take time steps using Backward Euler.
     *
     *  nsteps   -- number of steps
     *  dt       -- initial step size
     *  loglevel -- controls amount of printed diagnostics
     */
     doublereal OneDim::timeStep(int nsteps, doublereal dt, doublereal* x, 
         doublereal* r, int loglevel) {

         newton().setOptions(m_ts_jac_age);
                
         if (loglevel > 0) {
             writelog("Begin time integration.\n\n");
             writelog(" step    size (s)    log10(ss) ");
             writelog("===============================");
         }
         int n = 0, m;
         doublereal ss;
         char str[80];
         while (n < nsteps) {
             if (loglevel > 0) {
                 ss = ssnorm(x, r);
                 sprintf(str, " %4d  %10.4g  %10.4g" , n,dt,log10(ss));
                 writelog(str);
             }
             initTimeInteg(dt,x);
             m = solve(x, r, loglevel-1);
             if (m >= 0) {
                 n += 1;
                 if (loglevel > 0) writelog("\n");
                 copy(r, r + m_size, x);
             }
             else {
                 if (loglevel > 0) writelog("...failure.");
                 dt *= 0.5;
                 if (dt < 1.e-14)
                     throw CanteraError("OneDim::timeStep",
                         "Time integration failed.");
             }
         }
         setSteadyMode();
         newton().setOptions(m_ss_jac_age);        
         return dt;
     }

    void OneDim::save(string fname, string id, string desc, doublereal* sol) {

        struct tm *newtime;
        time_t aclock;
        ::time( &aclock );                  /* Get time in seconds */
        newtime = localtime( &aclock );   /* Convert time to struct tm form */

        XML_Node root("doc");
        ifstream fin(fname.c_str());
        XML_Node* ct;
        if (fin) {
            root.build(fin);
            XML_Node* same_ID = root.findID(id);
            int jid = 1;
            string idnew = id;
            while (same_ID != 0) {
                idnew = id + "_" + int2str(jid);
                jid++;
                same_ID = root.findID(idnew);
            }
            id = idnew;
            fin.close();
            ct = &root.child("ctml");
        }
        else {
            ct = &root.addChild("ctml");
        }
        XML_Node& sim = (XML_Node&)ct->addChild("simulation");
        sim.addAttribute("id",id);
        addString(sim,"timestamp",asctime(newtime));
        if (desc != "") addString(sim,"description",desc);
        
        Resid1D* d = left();
        while (d) {
            d->save(sim, sol);
            d = d->right();
        }
        ofstream s(fname.c_str());
        if (!s) 
            throw CanteraError("save","could not open file "+fname);
        ct->write(s);
        s.close();
        writelog("Solution saved to file "+fname+" as solution '"+id+"'.\n");
    }
}
