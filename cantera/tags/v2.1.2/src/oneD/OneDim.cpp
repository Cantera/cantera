//! @file OneDim.cpp
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/MultiNewton.h"
#include "cantera/oneD/OneDim.h"

#include "cantera/numerics/Func1.h"
#include "cantera/base/ctml.h"

#include <fstream>

using namespace ctml;
using namespace std;

namespace Cantera
{

OneDim::OneDim()
    : m_tmin(1.0e-16), m_tmax(10.0), m_tfactor(0.5),
      m_jac(0), m_newt(0),
      m_rdt(0.0), m_jac_ok(false),
      m_nd(0), m_bw(0), m_size(0),
      m_init(false),
      m_ss_jac_age(10), m_ts_jac_age(20),
      m_interrupt(0), m_nevals(0), m_evaltime(0.0)
{
    //writelog("OneDim default constructor\n");
    m_newt = new MultiNewton(1);
    //m_solve_time = 0.0;
}

OneDim::OneDim(vector<Domain1D*> domains) :
    m_tmin(1.0e-16), m_tmax(10.0), m_tfactor(0.5),
    m_jac(0), m_newt(0),
    m_rdt(0.0), m_jac_ok(false),
    m_nd(0), m_bw(0), m_size(0),
    m_init(false),
    m_ss_jac_age(10), m_ts_jac_age(20),
    m_interrupt(0), m_nevals(0), m_evaltime(0.0)
{
    //writelog("OneDim constructor\n");

    // create a Newton iterator, and add each domain.
    m_newt = new MultiNewton(1);
    int nd = static_cast<int>(domains.size());
    int i;
    for (i = 0; i < nd; i++) {
        addDomain(domains[i]);
    }
    init();
    resize();
}

size_t OneDim::domainIndex(const std::string& name)
{
    for (size_t n = 0; n < m_nd; n++) {
        if (domain(n).id() == name) {
            return n;
        }
    }
    throw CanteraError("OneDim::domainIndex","no domain named >>"+name+"<<");
    return npos;
}

void OneDim::addDomain(Domain1D* d)
{
    // if 'd' is not the first domain, link it to the last domain
    // added (the rightmost one)
    int n = static_cast<int>(m_dom.size());
    if (n > 0) {
        m_dom.back()->append(d);
    }

    // every other domain is a connector
    if (2*(n/2) == n) {
        m_connect.push_back(d);
    } else {
        m_bulk.push_back(d);
    }

    // add it also to the global domain list, and set its
    // container and position
    m_dom.push_back(d);
    d->setContainer(this, m_nd);
    m_nd++;
    resize();
}

OneDim::~OneDim()
{
    delete m_jac;
    delete m_newt;
}

MultiJac& OneDim::jacobian()
{
    return *m_jac;
}
MultiNewton& OneDim::newton()
{
    return *m_newt;
}

void OneDim::writeStats(int printTime)
{
    saveStats();
    char buf[100];
    sprintf(buf,"\nStatistics:\n\n Grid   Functions   Time      Jacobians   Time \n");
    writelog(buf);
    size_t n = m_gridpts.size();
    for (size_t i = 0; i < n; i++) {
        if (printTime) {
            sprintf(buf,"%5s   %5i    %9.4f    %5i    %9.4f \n",
                    int2str(m_gridpts[i]).c_str(), m_funcEvals[i], m_funcElapsed[i],
                    m_jacEvals[i], m_jacElapsed[i]);
        } else {
            sprintf(buf,"%5s   %5i       NA        %5i        NA    \n",
                    int2str(m_gridpts[i]).c_str(), m_funcEvals[i], m_jacEvals[i]);
        }
        writelog(buf);
    }
}

void OneDim::saveStats()
{
    if (m_jac) {
        int nev = m_jac->nEvals();
        if (nev > 0 && m_nevals > 0) {
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

void OneDim::resize()
{
    m_bw = 0;
    std::vector<size_t> nvars, loc;
    size_t lc = 0;

    // save the statistics for the last grid
    saveStats();
    m_pts = 0;
    for (size_t i = 0; i < m_nd; i++) {
        Domain1D* d = m_dom[i];

        size_t np = d->nPoints();
        size_t nv = d->nComponents();
        for (size_t n = 0; n < np; n++) {
            nvars.push_back(nv);
            loc.push_back(lc);
            lc += nv;
            m_pts++;
        }

        // update the Jacobian bandwidth
        size_t bw1, bw2 = 0;

        // bandwidth of the local block
        bw1 = d->bandwidth();
        if (bw1 == npos) {
            bw1 = 2*d->nComponents() - 1;
        }

        // bandwidth of the block coupling the first point of this
        // domain to the last point of the previous domain
        if (i > 0) {
            bw2 = m_dom[i-1]->bandwidth();
            if (bw2 == npos) {
                bw2 = m_dom[i-1]->nComponents();
            }
            bw2 += d->nComponents() - 1;
        }
        if (bw1 > m_bw) {
            m_bw = bw1;
        }
        if (bw2 > m_bw) {
            m_bw = bw2;
        }

        m_size = d->loc() + d->size();
    }
    m_nvars = nvars;
    m_loc = loc;

    m_newt->resize(size());
    m_mask.resize(size());

    // delete the current Jacobian evaluator and create a new one
    delete m_jac;
    m_jac = new MultiJac(*this);
    m_jac_ok = false;

    for (size_t i = 0; i < m_nd; i++) {
        m_dom[i]->setJac(m_jac);
    }
}

int OneDim::solve(doublereal* x, doublereal* xnew, int loglevel)
{
    if (!m_jac_ok) {
        eval(npos, x, xnew, 0.0, 0);
        m_jac->eval(x, xnew, 0.0);
        m_jac->updateTransient(m_rdt, DATA_PTR(m_mask));
        m_jac_ok = true;
    }
    return m_newt->solve(x, xnew, *this, *m_jac, loglevel);
}

void OneDim::evalSSJacobian(doublereal* x, doublereal* xnew)
{
    doublereal rdt_save = m_rdt;
    m_jac_ok = false;
    setSteadyMode();
    eval(npos, x, xnew, 0.0, 0);
    m_jac->eval(x, xnew, 0.0);
    m_rdt = rdt_save;
}

Domain1D* OneDim::pointDomain(size_t i)
{
    Domain1D* d = right();
    while (d) {
        if (d->loc() <= i) {
            return d;
        }
        d = d->left();
    }
    return 0;
}

void OneDim::eval(size_t j, double* x, double* r, doublereal rdt, int count)
{
    clock_t t0 = clock();
    if (m_interrupt) {
        m_interrupt->eval(m_nevals);
    }
    fill(r, r + m_size, 0.0);
    fill(m_mask.begin(), m_mask.end(), 0);
    if (rdt < 0.0) {
        rdt = m_rdt;
    }
    //        int nn;
    vector<Domain1D*>::iterator d;

    // iterate over the bulk domains first
    for (d = m_bulk.begin(); d != m_bulk.end(); ++d) {
        (*d)->eval(j, x, r, DATA_PTR(m_mask), rdt);
    }

    // then over the connector domains
    for (d = m_connect.begin(); d != m_connect.end(); ++d) {
        (*d)->eval(j, x, r, DATA_PTR(m_mask), rdt);
    }

    // increment counter and time
    if (count) {
        clock_t t1 = clock();
        m_evaltime += double(t1 - t0)/CLOCKS_PER_SEC;
        m_nevals++;
    }
}

doublereal OneDim::ssnorm(doublereal* x, doublereal* r)
{
    eval(npos, x, r, 0.0, 0);
    doublereal ss = 0.0;
    for (size_t i = 0; i < m_size; i++) {
        ss = std::max(fabs(r[i]),ss);
    }
    return ss;
}

void OneDim::initTimeInteg(doublereal dt, doublereal* x)
{
    doublereal rdt_old = m_rdt;
    m_rdt = 1.0/dt;

    // if the stepsize has changed, then update the transient
    // part of the Jacobian
    if (fabs(rdt_old - m_rdt) > Tiny) {
        m_jac->updateTransient(m_rdt, DATA_PTR(m_mask));
    }

    // iterate over all domains, preparing each one to begin
    // time stepping
    Domain1D* d = left();
    while (d) {
        d->initTimeInteg(dt, x);
        d = d->right();
    }
}

void OneDim::setSteadyMode()
{
    m_rdt = 0.0;
    m_jac->updateTransient(m_rdt, DATA_PTR(m_mask));

    // iterate over all domains, preparing them for steady-state solution
    Domain1D* d = left();
    while (d) {
        d->setSteadyMode();
        d = d->right();
    }
}

void OneDim::init()
{
    if (!m_init) {
        Domain1D* d = left();
        while (d) {
            d->init();
            d = d->right();
        }
    }
    m_init = true;
}

void Domain1D::needJacUpdate()
{
    if (m_container) {
        m_container->jacobian().setAge(10000);
        m_container->saveStats();
    }
}

doublereal OneDim::timeStep(int nsteps, doublereal dt, doublereal* x,
                            doublereal* r, int loglevel)
{
    // set the Jacobian age parameter to the transient value
    newton().setOptions(m_ts_jac_age);

    writelog("\n\n step    size (s)    log10(ss) \n", loglevel);
    writelog("===============================\n", loglevel);

    int n = 0, m;
    doublereal ss;
    char str[80];
    while (n < nsteps) {
        if (loglevel > 0) {
            ss = ssnorm(x, r);
            sprintf(str, " %4d  %10.4g  %10.4g" , n,dt,log10(ss));
            writelog(str);
        }

        // set up for time stepping with stepsize dt
        initTimeInteg(dt,x);

        // solve the transient problem
        m = solve(x, r, loglevel-1);

        // successful time step. Copy the new solution in r to
        // the current solution in x.
        if (m >= 0) {
            n += 1;
            writelog("\n", loglevel);
            copy(r, r + m_size, x);
            if (m == 100) {
                dt *= 1.5;
            }
            //                 else dt /= 1.5;
            if (dt > m_tmax) {
                dt = m_tmax;
            }
        }

        // No solution could be found with this time step.
        // Decrease the stepsize and try again.
        else {
            writelog("...failure.\n", loglevel);
            dt *= m_tfactor;
            if (dt < m_tmin)
                throw CanteraError("OneDim::timeStep",
                                   "Time integration failed.");
        }
    }

    // Prepare to solve the steady problem.
    setSteadyMode();
    newton().setOptions(m_ss_jac_age);

    // return the value of the last stepsize, which may be smaller
    // than the initial stepsize
    return dt;
}

void OneDim::save(const std::string& fname, std::string id,
                  const std::string& desc, doublereal* sol,
                  int loglevel)
{
    struct tm* newtime;
    time_t aclock;
    ::time(&aclock);                /* Get time in seconds */
    newtime = localtime(&aclock);   /* Convert time to struct tm form */

    XML_Node root("doc");
    ifstream fin(fname.c_str());
    XML_Node* ct;
    if (fin) {
        root.build(fin);
        const XML_Node* same_ID = root.findID(id);
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
    } else {
        ct = &root.addChild("ctml");
    }
    XML_Node& sim = (XML_Node&)ct->addChild("simulation");
    sim.addAttribute("id",id);
    addString(sim,"timestamp",asctime(newtime));
    if (desc != "") {
        addString(sim,"description",desc);
    }

    Domain1D* d = left();
    while (d) {
        d->save(sim, sol);
        d = d->right();
    }
    ofstream s(fname.c_str());
    if (!s) {
        throw CanteraError("save","could not open file "+fname);
    }
    ct->write(s);
    s.close();
    writelog("Solution saved to file "+fname+" as solution "+id+".\n", loglevel);
}

void Domain1D::setGrid(size_t n, const doublereal* z)
{
    m_z.resize(n);
    m_points = n;
    for (size_t j = 0; j < m_points; j++) {
        m_z[j] = z[j];
    }
}

}
