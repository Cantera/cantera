//! @file ReactorNet.cpp
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

ReactorNet::ReactorNet() : Cantera::FuncEval(), m_nr(0), m_nreactors(0),
    m_integ(0), m_time(0.0), m_init(false),
    m_nv(0), m_rtol(1.0e-9), m_rtolsens(1.0e-4),
    m_atols(1.0e-15), m_atolsens(1.0e-4),
    m_maxstep(-1.0), m_maxErrTestFails(0),
    m_verbose(false), m_ntotpar(0)
{
#ifdef DEBUG_MODE
    m_verbose = true;
#endif
    m_integ = newIntegrator("CVODE");// CVodeInt;

    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator

    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_integ->setIterator(Newton_Iter);
}

ReactorNet::~ReactorNet()
{
    for (size_t n = 0; n < m_nr; n++) {
        if (m_iown[n]) {
            delete m_r[n];
        }
        m_r[n] = 0;
    }
    m_r.clear();
    m_reactors.clear();
    delete m_integ;
}

void ReactorNet::initialize()
{
    size_t n, nv;
    char buf[100];
    m_nv = 0;
    m_reactors.clear();
    m_nreactors = 0;
    writelog("Initializing reactor network.\n", m_verbose);
    if (m_nr == 0)
        throw CanteraError("ReactorNet::initialize",
                           "no reactors in network!");
    size_t sensParamNumber = 0;
    for (n = 0; n < m_nr; n++) {
        if (m_r[n]->type() >= ReactorType) {
            m_r[n]->initialize(m_time);
            Reactor* r = (Reactor*)m_r[n];
            m_reactors.push_back(r);
            nv = r->neq();
            m_size.push_back(nv);
            m_nparams.push_back(r->nSensParams());
            std::vector<std::pair<void*, int> > sens_objs = r->getSensitivityOrder();
            for (size_t i = 0; i < sens_objs.size(); i++) {
                std::map<size_t, size_t>& s = m_sensOrder[sens_objs[i]];
                for (std::map<size_t, size_t>::iterator iter = s.begin();
                        iter != s.end();
                        ++iter) {
                    m_sensIndex.resize(std::max(iter->second + 1, m_sensIndex.size()));
                    m_sensIndex[iter->second] = sensParamNumber++;
                }
            }
            m_nv += nv;
            m_nreactors++;

            if (m_verbose) {
                sprintf(buf,"Reactor %s: %s variables.\n",
                        int2str(n).c_str(), int2str(nv).c_str());
                writelog(buf);
                sprintf(buf,"            %s sensitivity params.\n",
                        int2str(r->nSensParams()).c_str());
                writelog(buf);
            }
            if (m_r[n]->type() == FlowReactorType && m_nr > 1) {
                throw CanteraError("ReactorNet::initialize",
                                   "FlowReactors must be used alone.");
            }
        }
    }

    m_connect.resize(m_nr*m_nr,0);
    m_ydot.resize(m_nv,0.0);
    size_t i, j, nin, nout, nw;
    ReactorBase* r, *rj;
    for (i = 0; i < m_nr; i++) {
        r = m_reactors[i];
        for (j = 0; j < m_nr; j++) {
            if (i == j) {
                connect(i,j);
            } else {
                rj = m_reactors[j];
                nin = rj->nInlets();
                for (n = 0; n < nin; n++) {
                    if (&rj->inlet(n).out() == r) {
                        connect(i,j);
                    }
                }
                nout = rj->nOutlets();
                for (n = 0; n < nout; n++) {
                    if (&rj->outlet(n).in() == r) {
                        connect(i,j);
                    }
                }
                nw = rj->nWalls();
                for (n = 0; n < nw; n++) {
                    if (&rj->wall(n).left() == rj
                            && &rj->wall(n).right() == r) {
                        connect(i,j);
                    } else if (&rj->wall(n).left() == r
                               && &rj->wall(n).right() == rj) {
                        connect(i,j);
                    }
                }
            }
        }
    }

    m_atol.resize(neq());
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), DATA_PTR(m_atol));
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
    m_integ->setMaxStepSize(m_maxstep);
    m_integ->setMaxErrTestFails(m_maxErrTestFails);
    if (m_verbose) {
        sprintf(buf, "Number of equations: %s\n", int2str(neq()).c_str());
        writelog(buf);
        sprintf(buf, "Maximum time step:   %14.6g\n", m_maxstep);
        writelog(buf);
    }
    m_integ->initialize(m_time, *this);
    m_init = true;
}

void ReactorNet::advance(doublereal time)
{
    if (!m_init) {
        if (m_maxstep < 0.0) {
            m_maxstep = time - m_time;
        }
        initialize();
    }
    m_integ->integrate(time);
    m_time = time;
    updateState(m_integ->solution());
}

double ReactorNet::step(doublereal time)
{
    if (!m_init) {
        if (m_maxstep < 0.0) {
            m_maxstep = time - m_time;
        }
        initialize();
    }
    m_time = m_integ->step(time);
    updateState(m_integ->solution());
    return m_time;
}

void ReactorNet::addReactor(ReactorBase* r, bool iown)
{
    r->setNetwork(this);
    if (r->type() >= ReactorType) {
        m_r.push_back(r);
        m_iown.push_back(iown);
        m_nr++;
        writelog("Adding reactor "+r->name()+"\n", m_verbose);
    } else {
        writelog("Not adding reactor "+r->name()+
                 ", since type = "+int2str(r->type())+"\n", m_verbose);
    }
}

void ReactorNet::eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p)
{
    size_t n;
    size_t start = 0;
    size_t pstart = 0;

    updateState(y);
    for (n = 0; n < m_nreactors; n++) {
        m_reactors[n]->evalEqs(t, y + start,
                               ydot + start, p + pstart);
        start += m_size[n];
        pstart += m_nparams[n];
    }
}

void ReactorNet::evalJacobian(doublereal t, doublereal* y,
                              doublereal* ydot, doublereal* p, Array2D* j)
{
    doublereal ysave, dy;
    Array2D& jac = *j;

    // use a try... catch block, since exceptions are not passed
    // through CVODE, since it is C code
    try {
        //evaluate the unperturbed ydot
        eval(t, y, ydot, p);
        for (size_t n = 0; n < m_nv; n++) {

            // perturb x(n)
            ysave = y[n];
            dy = m_atol[n] + fabs(ysave)*m_rtol;
            y[n] = ysave + dy;
            dy = y[n] - ysave;

            // calculate perturbed residual
            eval(t, y, DATA_PTR(m_ydot), p);

            // compute nth column of Jacobian
            for (size_t m = 0; m < m_nv; m++) {
                jac(m,n) = (m_ydot[m] - ydot[m])/dy;
            }
            y[n] = ysave;
        }
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        error("Terminating execution.");
    }
}

void ReactorNet::updateState(doublereal* y)
{
    size_t start = 0;
    for (size_t n = 0; n < m_nreactors; n++) {
        m_reactors[n]->updateState(y + start);
        start += m_size[n];
    }
}

void ReactorNet::getInitialConditions(doublereal t0,
                                      size_t leny, doublereal* y)
{
    size_t start = 0;
    for (size_t n = 0; n < m_nreactors; n++) {
        m_reactors[n]->getInitialConditions(t0, m_size[n], y + start);
        start += m_size[n];
    }
}

size_t ReactorNet::globalComponentIndex(const string& species, size_t reactor)
{
    if (!m_init) {
        initialize();
    }
    size_t start = 0;
    size_t n;
    for (n = 0; n < reactor; n++) {
        start += m_size[n];
    }
    return start + m_reactors[n]->componentIndex(species);
}

void ReactorNet::registerSensitivityReaction(void* reactor,
        size_t reactionIndex, const std::string& name, int leftright)
{
    std::pair<void*, int> R = std::make_pair(reactor, leftright);
    if (m_sensOrder.count(R) &&
            m_sensOrder[R].count(reactionIndex)) {
        throw CanteraError("ReactorNet::registerSensitivityReaction",
                           "Attempted to register duplicate sensitivity reaction");
    }
    m_paramNames.push_back(name);
    m_sensOrder[R][reactionIndex] = m_ntotpar;
    m_ntotpar++;
}

}
