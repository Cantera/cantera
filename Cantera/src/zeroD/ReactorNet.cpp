#include "ReactorNet.h"
#include "../CVode.h"

namespace Cantera {

    ReactorNet::ReactorNet() : FuncEval(), m_nr(0), m_nreactors(0),
                               m_integ(0), m_init(false), m_nv(0), m_rtol(1.0e-6) {
        m_integ = new CVodeInt;

        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator

        m_integ->setMethod(BDF_Method);
        m_integ->setProblemType(DENSE + NOJAC);
        m_integ->setIterator(Newton_Iter);        
    }

    void ReactorNet::initialize(doublereal t0) {
        int n, nv;
        m_nv = 0;
        m_reactors.clear();
        for (n = 0; n < m_nr; n++) {
            m_r[n]->initialize(t0);
            if (m_r[n]->type() == ReactorType) {
                Reactor* r = (Reactor*)m_r[n];
                m_reactors.push_back(r);
                nv = m_reactors[n]->neq();
                m_size.push_back(nv);
                m_nv += nv;
                m_nreactors++;
            }
        }
        m_atol.resize(neq());
        fill(m_atol.begin(), m_atol.end(), 1.e-15);
        m_integ->setTolerances(m_rtol, neq(), m_atol.begin());
        m_integ->setMaxStep(m_maxstep);
        m_integ->initialize(t0, *this);
        m_init = true;
    }

    void ReactorNet::advance(doublereal time) {
        if (!m_init) {
            m_maxstep = time;
            initialize();
        }
        m_integ->integrate(time);
        m_time = time;
        updateState(m_integ->solution());
    }

    double ReactorNet::step(doublereal time) {
        if (!m_init) {
            m_maxstep = time;
            initialize();
        }
        m_time = m_integ->step(time);
        updateState(m_integ->solution());
        return m_time;
    }

    void ReactorNet::eval(doublereal t, doublereal* y, doublereal* ydot) {
        int n;
        int start = 0;
        updateState(y);
        for (n = 0; n < m_nreactors; n++) {
            m_reactors[n]->evalEqs(t, y + start, ydot + start);
            start += m_size[n];
        }
    }

    void ReactorNet::updateState(doublereal* y) {
        int n;
        int start = 0;
        for (n = 0; n < m_nreactors; n++) {
            m_reactors[n]->updateState(y + start);
            start += m_size[n];
        }
    }

    void ReactorNet::getInitialConditions(doublereal t0, 
        size_t leny, doublereal* y) {
        int n;
        int start = 0;
        for (n = 0; n < m_nreactors; n++) {
            m_reactors[n]->getInitialConditions(t0, m_size[n], y + start);
            start += m_size[n];
        }
    }
}

