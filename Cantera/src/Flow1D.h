
deprecated


#ifndef CT_FLOW1D_H
#define CT_FLOW1D_H

#include "MultiDomain.h"
#include "MultiJac.h"
#include "MultiNewton.h"

#include "Jac1D.h"
#include "Newton1D.h"

#include "Surf1D.h"

namespace Cantera {

    /**
     * Container class for multiple-domain 1D problems.
     */
    class OneDim : public MultiDomain {
    public: 
        OneDim() : m_jac(0), m_newt(0) {
            m_newt = new MultiNewton(1);
        }

        virtual void addDomain(Resid1D* d) {
            MultiDomain::addDomain(d);
            m_newt->resize(size());
            delete m_jac;
            m_jac = 0;
            m_jac = new MultiJac(*this);
            m_jac_ok = false;
            int nd = m_dom.size();
            for (int i = 0; i < nd; i++)
                m_dom[i]->setJac(m_jac);
        }

        virtual ~OneDim() {
            delete m_jac;
            delete m_newt;
        }

        MultiJac& jacobian() { return *m_jac; }
        MultiNewton& newton() { return *m_newt; }

        
        virtual int solve(doublereal* x, doublereal* xnew, int loglevel) {
            if (!m_jac) {
                cout << "creating new jac.." << endl;
                m_jac = new MultiJac(*this);
                m_jac_ok = false;
                int nd = m_dom.size();
                for (int i = 0; i < nd; i++)
                    m_dom[i]->setJac(m_jac);
                cout << "done" << endl;
            }

            if (!m_jac_ok) {
                cout << "eval jac" << endl;
                eval(-1, x, xnew, m_rdt);
                m_jac->eval(x, xnew, m_rdt);
                m_jac_ok = true;
            }
            return m_newt->solve(x, xnew, *this, *m_jac, loglevel);
        }



    protected:
        MultiJac* m_jac;
        MultiNewton* m_newt;
        Jac1D*  m_jac1;
        Newton1D* m_newt1;
    };

}

#endif


