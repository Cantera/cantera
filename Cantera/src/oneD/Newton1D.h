/**
 *
 *  @file Newton1D.h
 *
 *  Newton solver >>> under construction! <<<<
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_NEWTON1D_H
#define CT_NEWTON1D_H

#include "Jac1D.h"

namespace Cantera {

    class Newton1D {

    public:

        Newton1D(int nv, int np);
        virtual ~Newton1D();

        doublereal norm(const doublereal* step);
        void step(doublereal* x, doublereal* step, 
            Domain1D& r, Jac1D& jac, int loglevel, int update=1);
        doublereal boundStep(const doublereal* x0, const doublereal* step0,
            const Domain1D& r, int loglevel);
        int dampStep(const doublereal* x0, const doublereal* step0, 
            doublereal* x1, doublereal* step1, doublereal& s1, 
            Domain1D& r, Jac1D& jac, int loglevel, bool writetitle);
        void getErrorWeights(const doublereal* x, doublereal* ewt, Domain1D& r);
        doublereal norm2(const doublereal* step, doublereal* ewt);
        doublereal norm_infty(const doublereal* step, doublereal* ewt);
        int solve(doublereal* x0, doublereal* x1, Domain1D& r, Jac1D& jac,
            int loglevel);
        int timeIntegrate(int n, doublereal dt, 
            doublereal* x0, doublereal* x1, 
            Domain1D& r, Jac1D& jac, int loglevel);
        doublereal ssnorm(doublereal* x, doublereal* resid, Domain1D& r);

        void setOptions(int maxJacAge = 5, doublereal maxNormRatio = 0.001) {
            m_maxAge = maxJacAge;
            m_maxRatio = maxNormRatio;
        }
        void resize(int points);

    protected:

        doublereal* getWorkArray();
        void releaseWorkArray(doublereal* work);
        vector<doublereal*> m_workarrays;
        vector_fp m_ewt;
        int m_maxAge;
        int m_maxRatio;
        int m_nv, m_np, m_n;

    private:

        size_t index(int n, int j) {
            return m_nv * j + n;
        }
    };
}

#endif


