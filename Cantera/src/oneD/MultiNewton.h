/**
 *
 *  @file MultiNewton.h
 */

/*
 *  $Author: dggoodwin $
 *  $Date: 2006/11/27 21:43:34 $
 *  $Revision: 1.2 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_MULTINEWTON_H
#define CT_MULTINEWTON_H

#include "MultiJac.h"

namespace Cantera {

    /**
     * Newton iterator for multi-domain, one-dimensional problems.
     * Used by class OneDim.
     */
    class MultiNewton {

    public:

        MultiNewton(int sz);
        virtual ~MultiNewton();

        int size() { return m_n; }

        /// Compute undamped step
        void step(doublereal* x, doublereal* step, 
            OneDim& r, MultiJac& jac, int loglevel);

        /// Compute factor to keep all components in bounds.
        doublereal boundStep(const doublereal* x0, const doublereal* step0,
            const OneDim& r, int loglevel);

        int dampStep(const doublereal* x0, const doublereal* step0, 
            doublereal* x1, doublereal* step1, doublereal& s1, 
            OneDim& r, MultiJac& jac, int loglevel, bool writetitle);

        doublereal norm2(const doublereal* x, const doublereal* step, 
            OneDim& r) const;

        int solve(doublereal* x0, doublereal* x1, OneDim& r, MultiJac& jac,
            int loglevel);

        /// Set options.
        void setOptions(int maxJacAge = 5) {m_maxAge = maxJacAge;}

        /// Change the problem size.
        void resize(int points);


    protected:

        doublereal* getWorkArray();
        void releaseWorkArray(doublereal* work);
        std::vector<doublereal*> m_workarrays;
        int m_maxAge;
        int m_nv, m_np, m_n;
        doublereal m_elapsed;

    private:

        char m_buf[100];
    };
}

#endif


