
#ifndef CT_FTNODESYS_H
#define CT_FTNODESYS_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

namespace Cantera {


    typedef void (*Ftn_RHS_Func) (integer* n, doublereal* t, doublereal* y,
        doublereal* ydot);

    /**
     * A class to integrate a system of ODE's defined by a Fortran function.
     * Not currently used.
     */
    class FtnODESys : public FuncEval {

    public:
        
        FtnODESys(int n, doublereal* y0, Ftn_RHS_Func f) {
            m_func = f;
            m_n = n;
            m_y0.resize(n);
            int i;
            for (i = 0; i < n; i++) m_y0[i] = y0[i];
        }

        virtual ~FtnODESys(){}

	virtual void eval(double t, double* y, double* ydot) {
            m_func(&m_n, &t, y, ydot);
        }

        virtual void getInitialConditions(double t0, double* y) {
            int i;
            for (i = 0; i < m_n; i++) y[i] = m_y0[i];
        }

        virtual int neq() { return m_n; }

    protected:
        int m_n;
        vector_fp m_y0;
        Ftn_RHS_Func m_func;

    private:

    };

}

#endif
