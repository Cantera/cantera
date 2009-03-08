/**
 *  @file Func1.h
 *
 *  $Author: dggoodwin $
 *  $Date: 2005/11/10 15:06:33 $
 *  $Revision: 1.6 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FUNC1_H
#define CT_FUNC1_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"

namespace Cantera {

    const int FourierFuncType = 1;
    const int PolyFuncType = 2;
    const int ArrheniusFuncType = 3;
    const int GaussianFuncType = 4;
    const int SumFuncType = 20;
    const int DiffFuncType = 25;
    const int ProdFuncType = 30;
    const int RatioFuncType = 40;
    const int PeriodicFuncType = 50;

    /**
     * Base class for 'functor' classes that evaluate a function of
     * one variable. 
     */
    class Func1 {
    public:
        Func1() {}
        virtual ~Func1() {}
        /// Calls method eval to evaluate the function
        doublereal operator()(doublereal t) { return eval(t); }
        /// Evaluate the function.
        virtual doublereal eval(doublereal t) { return 0.0; }
    private:
    };


    /**
     * A Gaussian. 
     * \f[
     * f(t) = A e^{-[(t - t_0)/\tau]^2}
     * \f]
     * where \f[ \tau = \frac{fwhm}{2\sqrt{\ln 2}} \f]
     * @param A peak value
     * @param t0 offset
     * @param fwhm full width at half max  
     */
    class Gaussian : public Func1 {
    public:
        Gaussian(double A, double t0, double fwhm) {
            m_A = A;
            m_t0 = t0;
            m_tau = fwhm/(2.0*sqrt(log(2.0)));
        }
        virtual ~Gaussian() {}
        virtual doublereal eval(doublereal t) {
            doublereal x = (t - m_t0)/m_tau;
            return m_A*exp(-x*x);
        }
    protected:
        doublereal m_A, m_t0, m_tau;
    private:
    };


    /**
     * Polynomial of degree n.  
     */    
    class Poly1 : public Func1 {
    public:
        Poly1(int n, doublereal* c) {
            m_n = n+1;
            m_c.resize(n+1);
            copy(c, c+m_n, m_c.begin());
        }
        virtual ~Poly1() {}

        virtual doublereal eval(doublereal t) {
            int n;
            doublereal r = m_c[m_n-1];
            for (n = 1; n < m_n; n++) {
                r *= t;
                r += m_c[m_n - n - 1];
            }
            return r;
        }

    protected:
        int m_n;
        vector_fp m_c;
    };


    /**
     * Fourier cosine/sine series.
     * 
     * \f[
     * f(t) = \frac{A_0}{2} + 
     * \sum_{n=1}^N A_n \cos (n \omega t) + B_n \sin (n \omega t)
     * \f]
     */
    class Fourier1 : public Func1 {
    public:
        Fourier1(int n, doublereal omega, doublereal a0, 
            doublereal* a, doublereal* b) {
            m_n = n;
            m_omega = omega;
            m_a0_2 = 0.5*a0;
            m_ccos.resize(n);
            m_csin.resize(n);
            copy(a, a+n, m_ccos.begin());
            copy(b, b+n, m_csin.begin());
        }
        virtual ~Fourier1() {}

        virtual doublereal eval(doublereal t) {
            int n, nn;
            doublereal sum = m_a0_2;
            for (n = 0; n < m_n; n++) {
                nn = n + 1;
                sum += m_ccos[n]*cos(m_omega*nn*t) 
                       + m_csin[n]*sin(m_omega*nn*t);
            }
            return sum;
        }

    protected:
        int m_n;
        doublereal m_omega, m_a0_2;
        vector_fp m_ccos, m_csin;
    };


    /**
     * Sum of Arrhenius terms.
     * \f[
     * f(T) = \sum_{n=1}^N A_n T^b_n \exp(-E_n/T)
     * \f]
     */
    class Arrhenius1 : public Func1 {
    public:
        Arrhenius1(int n, doublereal* c) {
            m_n = n;
            m_A.resize(n);
            m_b.resize(n);
            m_E.resize(n);
            int loc;
            for (int i = 0; i < n; i++) {
                loc = 3*i;
                m_A[i] = c[loc];
                m_b[i] = c[loc+1];
                m_E[i] = c[loc+2];
            }
        }
        virtual ~Arrhenius1() {}

        virtual doublereal eval(doublereal t) {
            int n;
            doublereal sum = 0.0;
            for (n = 0; n < m_n; n++) {
                sum += m_A[n]*pow(t,m_b[n])*exp(-m_E[n]/t);
            }
            return sum;
        }

    protected:
        int m_n;
        vector_fp m_A, m_b, m_E;
    };

    /**
     *  Periodic function. Takes any function and makes it
     *  periodic with period T.
     */
    class PeriodicFunc : public Func1 {
    public:
        PeriodicFunc(Func1& f, doublereal T) {
            m_func = &f;
            m_period = T;
        }
        virtual ~PeriodicFunc() {}
        virtual doublereal eval(doublereal t) {
            int np = int(t/m_period);
            doublereal time = t - np*m_period;
            return m_func->eval(time);
        }
    protected:
        Func1* m_func;
        doublereal m_period;

    private:
    };
            

    /**
     * Sum of two functions.
     */
    class Func1Sum : public Func1 {
    public:
        Func1Sum(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Func1Sum() {}
        virtual doublereal eval(doublereal t) {
            return m_f1->eval(t) + m_f2->eval(t);
        }
    protected:
        Func1 *m_f1, *m_f2;
    private:
    };


    /**
     * Difference of two functions.
     */
    class Func1Diff : public Func1 {
    public:
        Func1Diff(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Func1Diff() {}
        virtual doublereal eval(doublereal t) {
            return m_f1->eval(t) - m_f2->eval(t);
        }
    protected:
        Func1 *m_f1, *m_f2;
    private:
    };


    /**
     * Product of two functions.
     */
    class Func1Product : public Func1 {
    public:
        Func1Product(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Func1Product() {}
        virtual doublereal eval(doublereal t) {
            return m_f1->eval(t) * m_f2->eval(t);
        }
    protected:
        Func1 *m_f1, *m_f2;
    private:
    };


    /**
     * Ratio of two functions.
     */
    class Func1Ratio : public Func1 {
    public:
        Func1Ratio(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Func1Ratio() {}
        virtual doublereal eval(doublereal t) {
            return m_f1->eval(t) / m_f2->eval(t);
        }
    protected:
        Func1 *m_f1, *m_f2;
    private:
    };
}

#endif
