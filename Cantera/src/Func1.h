/**
 *  @file Func1.h
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FUNC1_H
#define CT_FUNC1_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"

#include <iostream>
#include <string>
using namespace std;

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
    const int CompositeFuncType = 60;
    const int TimesConstantFuncType = 70;
    const int PlusConstantFuncType = 80;
    const int SinFuncType = 100;
    const int CosFuncType = 102;
    const int ExpFuncType = 104;
    const int PowFuncType = 106;
    const int ConstFuncType = 110;

    class Sin1;
    class Cos1;
    class Exp1;
    class Pow1;
    class TimesConstant1;

    /**
     * Base class for 'functor' classes that evaluate a function of
     * one variable.
     */
    class Func1 {
    public:
        Func1() : m_c(0.0), m_f1(0), m_f2(0), m_parent(0) {}
        virtual ~Func1() {}
        virtual int ID() const { return 0; }

        virtual Func1& duplicate() { cout << "DUPL ERR: ID = " << ID() << endl; 
            return *(new Func1);}

        /// Calls method eval to evaluate the function
        doublereal operator()(doublereal t) const { return eval(t); }
        
        /// Evaluate the function.
        virtual doublereal eval(doublereal t) const { return 0.0; }

        virtual Func1& derivative() const {
            cout << "derivative error... ERR: ID = " << ID() << endl;
            cout << write("x") << endl;
            return *(new Func1);
        }

        bool isIdentical(Func1& other) const {
            if ((ID() != other.ID()) || (m_c != other.m_c))
                return false;
            if (m_f1) {
                if (!other.m_f1) return false;
                if (!m_f1->isIdentical(*other.m_f1)) return false;
            }
            if (m_f2) {
                if (!other.m_f2) return false;
                if (!m_f2->isIdentical(*other.m_f2)) return false;
            }
            return true;
        }

        virtual doublereal isProportional(TimesConstant1& other);
        virtual doublereal isProportional(Func1& other);

        virtual std::string write(std::string arg) const;

        doublereal c() const { return m_c; }
        void setC(doublereal c) { m_c = c; }
        Func1& func1() { return *m_f1; }
        Func1& func2() { return *m_f2; }
        virtual int order() const { return 3; }
        Func1& func1_dup() const { return m_f1->duplicate(); }
        Func1& func2_dup() const { return m_f2->duplicate(); }
        Func1* parent() { return m_parent; }
        void setParent(Func1* p) { m_parent = p; }

    protected:
        doublereal m_c;
        Func1 *m_f1, *m_f2;
        Func1* m_parent;

    private:
    };


    Func1& newSumFunction(Func1& f1, Func1& f2);
    Func1& newDiffFunction(Func1& f1, Func1& f2);
    Func1& newProdFunction(Func1& f1, Func1& f2);
    Func1& newRatioFunction(Func1& f1, Func1& f2);
    Func1& newCompositeFunction(Func1& f1, Func1& f2);
    Func1& newTimesConstFunction(Func1& f1, doublereal c);
    Func1& newPlusConstFunction(Func1& f1, doublereal c);

    /// sin
    class Sin1 : public Func1 {
    public:
        Sin1(doublereal omega = 1.0) {
            m_c = omega;
        }
        virtual ~Sin1() {}
        virtual std::string write(std::string arg) const;
        virtual Func1& duplicate() { return *(new Sin1(m_c)); }
        virtual int ID() const { return SinFuncType; }
        virtual doublereal eval(doublereal t) const {
            return sin(m_c*t);
        }
        virtual Func1& derivative() const;

    protected:

    };

    /// cos
    class Cos1 : public Func1 {
    public:
        Cos1(doublereal omega = 1.0) {
            m_c = omega;
        }
        virtual ~Cos1() {}
        virtual std::string write(std::string arg) const;
        virtual Func1& duplicate() { return *(new Cos1(m_c)); }
        virtual int ID() const { return CosFuncType; }
        virtual doublereal eval(doublereal t) const {
            return cos(m_c * t);
        }
        virtual Func1& derivative() const;

    protected:
    };

    /// exp
    class Exp1 : public Func1 {
    public:
        Exp1(doublereal A = 1.0) {m_c = A;}
        virtual ~Exp1() {}
        virtual std::string write(std::string arg) const;
        virtual int ID() const { return ExpFuncType; }
        virtual Func1& duplicate() { return *(new Exp1(m_c)); }
        virtual doublereal eval(doublereal t) const {
            return exp(m_c*t);
        }

        virtual Func1& derivative() const;

    protected:

    };

    /// pow
    class Pow1 : public Func1 {
    public:
        Pow1(doublereal n) {m_c = n;}
        virtual ~Pow1() {}
        virtual std::string write(std::string arg) const;
        virtual int ID() const { return PowFuncType; }
        virtual Func1& duplicate() { return *(new Pow1(m_c)); }
        virtual doublereal eval(doublereal t) const {
            return pow(t, m_c);
        }
        virtual Func1& derivative() const;

    protected:

    };

    /**
     * Constant.
     */    
    class Const1 : public Func1 {
    public:
        Const1(doublereal A) {
            m_c = A;
        }
        virtual ~Const1() {}
        virtual std::string write(std::string arg) const;
        virtual int ID() const { return ConstFuncType; }
        virtual doublereal eval(doublereal t) const {
            return m_c;
        }
        virtual Func1& duplicate() { return *(new Const1(m_c)); }
        virtual Func1& derivative() const {
            Func1* z = new Const1(0.0);
            return *z;
        }

    protected:
    };



    /**
     * Sum of two functions.
     */
    class Sum1 : public Func1 {
    public:
        Sum1(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
            if (m_f1 == m_f2) 
                cout << "Same functions!" << endl;
            m_f1->setParent(this);
            m_f2->setParent(this);
        }
        virtual ~Sum1() {
            delete m_f1;
            delete m_f2;
        }
        virtual int ID() const { return SumFuncType; }
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) + m_f2->eval(t);
        }
        virtual Func1& duplicate() {
            Func1& f1d = m_f1->duplicate();
            Func1& f2d = m_f2->duplicate();
            Func1& dup = newSumFunction(f1d, f2d);
            return dup;
        }
        virtual Func1& derivative() const {
            Func1& d1 = m_f1->derivative();
            Func1& d2 = m_f2->derivative();
            Func1& d = newSumFunction(d1, d2);
            return d;
        }
        virtual int order() const { return 0; }
        virtual std::string write(std::string arg) const;            
    protected:
    };


    /**
     * Difference of two functions.
     */
    class Diff1 : public Func1 {
    public:
        Diff1(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Diff1() {
            delete m_f1;
            delete m_f2;
        }
        virtual int ID() const { return DiffFuncType; }
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) - m_f2->eval(t);
        }
        virtual Func1& duplicate() {
            Func1& f1d = m_f1->duplicate();
            Func1& f2d = m_f2->duplicate();
            Func1& dup = newDiffFunction(f1d, f2d);
            return dup;
        }
        virtual Func1& derivative() const {
            Func1& d = newDiffFunction(m_f1->derivative(), m_f2->derivative());
            return d;
        }
        virtual int order() const { return 0; }
        virtual std::string write(std::string arg) const;

    protected:
    };


    /**
     * Product of two functions.
     */
    class Product1 : public Func1 {
    public:
        Product1(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }

        virtual ~Product1() {
            //cout << "In Product1 destructor, deleting" << m_f1 << " " << m_f2 << endl;
            delete m_f1;
            delete m_f2;
        }
        virtual int ID() const { return ProdFuncType; }
        virtual Func1& duplicate() {
            Func1& f1d = m_f1->duplicate();
            Func1& f2d = m_f2->duplicate();
            Func1& dup = newProdFunction(f1d, f2d);
            return dup;
        }
        virtual std::string write(std::string arg) const;
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) * m_f2->eval(t);
        }
        virtual Func1& derivative() const {
            Func1& a1 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
            Func1& a2 = newProdFunction(m_f2->duplicate(), m_f1->derivative());
            Func1& s = newSumFunction(a1, a2);
            return s;
        }
        virtual int order() const { return 1; }
                
    protected:
    };

    /**
     * Product of two functions.
     */
    class TimesConstant1 : public Func1 {
    public:
        TimesConstant1(Func1& f1, doublereal A) {
            m_f1 = &f1;
            m_c = A;
        }

        virtual ~TimesConstant1() {
            delete m_f1;
        }
        virtual int ID() const { return TimesConstantFuncType; }
        virtual Func1& duplicate() {
            Func1& f1 = m_f1->duplicate();
            Func1* dup = new TimesConstant1(f1, m_c);
            return *dup;
        }
        virtual doublereal isProportional(TimesConstant1& other) {
            if (func1().isIdentical(other.func1())) 
                return (other.c()/c());
            else
                return 0.0;
        }
        virtual doublereal isProportional(Func1& other) {
            if (func1().isIdentical(other)) return 1.0/c();
            else return 0.0;
        }
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) * m_c;
        }
        virtual Func1& derivative() const {
            Func1& f1d = m_f1->derivative();
            Func1* d = new TimesConstant1(f1d, m_c);
            return *d;
        }
        virtual std::string write(std::string arg) const;
        virtual int order() const { return 0; }                
    protected:
    };

    /**
     * A function plus a constant.
     */
    class PlusConstant1 : public Func1 {
    public:
        PlusConstant1(Func1& f1, doublereal A) {
            m_f1 = &f1;
            m_c = A;
        }

        virtual ~PlusConstant1() {
            //cout << "PlusConstant1: deleting " << m_f1 << endl;
            delete m_f1;
        }
        virtual int ID() const { return PlusConstantFuncType; }
        virtual Func1& duplicate() {
            Func1& f1 = m_f1->duplicate();
            Func1* dup = new PlusConstant1(f1, m_c);
            return *dup;
        }

        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) + m_c;
        }
        virtual Func1& derivative() const {
            Func1& f1d = m_f1->derivative();
            return f1d;
        }
        virtual std::string write(std::string arg) const;
        virtual int order() const { return 0; }
                
    protected:
    };


    /**
     * Ratio of two functions.
     */
    class Ratio1 : public Func1 {
    public:
        Ratio1(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Ratio1() {
            //cout << "Ratio1: deleting " << m_f1 << " " << m_f2 << endl;
            delete m_f1;
            delete m_f2;
        }
        virtual int ID() const { return RatioFuncType; }
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval(t) / m_f2->eval(t);
        }
        virtual Func1& duplicate() {
            Func1& f1d = m_f1->duplicate();
            Func1& f2d = m_f2->duplicate();
            Func1& dup = newRatioFunction(f1d, f2d);
            return dup;
        }
        virtual Func1& derivative() const {
            Func1& a1 = newProdFunction(m_f1->derivative(), m_f2->duplicate());
            Func1& a2 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
            Func1& s = newDiffFunction(a1, a2);
            Func1& p = newProdFunction(m_f2->duplicate(), m_f2->duplicate());
            Func1& r = newRatioFunction(s, p);
            return r;
        }
        virtual std::string write(std::string arg) const;
        virtual int order() const { return 1; }

    protected:
    };

    /**
     * Composite function.
     */
    class Composite1 : public Func1 {
    public:
        Composite1(Func1& f1, Func1& f2) {
            m_f1 = &f1;
            m_f2 = &f2;
        }
        virtual ~Composite1() {
            delete m_f1;
            delete m_f2;
        }
        virtual int ID() const { return CompositeFuncType; }
        virtual doublereal eval(doublereal t) const {
            return m_f1->eval( m_f2->eval(t) );
        }
        virtual Func1& duplicate() {
            Func1& f1d = m_f1->duplicate();
            Func1& f2d = m_f2->duplicate();
            Func1& dup = newCompositeFunction(f1d, f2d);
            return dup;
        }
        virtual Func1& derivative() const {
            Func1& d1 = m_f1->derivative();
            Func1& d3 = newCompositeFunction(d1, m_f2->duplicate());
            Func1& d2 = m_f2->derivative();
            Func1& p = newProdFunction(d3, d2);
            return p;
        }
        virtual std::string write(std::string arg) const;
        virtual int order() const { return 2; }
    protected:
    };

    //
    // The functors below are the old-style ones. They still work, 
    // but can't do derivatives.
    //

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
            m_tau = fwhm/(2.0*std::sqrt(std::log(2.0)));
        }
        virtual ~Gaussian() {}
        virtual doublereal eval(doublereal t) const {
            doublereal x = (t - m_t0)/m_tau;
            return m_A*std::exp(-x*x);
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
            std::copy(c, c+m_n, m_c.begin());
        }
        virtual ~Poly1() {}

        virtual doublereal eval(doublereal t) const {
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
            std::copy(a, a+n, m_ccos.begin());
            std::copy(b, b+n, m_csin.begin());
        }
        virtual ~Fourier1() {}

        virtual doublereal eval(doublereal t) const {
            int n, nn;
            doublereal sum = m_a0_2;
            for (n = 0; n < m_n; n++) {
                nn = n + 1;
                sum += m_ccos[n]*std::cos(m_omega*nn*t)
                       + m_csin[n]*std::sin(m_omega*nn*t);
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

        virtual doublereal eval(doublereal t) const {
            int n;
            doublereal sum = 0.0;
            for (n = 0; n < m_n; n++) {
                sum += m_A[n]*std::pow(t,m_b[n])*std::exp(-m_E[n]/t);
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
    class Periodic1 : public Func1 {
    public:
        Periodic1(Func1& f, doublereal T) {
            m_func = &f;
            m_c = T;
        }
        virtual ~Periodic1() { delete m_func; }
        virtual doublereal eval(doublereal t) const {
            int np = int(t/m_c);
            doublereal time = t - np*m_c;
            return m_func->eval(time);
        }
    protected:
        Func1* m_func;

    private:
    };

}


#endif
