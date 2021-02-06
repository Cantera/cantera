/**
 *  @file Func1.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FUNC1_H
#define CT_FUNC1_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"

#include <iostream>

namespace Cantera
{

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
const int TabulatedFuncType = 120;

class TimesConstant1;

/**
 * Base class for 'functor' classes that evaluate a function of one variable.
 */
class Func1
{
public:
    Func1();

    virtual ~Func1() {}

    Func1(const Func1& right);

    Func1& operator=(const Func1& right);

    //! Duplicate the current function.
    /*!
     * This duplicates the current function, returning a reference to the newly
     * created function.
     */
    virtual Func1& duplicate() const;

    virtual int ID() const;

    //! Calls method eval to evaluate the function
    doublereal operator()(doublereal t) const;

    /// Evaluate the function.
    virtual doublereal eval(doublereal t) const;

    //! Creates a derivative to the current function
    /*!
     * This will create a new derivative function and return a reference to the
     * function.
     */
    virtual Func1& derivative() const;

    //! Routine to determine if two functions are the same.
    /*!
     * Two functions are the same if they are the same function. This means
     * that the ID and stored constant is the same. This means that the m_f1
     * and m_f2 are identical if they are non-null.
     */
    bool isIdentical(Func1& other) const;

    virtual doublereal isProportional(TimesConstant1& other);
    virtual doublereal isProportional(Func1& other);

    virtual std::string write(const std::string& arg) const;

    //! accessor function for the stored constant
    doublereal c() const;

    //! Function to set the stored constant
    void setC(doublereal c);

    //! accessor function for m_f1
    Func1& func1() const;

    //! accessor function for m_f2
    Func1& func2() const;

    //! Return the order of the function, if it makes sense
    virtual int order() const;

    Func1& func1_dup() const;

    Func1& func2_dup() const;

    Func1* parent() const;

    void setParent(Func1* p);

protected:
    doublereal m_c;
    Func1* m_f1;
    Func1* m_f2;
    Func1* m_parent;
};


Func1& newSumFunction(Func1& f1, Func1& f2);
Func1& newDiffFunction(Func1& f1, Func1& f2);
Func1& newProdFunction(Func1& f1, Func1& f2);
Func1& newRatioFunction(Func1& f1, Func1& f2);
Func1& newCompositeFunction(Func1& f1, Func1& f2);
Func1& newTimesConstFunction(Func1& f1, doublereal c);
Func1& newPlusConstFunction(Func1& f1, doublereal c);


//! implements the sin() function
/*!
 * The argument to sin() is in radians
 */
class Sin1 : public Func1
{
public:
    Sin1(doublereal omega = 1.0) :
        Func1() {
        m_c = omega;
    }

    Sin1(const Sin1& b) :
        Func1(b) {
    }

    Sin1& operator=(const Sin1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        return *this;
    }

    virtual Func1& duplicate() const {
        Sin1* nfunc = new Sin1(*this);
        return (Func1&) *nfunc;
    }

    virtual std::string write(const std::string& arg) const;

    virtual int ID() const {
        return SinFuncType;
    }

    virtual doublereal eval(doublereal t) const {
        return sin(m_c*t);
    }

    virtual Func1& derivative() const;
};


//! implements the cos() function
/*!
 * The argument to cos() is in radians
 */
class Cos1 : public Func1
{
public:
    Cos1(doublereal omega = 1.0) :
        Func1() {
        m_c = omega;
    }

    Cos1(const Cos1& b) :
        Func1(b) {
    }

    Cos1& operator=(const Cos1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        return *this;
    }

    virtual Func1& duplicate() const {
        Cos1* nfunc = new Cos1(*this);
        return (Func1&) *nfunc;
    }
    virtual std::string write(const std::string& arg) const;
    virtual int ID() const {
        return CosFuncType;
    }
    virtual doublereal eval(doublereal t) const {
        return cos(m_c * t);
    }
    virtual Func1& derivative() const;
};


//! implements the exponential function
class Exp1 : public Func1
{
public:
    Exp1(doublereal A = 1.0) :
        Func1() {
        m_c = A;
    }

    Exp1(const Exp1& b) :
        Func1(b) {
    }
    Exp1& operator=(const Exp1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        return *this;
    }
    virtual std::string write(const std::string& arg) const;
    virtual int ID() const {
        return ExpFuncType;
    }
    virtual Func1& duplicate() const {
        return *(new Exp1(m_c));
    }
    virtual doublereal eval(doublereal t) const {
        return exp(m_c*t);
    }

    virtual Func1& derivative() const;
};


//! implements the power function (pow)
class Pow1 : public Func1
{
public:
    Pow1(doublereal n) :
        Func1() {
        m_c = n;
    }

    Pow1(const Pow1& b) :
        Func1(b) {
    }
    Pow1& operator=(const Pow1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        return *this;
    }
    virtual std::string write(const std::string& arg) const;
    virtual int ID() const {
        return PowFuncType;
    }
    virtual Func1& duplicate() const {
        return *(new Pow1(m_c));
    }
    virtual doublereal eval(doublereal t) const {
        return pow(t, m_c);
    }
    virtual Func1& derivative() const;
};


//! The Tabulated1 class implements a tabulated function
class Tabulated1 : public Func1
{
public:
    //! Constructor.
    /*!
     * @param n      Size of tabulated value arrays
     * @param tvals   Pointer to time value array
     * @param fvals   Pointer to function value array
     * @param method Interpolation method ('linear' or 'previous')
     */
    Tabulated1(size_t n, const double* tvals, const double* fvals,
               const std::string& method = "linear");

    virtual std::string write(const std::string& arg) const;
    virtual int ID() const {
        return TabulatedFuncType;
    }
    virtual double eval(double t) const;
    virtual Func1& duplicate() const {
        if (m_isLinear) {
            return *(new Tabulated1(m_tvec.size(), &m_tvec[0], &m_fvec[0],
                                    "linear"));
        } else {
            return *(new Tabulated1(m_tvec.size(), &m_tvec[0], &m_fvec[0],
                                    "previous"));
        }
    }

    virtual Func1& derivative() const;
private:
    vector_fp m_tvec; //!< Vector of time values
    vector_fp m_fvec; //!< Vector of function values
    bool m_isLinear; //!< Boolean indicating interpolation method
};


//! The Const1 class implements a constant
class Const1 : public Func1
{
public:
    //! Constructor.
    /*!
     * @param A   Constant
     */
    Const1(double A) :
        Func1() {
        m_c = A;
    }

    Const1(const Const1& b) :
        Func1(b) {
    }

    Const1& operator=(const Const1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        return *this;
    }

    virtual std::string write(const std::string& arg) const;
    virtual int ID() const {
        return ConstFuncType;
    }
    virtual doublereal eval(doublereal t) const {
        return m_c;
    }
    virtual Func1& duplicate() const {
        return *(new Const1(m_c));
    }

    virtual Func1& derivative() const {
        Func1* z = new Const1(0.0);
        return *z;
    }
};


/**
 * Sum of two functions.
 */
class Sum1 : public Func1
{
public:
    Sum1(Func1& f1, Func1& f2) :
        Func1() {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    virtual ~Sum1() {
        delete m_f1;
        delete m_f2;
    }

    Sum1(const Sum1& b) :
        Func1(b) {
        *this = Sum1::operator=(b);
    }

    Sum1& operator=(const Sum1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f2 = &m_f2->duplicate();
        m_f1->setParent(this);
        m_f2->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return SumFuncType;
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) + m_f2->eval(t);
    }

    virtual Func1& duplicate() const {
        Func1& f1d = m_f1->duplicate();
        Func1& f2d = m_f2->duplicate();
        return newSumFunction(f1d, f2d);
    }

    virtual Func1& derivative() const {
        Func1& d1 = m_f1->derivative();
        Func1& d2 = m_f2->derivative();
        return newSumFunction(d1, d2);
    }
    virtual int order() const {
        return 0;
    }

    virtual std::string write(const std::string& arg) const;
};


/**
 * Difference of two functions.
 */
class Diff1 : public Func1
{
public:
    Diff1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    virtual ~Diff1() {
        delete m_f1;
        delete m_f2;
    }

    Diff1(const Diff1& b) :
        Func1(b) {
        *this = Diff1::operator=(b);
    }

    Diff1& operator=(const Diff1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f2 = &m_f2->duplicate();
        m_f1->setParent(this);
        m_f2->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return DiffFuncType;
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) - m_f2->eval(t);
    }

    virtual Func1& duplicate() const {
        Func1& f1d = m_f1->duplicate();
        Func1& f2d = m_f2->duplicate();
        return newDiffFunction(f1d, f2d);
    }
    virtual Func1& derivative() const {
        return newDiffFunction(m_f1->derivative(), m_f2->derivative());
    }
    virtual int order() const {
        return 0;
    }

    virtual std::string write(const std::string& arg) const;
};


/**
 * Product of two functions.
 */
class Product1 : public Func1
{
public:
    Product1(Func1& f1, Func1& f2) :
        Func1() {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    virtual ~Product1() {
        delete m_f1;
        delete m_f2;
    }

    Product1(const Product1& b) :
        Func1(b) {
        *this = Product1::operator=(b);
    }

    Product1& operator=(const Product1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f2 = &m_f2->duplicate();
        m_f1->setParent(this);
        m_f2->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return ProdFuncType;
    }

    virtual Func1& duplicate() const {
        Func1& f1d = m_f1->duplicate();
        Func1& f2d = m_f2->duplicate();
        return newProdFunction(f1d, f2d);
    }

    virtual std::string write(const std::string& arg) const;

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) * m_f2->eval(t);
    }

    virtual Func1& derivative() const {
        Func1& a1 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
        Func1& a2 = newProdFunction(m_f2->duplicate(), m_f1->derivative());
        return newSumFunction(a1, a2);
    }
    virtual int order() const {
        return 1;
    }
};

/**
 * Product of two functions.
 */
class TimesConstant1 : public Func1
{
public:
    TimesConstant1(Func1& f1, doublereal A) :
        Func1() {
        m_f1 = &f1;
        m_c = A;
        m_f1->setParent(this);
    }

    virtual ~TimesConstant1() {
        delete m_f1;
    }

    TimesConstant1(const TimesConstant1& b) :
        Func1(b) {
        *this = TimesConstant1::operator=(b);
    }

    TimesConstant1& operator=(const TimesConstant1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f1->setParent(this);
        m_parent = 0;
        return *this;
    }
    virtual int ID() const {
        return TimesConstantFuncType;
    }

    virtual Func1& duplicate() const {
        Func1& f1 = m_f1->duplicate();
        Func1* dup = new TimesConstant1(f1, m_c);
        return *dup;
    }

    virtual doublereal isProportional(TimesConstant1& other) {
        if (func1().isIdentical(other.func1())) {
            return (other.c()/c());
        } else {
            return 0.0;
        }
    }

    virtual doublereal isProportional(Func1& other) {
        if (func1().isIdentical(other)) {
            return 1.0/c();
        } else {
            return 0.0;
        }
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) * m_c;
    }

    virtual Func1& derivative() const {
        Func1& f1d = m_f1->derivative();
        Func1* d = &newTimesConstFunction(f1d, m_c);
        return *d;
    }

    virtual std::string write(const std::string& arg) const;

    virtual int order() const {
        return 0;
    }
};

/**
 * A function plus a constant.
 */
class PlusConstant1 : public Func1
{
public:
    PlusConstant1(Func1& f1, doublereal A) :
        Func1() {
        m_f1 = &f1;
        m_c = A;
        m_f1->setParent(this);
    }

    virtual ~PlusConstant1() {
        delete m_f1;
    }

    PlusConstant1(const PlusConstant1& b) :
        Func1(b) {
        *this = PlusConstant1::operator=(b);
    }

    PlusConstant1& operator=(const PlusConstant1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f1->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return PlusConstantFuncType;
    }

    virtual Func1& duplicate() const {
        Func1& f1 = m_f1->duplicate();
        Func1* dup = new PlusConstant1(f1, m_c);
        return *dup;
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) + m_c;
    }
    virtual Func1& derivative() const {
        return m_f1->derivative();
    }
    virtual std::string write(const std::string& arg) const;

    virtual int order() const {
        return 0;
    }
};


/**
 * Ratio of two functions.
 */
class Ratio1 : public Func1
{
public:
    Ratio1(Func1& f1, Func1& f2) :
        Func1() {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    virtual ~Ratio1() {
        delete m_f1;
        delete m_f2;
    }

    Ratio1(const Ratio1& b) :
        Func1(b) {
        *this = Ratio1::operator=(b);
    }

    Ratio1& operator=(const Ratio1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f2 = &m_f2->duplicate();
        m_f1->setParent(this);
        m_f2->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return RatioFuncType;
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(t) / m_f2->eval(t);
    }

    virtual Func1& duplicate() const {
        Func1& f1d = m_f1->duplicate();
        Func1& f2d = m_f2->duplicate();
        return newRatioFunction(f1d, f2d);
    }

    virtual Func1& derivative() const {
        Func1& a1 = newProdFunction(m_f1->derivative(), m_f2->duplicate());
        Func1& a2 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
        Func1& s = newDiffFunction(a1, a2);
        Func1& p = newProdFunction(m_f2->duplicate(), m_f2->duplicate());
        return newRatioFunction(s, p);
    }

    virtual std::string write(const std::string& arg) const;

    virtual int order() const {
        return 1;
    }
};

/**
 * Composite function.
 */
class Composite1 : public Func1
{
public:
    Composite1(Func1& f1, Func1& f2) :
        Func1() {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    virtual ~Composite1() {
        delete m_f1;
        delete m_f2;
    }

    Composite1(const Composite1& b) :
        Func1(b) {
        *this = Composite1::operator=(b);
    }

    Composite1& operator=(const Composite1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &m_f1->duplicate();
        m_f2 = &m_f2->duplicate();
        m_f1->setParent(this);
        m_f2->setParent(this);
        m_parent = 0;
        return *this;
    }

    virtual int ID() const {
        return CompositeFuncType;
    }

    virtual doublereal eval(doublereal t) const {
        return m_f1->eval(m_f2->eval(t));
    }

    virtual Func1& duplicate() const {
        Func1& f1d = m_f1->duplicate();
        Func1& f2d = m_f2->duplicate();
        return newCompositeFunction(f1d, f2d);
    }

    virtual Func1& derivative() const {
        Func1* d1 = &m_f1->derivative();

        Func1* d3 = &newCompositeFunction(*d1, m_f2->duplicate());
        Func1* d2 = &m_f2->derivative();
        Func1* p = &newProdFunction(*d3, *d2);
        return *p;
    }

    virtual std::string write(const std::string& arg) const;

    virtual int order() const {
        return 2;
    }
};

// The functors below are the old-style ones. They still work,
// but can't do derivatives.

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
class Gaussian : public Func1
{
public:
    Gaussian(double A, double t0, double fwhm) :
        Func1() {
        m_A = A;
        m_t0 = t0;
        m_tau = fwhm/(2.0*std::sqrt(std::log(2.0)));
    }

    Gaussian(const Gaussian& b) :
        Func1(b) {
        *this = Gaussian::operator=(b);
    }

    Gaussian& operator=(const Gaussian& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_A = right.m_A;
        m_t0 = right.m_t0;
        m_tau = right.m_tau;
        m_parent = 0;
        return *this;
    }

    virtual Func1& duplicate() const {
        Gaussian* np = new Gaussian(*this);
        return *((Func1*)np);
    }

    virtual doublereal eval(doublereal t) const {
        doublereal x = (t - m_t0)/m_tau;
        return m_A * std::exp(-x*x);
    }

protected:
    doublereal m_A, m_t0, m_tau;
};


/**
 * Polynomial of degree n.
 */
class Poly1 : public Func1
{
public:
    Poly1(size_t n, const double* c) :
        Func1() {
        m_cpoly.resize(n+1);
        std::copy(c, c+m_cpoly.size(), m_cpoly.begin());
    }

    Poly1(const Poly1& b) :
        Func1(b) {
        *this = Poly1::operator=(b);
    }

    Poly1& operator=(const Poly1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_cpoly = right.m_cpoly;
        m_parent = 0;
        return *this;
    }

    virtual Func1& duplicate() const {
        Poly1* np = new Poly1(*this);
        return *((Func1*)np);
    }

    virtual doublereal eval(doublereal t) const {
        doublereal r = m_cpoly[m_cpoly.size()-1];
        for (size_t n = 1; n < m_cpoly.size(); n++) {
            r *= t;
            r += m_cpoly[m_cpoly.size() - n - 1];
        }
        return r;
    }

protected:
    vector_fp m_cpoly;
};


/**
 * Fourier cosine/sine series.
 *
 * \f[
 * f(t) = \frac{A_0}{2} +
 * \sum_{n=1}^N A_n \cos (n \omega t) + B_n \sin (n \omega t)
 * \f]
 */
class Fourier1 : public Func1
{
public:
    Fourier1(size_t n, double omega, double a0,
             const double* a, const double* b) :
        Func1() {
        m_omega = omega;
        m_a0_2 = 0.5*a0;
        m_ccos.resize(n);
        m_csin.resize(n);
        std::copy(a, a+n, m_ccos.begin());
        std::copy(b, b+n, m_csin.begin());
    }

    Fourier1(const Fourier1& b) :
        Func1(b) {
        *this = Fourier1::operator=(b);
    }

    Fourier1& operator=(const Fourier1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_omega = right.m_omega;
        m_a0_2 = right.m_a0_2;
        m_ccos = right.m_ccos;
        m_csin = right.m_csin;
        m_parent = 0;
        return *this;
    }

    virtual Func1& duplicate() const {
        Fourier1* np = new Fourier1(*this);
        return *((Func1*)np);
    }

    virtual doublereal eval(doublereal t) const {
        size_t n, nn;
        doublereal sum = m_a0_2;
        for (n = 0; n < m_ccos.size(); n++) {
            nn = n + 1;
            sum += m_ccos[n]*std::cos(m_omega*nn*t)
                   + m_csin[n]*std::sin(m_omega*nn*t);
        }
        return sum;
    }

protected:
    doublereal m_omega, m_a0_2;
    vector_fp m_ccos, m_csin;
};


/**
 * Sum of Arrhenius terms.
 * \f[
 * f(T) = \sum_{n=1}^N A_n T^b_n \exp(-E_n/T)
 * \f]
 */
class Arrhenius1 : public Func1
{
public:
    Arrhenius1(size_t n, const double* c) :
        Func1() {
        m_A.resize(n);
        m_b.resize(n);
        m_E.resize(n);
        for (size_t i = 0; i < n; i++) {
            size_t loc = 3*i;
            m_A[i] = c[loc];
            m_b[i] = c[loc+1];
            m_E[i] = c[loc+2];
        }
    }

    Arrhenius1(const Arrhenius1& b) :
        Func1() {
        *this = Arrhenius1::operator=(b);
    }

    Arrhenius1& operator=(const Arrhenius1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_A = right.m_A;
        m_b = right.m_b;
        m_E = right.m_E;
        m_parent = 0;
        return *this;
    }

    virtual Func1& duplicate() const {
        Arrhenius1* np = new Arrhenius1(*this);
        return *((Func1*)np);
    }

    virtual doublereal eval(doublereal t) const {
        doublereal sum = 0.0;
        for (size_t n = 0; n < m_A.size(); n++) {
            sum += m_A[n]*std::pow(t,m_b[n])*std::exp(-m_E[n]/t);
        }
        return sum;
    }

protected:
    vector_fp m_A, m_b, m_E;
};

/**
 *  Periodic function. Takes any function and makes it periodic with period T.
 */
class Periodic1 : public Func1
{
public:
    Periodic1(Func1& f, doublereal T) :
        Func1() {
        m_func = &f;
        m_c = T;
    }

    Periodic1(const Periodic1& b) :
        Func1() {
        *this = Periodic1::operator=(b);
    }

    Periodic1& operator=(const Periodic1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_func = &right.m_func->duplicate();
        return *this;
    }

    virtual Func1& duplicate() const {
        Periodic1* np = new Periodic1(*this);
        return *((Func1*)np);
    }

    virtual ~Periodic1() {
        delete m_func;
    }

    virtual doublereal eval(doublereal t) const {
        int np = int(t/m_c);
        doublereal time = t - np*m_c;
        return m_func->eval(time);
    }

protected:
    Func1* m_func;
};

}

#endif
