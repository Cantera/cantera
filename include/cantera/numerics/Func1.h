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

class TimesConstant1;

//! @defgroup func1 Functor Objects
//! Functors implement functions of a single variable @f$ f(x) @f$.
//! Functor objects can be combined to form compound expressions, which allows for
//! the implementation of generic mathematical expressions.
//! @ingroup numerics

//! @defgroup func1basic Basic Functors
//! Basic functors implement standard mathematical expressions with a single
//! parameter.
//! The following simple functor types are implemented:
//! - @c "sin" (class Sin1), @c "cos" (class Cos1),
//! - @c "exp" (class Exp1), @c "log" (class Log1),
//! - @c "pow" (class Pow1),
//! - @c "constant" (class Const1).
//! @ingroup func1

//! @defgroup func1advanced Advanced Functors
//! Advanced functors implement expressions that require multiple parameters.
//! The following advanced functor types are implemented:
//! - @c "tabulated-linear" and @c "tabulated-previous" (class Tabulated1),
//! - @c "polynomial3" (class Poly13),
//! - @c "Fourier" (class Fourier1),
//! - @c "Gaussian" (class Gaussian1),
//! - @c "Arrhenius" (class Arrhenius1).
//! @ingroup func1

//! @defgroup func1compound Compound Functors
//! Compound functors implement expressions that are composed of other functors.
//! The following compound functor types are implemented:
//! - @c "sum" (class Sum1),
//! - @c "diff" (class Diff1),
//! - @c "product" (class Product1),
//! - @c "ratio" (class Ratio1),
//! - @c "composite" (class Composite1),
//! @ingroup func1

//! @defgroup func1modified Modified Functors
//! Modified functors implement expressions that involve one functor and
//! a single parameter.
//! The following modified functor types are implemented:
//! - @c "times-constant" (class TimesConstant1),
//! - @c "plus-constant" (class PlusConstant1),
//! - @c "periodic" (class Periodic1).
//! @ingroup func1

//! @defgroup func1helper Helper Functions
//! Helper functions detect simplifications that can be made to compound expressions.
//! @ingroup func1

/**
 * Base class for 'functor' classes that evaluate a function of one variable.
 * @ingroup func1
 */
class Func1
{
public:
    Func1() = default;

    Func1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : m_f1(f1), m_f2(f2) {}

    Func1(shared_ptr<Func1> f1, double A) : m_c(A), m_f1(f1) {}

    virtual ~Func1() = default;

    // Func1(const Func1& right) = delete;  //! @todo Uncomment after %Cantera 3.1
    Func1& operator=(const Func1& right) = delete;

    //! Returns a string describing the type of the function
    //! @since New in %Cantera 3.0.
    virtual string type() const {
        return "functor";
    }

    //! Returns a string with the class name of the functor
    //! @since New in %Cantera 3.0.
    string typeName() const;

    //! Calls method eval to evaluate the function
    double operator()(double t) const;

    //! Evaluate the function.
    virtual double eval(double t) const;

    //! Creates a derivative to the current function
    /*!
     * @return  shared pointer to new derivative function.
     * @since Starting in %Cantera 3.1, the return type is a `shared_ptr`.
     */
    virtual shared_ptr<Func1> derivative() const;

    //! Routine to determine if two functions are the same.
    /*!
     * Two functions are the same if they are the same function. For example, either
     * ID and stored constant are the same, or the #m_f1 and #m_f2 are identical if they
     * are non-null. Functors of the base class Func1 are by default not identical, as
     * they are used by callback functions that cannot be differentiated. In instances
     * where exact comparisons are not implemented, `false` is returned to prevent false
     * positives that could lead to incorrect simplifications of compound functors.
     */
    virtual bool isIdentical(shared_ptr<Func1> other) const;

    //! Routine to determine if two functions are the same.
    /*!
     * @deprecated Deprecated in %Cantera 3.1 and removed thereafter; replaced by
     *      isIdentical(shared_ptr<Func1>&).
     * @todo Restore deleted copy constructor after removal.
     */
    virtual bool isIdentical(Func1& other) const;

    /**
     * @deprecated Deprecated in %Cantera 3.1 and removed thereafter; replaced by
     *      internal function.
     */
    virtual double isProportional(TimesConstant1& other);
    /**
     * @deprecated Deprecated in %Cantera 3.1 and removed thereafter; replaced by
     *      internal function.
     */
    virtual double isProportional(Func1& other);

    //! Write LaTeX string describing function.
    virtual string write(const string& arg) const;

    //! Accessor function for the stored constant #m_c.
    double c() const;

    //! Accessor function for #m_f1.
    //! @since New in %Cantera 3.0.
    shared_ptr<Func1> func1_shared() const {
        return m_f1;
    }

    //! Accessor function for #m_f2.
    //! @since New in %Cantera 3.0.
    shared_ptr<Func1> func2_shared() const {
        return m_f2;
    }

    //! Return the order of the function, if it makes sense
    virtual int order() const;

protected:
    double m_c = 0.0;
    shared_ptr<Func1> m_f1;
    shared_ptr<Func1> m_f2;
};

//! Sum of two functions.
//! @ingroup func1helper
shared_ptr<Func1> newSumFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2);

//! Difference of two functions.
//! @ingroup func1helper
shared_ptr<Func1> newDiffFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2);

//! Product of two functions.
//! @ingroup func1helper
shared_ptr<Func1> newProdFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2);

//! Ratio of two functions.
//! @ingroup func1helper
shared_ptr<Func1> newRatioFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2);

//! Composite of two functions.
//! @ingroup func1helper
shared_ptr<Func1> newCompositeFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2);

//! Product of function and constant.
//! @ingroup func1helper
shared_ptr<Func1> newTimesConstFunction(shared_ptr<Func1> f1, double c);

//! Sum of function and constant.
//! @ingroup func1helper
shared_ptr<Func1> newPlusConstFunction(shared_ptr<Func1> f1, double c);

//! Implements the @c sin() function.
/*!
 * The functor class with type @c "sin" returns @f$ f(x) = \sin(\omega x) @f$,
 * where the argument @f$ x @f$ is in radians.
 * @param omega  Frequency @f$ \omega @f$ (default=1.0)
 * @ingroup func1basic
 */
class Sin1 : public Func1
{
public:
    Sin1(double omega=1.0) {
        m_c = omega;
    }

    //! Constructor uses single parameter (frequency)
    Sin1(const vector<double>& params);

    string write(const string& arg) const override;

    string type() const override {
        return "sin";
    }

    double eval(double t) const override{
        return sin(m_c*t);
    }

    shared_ptr<Func1> derivative() const override;
};


//! Implements the @c cos() function.
/*!
 * The functor class with type @c "cos" returns @f$ f(x) = \cos(\omega x) @f$,
 * where the argument @f$ x @f$ is in radians.
 * @param omega  Frequency @f$ \omega @f$ (default=1.0)
 * @ingroup func1basic
 */
class Cos1 : public Func1
{
public:
    Cos1(double omega=1.0) {
        m_c = omega;
    }

    //! Constructor uses single parameter (frequency)
    Cos1(const vector<double>& params);

    string write(const string& arg) const override;

    string type() const override {
        return "cos";
    }

    double eval(double t) const override {
        return cos(m_c * t);
    }
    shared_ptr<Func1> derivative() const override;
};


//! Implements the @c exp() (exponential) function.
/*!
 * The functor class with type @c "exp" returns @f$ f(x) = \exp(a x) @f$.
 * @param a  Factor (default=1.0)
 * @ingroup func1basic
 */
class Exp1 : public Func1
{
public:
    Exp1(double a=1.0) {
        m_c = a;
    }

    //! Constructor uses single parameter (exponent factor)
    Exp1(const vector<double>& params);

    string write(const string& arg) const override;

    string type() const override {
        return "exp";
    }

    double eval(double t) const override {
        return exp(m_c*t);
    }

    shared_ptr<Func1> derivative() const override;
};


//! Implements the @c log() (natural logarithm) function.
/*!
 * The functor class with type @c "log" returns @f$ f(x) = \ln(a x) @f$.
 * @param a  Factor (default=1.0)
 * @ingroup func1basic
 * @since New in %Cantera 3.0
 */
class Log1 : public Func1
{
public:
    Log1(double a=1.0) {
        m_c = a;
    }

    //! Constructor uses single parameter (factor)
    Log1(const vector<double>& params);

    string type() const override {
        return "log";
    }

    double eval(double t) const override {
        return log(m_c * t);
    }

    shared_ptr<Func1> derivative() const override;

    string write(const string& arg) const override;
};

//! Implements the @c pow() (power) function.
/*!
 * The functor class with type @c "pow" returns @f$ f(x) = x^n @f$.
 * @param n  Exponent
 * @ingroup func1basic
 */
class Pow1 : public Func1
{
public:
    Pow1(double n) {
        m_c = n;
    }

    //! Constructor uses single parameter (exponent)
    Pow1(const vector<double>& params);

    string write(const string& arg) const override;

    string type() const override {
        return "pow";
    }

    double eval(double t) const override {
        return pow(t, m_c);
    }
    shared_ptr<Func1> derivative() const override;
};

//! Implements a tabulated function.
/*!
 * The functor class is based on tabulated arrays @c tvals and @c fvals, where
 * @c tvals contain independent variables and @c fvals are corresponding function
 * values. Depending on configuration, the function is either interpolated linearly
 * between the tabulated points (type @c "tabulated-linear" ; default), or yields
 * the last tabulated value until a new tabulated time value is reached (type
 * @c "tabulated-previous" ).
 * @ingroup func1advanced
 */
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
               const string& method="linear");

    //! Constructor uses @f$ 2 n @f$ parameters in the following order:
    //! @f$ [t_0, t_1, \dots, t_{n-1}, f_0, f_1, \dots, f_{n-1}] @f$
    Tabulated1(const vector<double>& params);

    //! Set the interpolation method
    //! @param method  Evaluation method. If @c "linear" (default), a linear
    //!     interpolation between tabulated values is used; if @c "previous", the
    //!     last tabulated value is held until a new tabulated time value is reached.
    //! @since New in %Cantera 3.0
    void setMethod(const string& method);

    bool isIdentical(shared_ptr<Func1> other) const override {
        return false;  // base class check is insufficient
    }

    string write(const string& arg) const override;

    string type() const override {
        if (m_isLinear) {
            return "tabulated-linear";
        }
        return "tabulated-previous";
    }

    double eval(double t) const override;
    shared_ptr<Func1> derivative() const override;
private:
    vector<double> m_tvec; //!< Vector of time values
    vector<double> m_fvec; //!< Vector of function values
    bool m_isLinear; //!< Boolean indicating interpolation method
};


//! Implements a constant.
/*!
 * The functor class with type @c "constant" returns @f$ f(x) = a @f$.
 * @param a  Constant
 * @ingroup func1basic
 */
class Const1 : public Func1
{
public:
    Const1(double a) {
        m_c = a;
    }

    //! Constructor uses single parameter (constant)
    Const1(const vector<double>& params);

    string write(const string& arg) const override;

    string type() const override {
        return "constant";
    }

    double eval(double t) const override {
        return m_c;
    }
    shared_ptr<Func1> derivative() const override {
        return make_shared<Const1>(0.0);
    }
};


/**
 * Implements the sum of two functions.
 * The functor class with type @c "sum" returns @f$ f(x) = f_1(x) + f_2(x) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param f2  Functor @f$ f_2(x) @f$
 * @ingroup func1compound
 */
class Sum1 : public Func1
{
public:
    Sum1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    string type() const override {
        return "sum";
    }

    double eval(double t) const override {
        return m_f1->eval(t) + m_f2->eval(t);
    }

    shared_ptr<Func1> derivative() const override {
        return newSumFunction(m_f1->derivative(), m_f2->derivative());
    }

    int order() const override {
        return 0;
    }

    string write(const string& arg) const override;
};

/**
 * Implements the difference of two functions.
 * The functor class with type @c "diff" returns @f$ f(x) = f_1(x) - f_2(x) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param f2  Functor @f$ f_2(x) @f$
 * @ingroup func1compound
 */
class Diff1 : public Func1
{
public:
    Diff1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    string type() const override {
        return "diff";
    }

    double eval(double t) const override {
        return m_f1->eval(t) - m_f2->eval(t);
    }

    shared_ptr<Func1> derivative() const override {
        return newDiffFunction(m_f1->derivative(), m_f2->derivative());
    }

    int order() const override {
        return 0;
    }

    string write(const string& arg) const override;
};


/**
 * Implements the product of two functions.
 * The functor class with type @c "product" returns @f$ f(x) = f_1(x) f_2(x) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param f2  Functor @f$ f_2(x) @f$
 * @ingroup func1compound
 */
class Product1 : public Func1
{
public:
    Product1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    string type() const override {
        return "product";
    }

    string write(const string& arg) const override;

    double eval(double t) const override {
        return m_f1->eval(t) * m_f2->eval(t);
    }

    shared_ptr<Func1> derivative() const override;

    int order() const override {
        return 1;
    }
};

/**
 * Implements the product of a function and a constant.
 * The functor class with type @c "times-constant" returns @f$ f(x) = a f_1(x) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param a   Constant @f$ a @f$
 * @ingroup func1modified
 */
class TimesConstant1 : public Func1
{
public:
    TimesConstant1(shared_ptr<Func1> f1, double a) : Func1(f1, a) {}

    string type() const override {
        return "times-constant";
    }

    double isProportional(TimesConstant1& other) override {
        if (func1_shared()->isIdentical(*other.func1_shared())) {
            return (other.c()/c());
        } else {
            return 0.0;
        }
    }

    double isProportional(Func1& other) override {
        if (func1_shared()->isIdentical(other)) {
            return 1.0/c();
        } else {
            return 0.0;
        }
    }

    double eval(double t) const override {
        return m_f1->eval(t) * m_c;
    }

    shared_ptr<Func1> derivative() const override {
        return newTimesConstFunction(m_f1->derivative(), m_c);
    }

    string write(const string& arg) const override;

    int order() const override {
        return 0;
    }
};

/**
 * Implements the sum of a function and a constant.
 * The functor class with type @c "plus-constant" returns @f$ f(x) = f_1(x) + a @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param a   Constant @f$ a @f$
 * @ingroup func1modified
 */
class PlusConstant1 : public Func1
{
public:
    PlusConstant1(shared_ptr<Func1> f1, double a) : Func1(f1, a) {}

    string type() const override {
        return "plus-constant";
    }

    double eval(double t) const override {
        return m_f1->eval(t) + m_c;
    }

    shared_ptr<Func1> derivative() const override {
        return m_f1->derivative();
    }

    string write(const string& arg) const override;

    int order() const override {
        return 0;
    }
};


/**
 * Implements the ratio of two functions.
 * The functor class with type @c "ratio" returns @f$ f(x) = f_1(x) / f_2(x) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param f2  Functor @f$ f_2(x) @f$
 * @ingroup func1compound
 */
class Ratio1 : public Func1
{
public:
    Ratio1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    string type() const override {
        return "ratio";
    }

    double eval(double t) const override {
        return m_f1->eval(t) / m_f2->eval(t);
    }

    shared_ptr<Func1> derivative() const override;

    string write(const string& arg) const override;

    int order() const override {
        return 1;
    }
};

/**
 * Implements a composite function.
 * The functor class with type @c "composite" returns @f$ f(x) = f_1\left(f_2(x)\right) @f$.
 * @param f1  Functor @f$ f_1(x) @f$
 * @param f2  Functor @f$ f_2(x) @f$
 * @ingroup func1compound
 */
class Composite1 : public Func1
{
public:
    Composite1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    string type() const override {
        return "composite";
    }

    double eval(double t) const override {
        return m_f1->eval(m_f2->eval(t));
    }

    shared_ptr<Func1> derivative() const override;

    string write(const string& arg) const override;

    int order() const override {
        return 2;
    }
};

// The functors below are the old-style ones. They still work,
// but can't do derivatives.

/**
 * Implements a Gaussian function.
 * The functor class with type @c "Gaussian" returns
 * @f[
 * f(t) = A e^{-[(t - t_0)/\tau]^2}
 * @f]
 * where @f$ \tau = \mathrm{fwhm} / (2 \sqrt{\ln 2}) @f$.
 * @param A peak value
 * @param t0 offset
 * @param fwhm full width at half max
 * @ingroup func1advanced
 * @since New in %Cantera 3.0.
 */
class Gaussian1 : public Func1
{
public:
    Gaussian1(double A, double t0, double fwhm) {
        m_A = A;
        m_t0 = t0;
        m_tau = fwhm/(2.0*std::sqrt(std::log(2.0)));
    }

    //! Constructor uses 3 parameters in the following order:
    //! @f$ [A, t_0, \mathrm{fwhm}] @f$
    Gaussian1(const vector<double>& params);

    string type() const override {
        return "Gaussian";
    }

    bool isIdentical(shared_ptr<Func1> other) const override {
        return false;  // base class check is insufficient
    }

    double eval(double t) const override {
        double x = (t - m_t0)/m_tau;
        return m_A * std::exp(-x*x);
    }

protected:
    double m_A, m_t0, m_tau;
};


/**
 * Implements a polynomial of degree @e n.
 * The functor class with type @c "polynomial3" returns
 * @f[
 * f(x) = a_n x^n + \dots + a_1 x + a_0
 * @f]
 * with coefficients provided in the order @f$ [a_n, \dots, a_1, a_0] @f$ (consistent
 * with MATLAB and NumPy conventions). Note that %Cantera 3.1 reversed the coefficient
 * order with respect to its earlier definition. A deprecation cycle is skipped as the
 * functor class is not expected to be widely used; the transitional name Poly13 ensures
 * that changed behavior does not go unnoticed. The class name will revert to @b Poly1
 * after %Cantera 3.1.
 * @since Changed in %Cantera 3.1.
 * @todo Rename to Poly1 after %Cantera 3.1
 * @ingroup func1advanced
 */
class Poly13 : public Func1
{
public:
    Poly13(size_t n, const double* c) {
        m_cpoly.resize(n+1);
        std::copy(c, c+m_cpoly.size(), m_cpoly.begin());
    }

    //! Constructor uses @f$ n + 1 @f$ parameters in the following order:
    //! @f$ [a_n, \dots, a_1, a_0] @f$
    Poly13(const vector<double>& params);

    string type() const override {
        return "polynomial3";
    }

    bool isIdentical(shared_ptr<Func1> other) const override {
        return false;  // base class check is insufficient
    }

    double eval(double t) const override {
        double r = m_cpoly[0];
        for (size_t n = 1; n < m_cpoly.size(); n++) {
            r *= t;
            r += m_cpoly[n];
        }
        return r;
    }

    string write(const string& arg) const override;

protected:
    vector<double> m_cpoly;
};


/**
 * Implements a Fourier cosine/sine series.
 * The functor class with type @c "Fourier" returns
 * @f[
 * f(t) = \frac{a_0}{2} + \sum_{n=1}^N a_n \cos (n \omega t) + b_n \sin (n \omega t)
 * @f]
 * @ingroup func1advanced
 */
class Fourier1 : public Func1
{
public:
    Fourier1(size_t n, double omega, double a0, const double* a, const double* b) {
        m_omega = omega;
        m_a0_2 = 0.5*a0;
        m_ccos.resize(n);
        m_csin.resize(n);
        std::copy(a, a+n, m_ccos.begin());
        std::copy(b, b+n, m_csin.begin());
    }

    //! Constructor uses @f$ 2 n + 2 @f$ parameters in the following order:
    //! @f$ [a_0, a_1, \dots, a_n, \omega, b_1, \dots, b_n] @f$
    Fourier1(const vector<double>& params);

    string type() const override {
        return "Fourier";
    }

    bool isIdentical(shared_ptr<Func1> other) const override {
        return false;  // base class check is insufficient
    }

    double eval(double t) const override {
        size_t n, nn;
        double sum = m_a0_2;
        for (n = 0; n < m_ccos.size(); n++) {
            nn = n + 1;
            sum += m_ccos[n]*std::cos(m_omega*nn*t)
                   + m_csin[n]*std::sin(m_omega*nn*t);
        }
        return sum;
    }

protected:
    double m_omega, m_a0_2;
    vector<double> m_ccos, m_csin;
};


/**
 * Implements a sum of Arrhenius terms.
 * The functor class with type @c "Arrhenius" returns
 * @f[
 * f(T) = \sum_{n=1}^N A_n T^b_n \exp(-E_n/T)
 * @f]
 * @ingroup func1advanced
 */
class Arrhenius1 : public Func1
{
public:
    Arrhenius1(size_t n, const double* c) {
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

    //! Constructor uses @f$ 3 n @f$ parameters in the following order:
    //! @f$ [A_1, b_1, E_1, A_2, b_2, E_2, \dots, A_n, b_n, E_n] @f$
    Arrhenius1(const vector<double>& params);

    string type() const override {
        return "Arrhenius";
    }

    bool isIdentical(shared_ptr<Func1> other) const override {
        return false;  // base class check is insufficient
    }

    double eval(double t) const override {
        double sum = 0.0;
        for (size_t n = 0; n < m_A.size(); n++) {
            sum += m_A[n]*std::pow(t,m_b[n])*std::exp(-m_E[n]/t);
        }
        return sum;
    }

protected:
    vector<double> m_A, m_b, m_E;
};

/**
 * Implements a periodic function.
 * Takes any function and makes it periodic with period @f$ T @f$.
 * @param f  Functor to be made periodic
 * @param T  Period
 * @ingroup func1modified
 */
class Periodic1 : public Func1
{
public:
    Periodic1(shared_ptr<Func1> f, double A) : Func1(f, A) {}

    string type() const override {
        return "periodic";
    }

    double eval(double t) const override {
        int np = int(t/m_c);
        double time = t - np*m_c;
        return m_f1->eval(time);
    }
};

}

#endif
