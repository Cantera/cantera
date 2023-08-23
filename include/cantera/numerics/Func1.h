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

// Magic numbers are only used by legacy C API methods
// Example: traditional MATLAB toolbox
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

//! @defgroup func1 Functor Objects
//! Functors implement functions of a single variable @f$ f(x) @f$.
//! Functor objects can be combined to form compound expressions, which allows for
//! the implementation of generic mathematical expressions.
//! @ingroup numerics

//! @defgroup func1simple Simple Functors
//! Simple functors implement standard mathematical expressions with a single
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
//! - @c "polynomial" (class Poly1),
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

    Func1(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
        : m_f1_shared(f1), m_f2_shared(f2)
    {
        m_f1 = f1.get();
        m_f2 = f2.get();
    }

    Func1(shared_ptr<Func1> f1, double A) : m_c(A), m_f1_shared(f1) {
        m_f1 = f1.get();
    }

    virtual ~Func1() = default;

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1(const Func1& right);

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1& operator=(const Func1& right);

    //! Duplicate the current function.
    /*!
     * This duplicates the current function, returning a reference to the newly
     * created function.
     * @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
     */
    virtual Func1& duplicate() const;

    //! @deprecated To be removed after %Cantera 3.0. Replaced by type.
    virtual int ID() const;

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
     * This will create a new derivative function and return a reference to the
     * function.
     * @deprecated To be changed after %Cantera 3.0; for new behavior, see derivative3.
     */
    virtual Func1& derivative() const;

    //! Creates a derivative to the current function
    /*!
     * This will create a new derivative function
     * @return  shared pointer to new derivative function.
     * @since New in %Cantera 3.0.
     */
    virtual shared_ptr<Func1> derivative3() const;

    //! Routine to determine if two functions are the same.
    /*!
     * Two functions are the same if they are the same function. This means
     * that the ID and stored constant is the same. This means that the m_f1
     * and m_f2 are identical if they are non-null.
     */
    bool isIdentical(Func1& other) const;

    virtual double isProportional(TimesConstant1& other);
    virtual double isProportional(Func1& other);

    //! Write LaTeX string describing function.
    virtual string write(const string& arg) const;

    //! Accessor function for the stored constant
    double c() const;

    //! Function to set the stored constant
    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    void setC(double c);

    //! accessor function for m_f1
    //! @deprecated To be removed after %Cantera 3.0; replaced by func1_shared().
    Func1& func1() const;

    //! Accessor function for m_f1_shared
    //! @since New in %Cantera 3.0.
    shared_ptr<Func1> func1_shared() const {
        return m_f1_shared;
    }

    //! accessor function for m_f2
    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1& func2() const;

    //! Accessor function for m_f2_shared
    //! @since New in %Cantera 3.0.
    shared_ptr<Func1> func2_shared() const {
        return m_f2_shared;
    }

    //! Return the order of the function, if it makes sense
    virtual int order() const;

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1& func1_dup() const;

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1& func2_dup() const;

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    Func1* parent() const;

    //! @deprecated To be removed after %Cantera 3.0. Only used by deprecated methods.
    void setParent(Func1* p);

protected:
    double m_c = 0.0;
    Func1* m_f1 = nullptr;
    Func1* m_f2 = nullptr;
    Func1* m_parent = nullptr;

    shared_ptr<Func1> m_f1_shared;
    shared_ptr<Func1> m_f2_shared;
};


// all functions using references are deprecated
Func1& newSumFunction(Func1& f1, Func1& f2);
Func1& newDiffFunction(Func1& f1, Func1& f2);
Func1& newProdFunction(Func1& f1, Func1& f2);
Func1& newRatioFunction(Func1& f1, Func1& f2);
Func1& newCompositeFunction(Func1& f1, Func1& f2);
Func1& newTimesConstFunction(Func1& f1, double c);
Func1& newPlusConstFunction(Func1& f1, double c);


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
 * The functor class with type @c "sin" returns @f$ f(x) = \cos(\omega x) @f$,
 * where the argument @f$ x @f$ is in radians.
 * @param omega  Frequency @f$ \omega @f$ (default=1.0)
 * @ingroup func1simple
 */
class Sin1 : public Func1
{
public:
    Sin1(double omega=1.0) {
        m_c = omega;
    }

    //! Constructor uses single parameter (frequency)
    Sin1(const vector<double>& params);

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

    string write(const string& arg) const override;

    int ID() const override {
        return SinFuncType;
    }

    string type() const override {
        return "sin";
    }

    double eval(double t) const override{
        return sin(m_c*t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override;
};


//! Implements the @c cos() function.
/*!
 * The functor class with type @c "cos" returns @f$ f(x) = \cos(\omega x) @f$,
 * where the argument @f$ x @f$ is in radians.
 * @param omega  Frequency @f$ \omega @f$ (default=1.0)
 * @ingroup func1simple
 */
class Cos1 : public Func1
{
public:
    Cos1(double omega=1.0) {
        m_c = omega;
    }

    //! Constructor uses single parameter (frequency)
    Cos1(const vector<double>& params);

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

    Func1& duplicate() const override;
    string write(const string& arg) const override;
    int ID() const override {
        return CosFuncType;
    }
    string type() const override {
        return "cos";
    }

    double eval(double t) const override {
        return cos(m_c * t);
    }
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override;
};


//! Implements the @c exp() (exponential) function.
/*!
 * The functor class with type @c "exp" returns @f$ f(x) = \exp(a x) @f$.
 * @param a  Factor (default=1.0)
 * @ingroup func1simple
 */
class Exp1 : public Func1
{
public:
    Exp1(double a=1.0) {
        m_c = a;
    }

    //! Constructor uses single parameter (exponent factor)
    Exp1(const vector<double>& params);

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
    string write(const string& arg) const override;
    int ID() const override {
        return ExpFuncType;
    }
    string type() const override {
        return "exp";
    }

    double eval(double t) const override {
        return exp(m_c*t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override;
};


//! Implements the @c log() (natural logarithm) function.
/*!
 * The functor class with type @c "log" returns @f$ f(x) = \ln(a x) @f$.
 * @param a  Factor (default=1.0)
 * @ingroup func1simple
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

    shared_ptr<Func1> derivative3() const override;

    string write(const string& arg) const override;
};

//! Implements the @c pow() (power) function.
/*!
 * The functor class with type @c "pow" returns @f$ f(x) = x^n @f$.
 * @param n  Exponent
 * @ingroup func1simple
 */
class Pow1 : public Func1
{
public:
    Pow1(double n) {
        m_c = n;
    }

    //! Constructor uses single parameter (exponent)
    Pow1(const vector<double>& params);

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
    string write(const string& arg) const override;
    int ID() const override {
        return PowFuncType;
    }
    string type() const override {
        return "pow";
    }

    double eval(double t) const override {
        return pow(t, m_c);
    }
    Func1& duplicate() const override;
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override;
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

    string write(const string& arg) const override;
    int ID() const override {
        return TabulatedFuncType;
    }
    string type() const override {
        if (m_isLinear) {
            return "tabulated-linear";
        }
        return "tabulated-previous";
    }

    double eval(double t) const override;
    Func1& duplicate() const override;
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override;
private:
    vector<double> m_tvec; //!< Vector of time values
    vector<double> m_fvec; //!< Vector of function values
    bool m_isLinear; //!< Boolean indicating interpolation method
};


//! Implements a constant.
/*!
 * The functor class with type @c "constant" returns @f$ f(x) = a @f$.
 * @param a  Constant
 * @ingroup func1simple
 */
class Const1 : public Func1
{
public:
    Const1(double a) {
        m_c = a;
    }

    //! Constructor uses single parameter (constant)
    Const1(const vector<double>& params);

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

    string write(const string& arg) const override;
    int ID() const override {
        return ConstFuncType;
    }
    string type() const override {
        return "constant";
    }

    double eval(double t) const override {
        return m_c;
    }
    Func1& duplicate() const override;
    Func1& derivative() const override;
    shared_ptr<Func1> derivative3() const override {
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
    Sum1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    Sum1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    ~Sum1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
        if (!m_f2_shared) {
            delete m_f2;
        }
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

    int ID() const override {
        return SumFuncType;
    }
    string type() const override {
        return "sum";
    }

    double eval(double t) const override {
        return m_f1->eval(t) + m_f2->eval(t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override {
        return newSumFunction(m_f1_shared->derivative3(), m_f2_shared->derivative3());
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
    Diff1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    Diff1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    ~Diff1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
        if (!m_f2_shared) {
            delete m_f2;
        }
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

    int ID() const override {
        return DiffFuncType;
    }

    string type() const override {
        return "diff";
    }

    double eval(double t) const override {
        return m_f1->eval(t) - m_f2->eval(t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override {
        return newDiffFunction(m_f1_shared->derivative3(), m_f2_shared->derivative3());
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
    Product1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    Product1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    ~Product1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
        if (!m_f2_shared) {
            delete m_f2;
        }
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

    int ID() const override {
        return ProdFuncType;
    }

    string type() const override {
        return "product";
    }

    string write(const string& arg) const override;

    double eval(double t) const override {
        return m_f1->eval(t) * m_f2->eval(t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override;

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
    TimesConstant1(Func1& f1, double a) {
        m_f1 = &f1;
        m_c = a;
        m_f1->setParent(this);
    }

    TimesConstant1(shared_ptr<Func1> f1, double a) : Func1(f1, a) {}

    ~TimesConstant1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
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
    int ID() const override {
        return TimesConstantFuncType;
    }
    string type() const override {
        return "times-constant";
    }

    double isProportional(TimesConstant1& other) override {
        if (func1().isIdentical(other.func1())) {
            return (other.c()/c());
        } else {
            return 0.0;
        }
    }

    double isProportional(Func1& other) override {
        if (func1().isIdentical(other)) {
            return 1.0/c();
        } else {
            return 0.0;
        }
    }

    double eval(double t) const override {
        return m_f1->eval(t) * m_c;
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override {
        return newTimesConstFunction(m_f1_shared->derivative3(), m_c);
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
    PlusConstant1(Func1& f1, double a) {
        m_f1 = &f1;
        m_c = a;
        m_f1->setParent(this);
    }

    PlusConstant1(shared_ptr<Func1> f1, double a) : Func1(f1, a) {}

    ~PlusConstant1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
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

    int ID() const override {
        return PlusConstantFuncType;
    }
    string type() const override {
        return "plus-constant";
    }

    double eval(double t) const override {
        return m_f1->eval(t) + m_c;
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override {
        return m_f1_shared->derivative3();
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
    Ratio1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    Ratio1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    ~Ratio1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
        if (!m_f2_shared) {
            delete m_f2;
        }
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

    int ID() const override {
        return RatioFuncType;
    }
    string type() const override {
        return "ratio";
    }

    double eval(double t) const override {
        return m_f1->eval(t) / m_f2->eval(t);
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override;

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
    Composite1(Func1& f1, Func1& f2) {
        m_f1 = &f1;
        m_f2 = &f2;
        m_f1->setParent(this);
        m_f2->setParent(this);
    }

    Composite1(shared_ptr<Func1> f1, shared_ptr<Func1> f2) : Func1(f1, f2) {}

    ~Composite1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
        if (!m_f2_shared) {
            delete m_f2;
        }
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

    int ID() const override {
        return CompositeFuncType;
    }
    string type() const override {
        return "composite";
    }

    double eval(double t) const override {
        return m_f1->eval(m_f2->eval(t));
    }

    Func1& duplicate() const override;
    Func1& derivative() const override;

    shared_ptr<Func1> derivative3() const override;

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

    Gaussian1(const Gaussian1& b) :
        Func1(b) {
        *this = Gaussian1::operator=(b);
    }

    Gaussian1& operator=(const Gaussian1& right) {
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

    string type() const override {
        return "Gaussian";
    }

    double eval(double t) const override {
        double x = (t - m_t0)/m_tau;
        return m_A * std::exp(-x*x);
    }

protected:
    double m_A, m_t0, m_tau;
};


/**
 * A Gaussian.
 * @f[
 * f(t) = A e^{-[(t - t_0)/\tau]^2}
 * @f]
 * where @f[ \tau = \frac{fwhm}{2\sqrt{\ln 2}} @f]
 * @param A peak value
 * @param t0 offset
 * @param fwhm full width at half max
 * @ingroup func1advanced
 * @deprecated To be removed after %Cantera 3.0; replaced by Gaussian1.
 */
class Gaussian : public Gaussian1
{
    Gaussian(double A, double t0, double fwhm);

    Gaussian(const Gaussian& b);

    Func1& duplicate() const override;
};


/**
 * Implements a polynomial of degree @e n.
 * The functor class with type @c "polynomial" returns
 * @f[
 * f(x) = a_n x^n + \dots + a_1 x + a_0
 * @f]
 * @ingroup func1advanced
 */
class Poly1 : public Func1
{
public:
    Poly1(size_t n, const double* c) {
        m_cpoly.resize(n+1);
        std::copy(c, c+m_cpoly.size(), m_cpoly.begin());
    }

    //! Constructor uses @f$ n + 1 @f$ parameters in the following order:
    //! @f$ [a_n, \dots, a_1, a_0] @f$
    Poly1(const vector<double>& params);

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

    string type() const override {
        return "polynomial";
    }

    Func1& duplicate() const override;

    double eval(double t) const override {
        double r = m_cpoly[m_cpoly.size()-1];
        for (size_t n = 1; n < m_cpoly.size(); n++) {
            r *= t;
            r += m_cpoly[m_cpoly.size() - n - 1];
        }
        return r;
    }

protected:
    vector<double> m_cpoly;
};


/**
 * Implements a Fourier cosine/sine series.
 * The functor class with type @c "Fourier" returns
 * @f[
 * f(t) = \frac{A_0}{2} +
 * \sum_{n=1}^N A_n \cos (n \omega t) + B_n \sin (n \omega t)
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

    string type() const override {
        return "Fourier";
    }

    Func1& duplicate() const override;

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

    string type() const override {
        return "Arrhenius";
    }

    Func1& duplicate() const override;

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
    Periodic1(Func1& f, double T) {
        m_f1 = &f;
        m_c = T;
    }

    Periodic1(const Periodic1& b) {
        *this = Periodic1::operator=(b);
    }

    Periodic1(shared_ptr<Func1> f, double A) : Func1(f, A) {}

    Periodic1& operator=(const Periodic1& right) {
        if (&right == this) {
            return *this;
        }
        Func1::operator=(right);
        m_f1 = &right.m_f1->duplicate();
        return *this;
    }

    string type() const override {
        return "periodic";
    }

    Func1& duplicate() const override;

    ~Periodic1() override {
        if (!m_f1_shared) {
            delete m_f1;
        }
    }

    double eval(double t) const override {
        int np = int(t/m_c);
        double time = t - np*m_c;
        return m_f1->eval(time);
    }
};

}

#endif
