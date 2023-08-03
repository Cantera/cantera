//! @file Func1.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/Func1.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

Func1::Func1(const Func1& right) :
    m_c(right.m_c),
    m_f1(right.m_f1),
    m_f2(right.m_f2),
    m_parent(right.m_parent)
{
}

Func1& Func1::operator=(const Func1& right)
{
    if (&right == this) {
        return *this;
    }
    m_c = right.m_c;
    m_f1 = right.m_f1;
    m_f2 = right.m_f2;
    m_parent = right.m_parent;
    return *this;
}

Func1& Func1::duplicate() const
{
    warn_deprecated("Func1::duplicate",
        "To be removed after Cantera 3.0. No longer needed.");
    Func1* nfunc = new Func1(*this);
    return *nfunc;
}

int Func1::ID() const
{
    return 0;
}

string Func1::typeName() const
{
    return demangle(typeid(*this));
}

// Calls method eval to evaluate the function
double Func1::operator()(double t) const
{
    return eval(t);
}

// Evaluate the function.
double Func1::eval(double t) const
{
    return 0.0;
}

Func1& Func1::derivative() const
{
    warn_deprecated("Func1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    cout << "derivative error... ERR: ID = " << ID() << endl;
    cout << write("x") << endl;
    return *(new Func1);
}

shared_ptr<Func1> Func1::derivative3() const
{
    throw CanteraError("Func1::derivative3",
        "Needs to be overloaded by Func1 specialization.");
}

bool Func1::isIdentical(Func1& other) const
{
    if (ID() != other.ID() || m_c != other.m_c) {
        return false;
    }
    if (m_f1) {
        if (!other.m_f1) {
            return false;
        }
        if (!m_f1->isIdentical(*other.m_f1)) {
            return false;
        }
    }
    if (m_f2) {
        if (!other.m_f2) {
            return false;
        }
        if (!m_f2->isIdentical(*other.m_f2)) {
            return false;
        }
    }
    return true;
}

//! accessor function for the returned constant
double Func1::c() const
{
    return m_c;
}

// Function to set the stored constant
void Func1::setC(double c)
{
    m_c = c;
}

//! accessor function for m_f1
Func1& Func1::func1() const
{
    return *m_f1;
}

Func1& Func1::func2() const
{
    return *m_f2;
}

int Func1::order() const
{
    return 3;
}

Func1& Func1::func1_dup() const
{
    return m_f1->duplicate();
}

Func1& Func1::func2_dup() const
{
    return m_f2->duplicate();
}

Func1* Func1::parent() const
{
    return m_parent;
}

void Func1::setParent(Func1* p)
{
    m_parent = p;
}

/*****************************************************************************/

Sin1::Sin1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Sin1::Sin1",
            "Constructor needs exactly one parameter (frequency).");
    }
    m_c = params[0];
}

string Sin1::write(const string& arg) const
{
    if (m_c == 1.0) {
        return fmt::format("\\sin({})", arg);
    } else {
        return fmt::format("\\sin({}{})", m_c, arg);
    }
}

Func1& Sin1::duplicate() const {
    warn_deprecated("Sin1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Sin1* nfunc = new Sin1(*this);
    return (Func1&) *nfunc;
}

Func1& Sin1::derivative() const
{
    warn_deprecated("Sin1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* c = new Cos1(m_c);
    Func1* r = &newTimesConstFunction(*c, m_c);
    return *r;
}

shared_ptr<Func1> Sin1::derivative3() const
{
    auto c = make_shared<Cos1>(m_c);
    return newTimesConstFunction(c, m_c);
}

/*****************************************************************************/

Cos1::Cos1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Cos1::Cos1",
            "Constructor needs exactly one parameter (frequency).");
    }
    m_c = params[0];
}

Func1& Cos1::duplicate() const {
    warn_deprecated("Cos1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Cos1* nfunc = new Cos1(*this);
    return (Func1&) *nfunc;
}

Func1& Cos1::derivative() const
{
    warn_deprecated("Cos1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* s = new Sin1(m_c);
    Func1* r = &newTimesConstFunction(*s, -m_c);
    return *r;
}

shared_ptr<Func1> Cos1::derivative3() const
{
    auto s = make_shared<Sin1>(m_c);
    return newTimesConstFunction(s, -m_c);
}

string Cos1::write(const string& arg) const
{
    if (m_c == 1.0) {
        return fmt::format("\\cos({})", arg);
    } else {
        return fmt::format("\\cos({}{})", m_c, arg);
    }
}

/**************************************************************************/

Exp1::Exp1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Exp1::Exp1",
            "Constructor needs exactly one parameter (exponent factor).");
    }
    m_c = params[0];
}

Func1& Exp1::duplicate() const {
    warn_deprecated("Exp1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    return *(new Exp1(m_c));
}

Func1& Exp1::derivative() const
{
    warn_deprecated("Exp1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* f = new Exp1(m_c);
    if (m_c != 1.0) {
        return newTimesConstFunction(*f, m_c);
    } else {
        return *f;
    }
}

shared_ptr<Func1> Exp1::derivative3() const
{
    auto f = make_shared<Exp1>(m_c);
    if (m_c != 1.0) {
        return newTimesConstFunction(f, m_c);
    }
    return f;
}

string Exp1::write(const string& arg) const
{
    if (m_c == 1.0) {
        return fmt::format("\\exp({})", arg);
    } else {
        return fmt::format("\\exp({}{})", m_c, arg);
    }
}

Log1::Log1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Log1::Log1",
            "Constructor needs exactly one parameter (factor).");
    }
    m_c = params[0];
}

shared_ptr<Func1> Log1::derivative3() const
{
    auto f = make_shared<Pow1>(-1.);
    if (m_c != 1.0) {
        return newTimesConstFunction(f, m_c);
    }
    return f;
}

string Log1::write(const string& arg) const
{
    if (m_c == 1.0) {
        return fmt::format("\\log({})", arg);
    }
    return fmt::format("\\log({}{})", m_c, arg);
}

/******************************************************************************/

Pow1::Pow1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Pow1::Pow1",
            "Constructor needs exactly one parameter (exponent).");
    }
    m_c = params[0];
}

Func1& Pow1::duplicate() const {
    warn_deprecated("Pow1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    return *(new Pow1(m_c));
}

Func1& Pow1::derivative() const
{
    warn_deprecated("Pow1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* r;
    if (m_c == 0.0) {
        r = new Const1(0.0);
    } else if (m_c == 1.0) {
        r = new Const1(1.0);
    } else {
        Func1* f = new Pow1(m_c - 1.0);
        r = &newTimesConstFunction(*f, m_c);
    }
    return *r;
}

shared_ptr<Func1> Pow1::derivative3() const
{
    if (m_c == 0.0) {
        return make_shared<Const1>(0.0);
    }
    if (m_c == 1.0) {
        return make_shared<Const1>(1.0);
    }
    auto f = make_shared<Pow1>(m_c - 1.);
    return newTimesConstFunction(f, m_c);
}

/******************************************************************************/

Const1::Const1(const vector<double>& params)
{
    if (params.size() != 1) {
        throw CanteraError("Const1::Const1",
            "Constructor needs exactly one parameter (constant).");
    }
    m_c = params[0];
}

Poly1::Poly1(const vector<double>& params)
{
    if (params.size() == 0) {
        throw CanteraError("Poly1::Poly1",
            "Constructor needs an array that is not empty.");
    }
    size_t n = params.size() - 1;
    m_cpoly.resize(n + 1);
    copy(params.data(), params.data() + m_cpoly.size(), m_cpoly.begin());
}

Fourier1::Fourier1(const vector<double>& params)
{
    if (params.size() < 4) {
        throw CanteraError("Fourier1::Fourier1",
            "Constructor needs an array with at least 4 entries.");
    }
    if (params.size() % 2 != 0) {
        throw CanteraError("Fourier1::Fourier1",
            "Constructor needs an array with an even number of entries.");
    }
    size_t n = params.size() / 2 - 1;
    m_omega = params[n + 1];
    m_a0_2 = 0.5 * params[0];
    m_ccos.resize(n);
    m_csin.resize(n);
    copy(params.data() + 1, params.data() + n + 1, m_ccos.begin());
    copy(params.data() + n + 2, params.data() + 2 * n + 2, m_csin.begin());
}

Gaussian1::Gaussian1(const vector<double>& params)
{
    if (params.size() != 3) {
        throw CanteraError("Gaussian1::Gaussian1",
            "Constructor needs exactly 3 parameters (amplitude, center, width).");
    }
    m_A = params[0];
    m_t0 = params[1];
    m_tau = params[2] / (2. * sqrt(log(2.)));
}

Arrhenius1::Arrhenius1(const vector<double>& params)
{
    if (params.size() < 3) {
        throw CanteraError("Arrhenius1::Arrhenius1",
            "Constructor needs an array with at least 3 entries.");
    }
    if (params.size() % 3 != 0) {
        throw CanteraError("Arrhenius1::Arrhenius1",
            "Constructor needs an array with multiples of 3 entries.");
    }
    size_t n = params.size() / 3;
    m_A.resize(n);
    m_b.resize(n);
    m_E.resize(n);
    for (size_t i = 0; i < n; i++) {
        size_t loc = 3 * i;
        m_A[i] = params[loc];
        m_b[i] = params[loc + 1];
        m_E[i] = params[loc + 2];
    }
}

Tabulated1::Tabulated1(size_t n, const double* tvals, const double* fvals,
                       const string& method)
{
    m_tvec.resize(n);
    std::copy(tvals, tvals + n, m_tvec.begin());
    for (auto it = std::begin(m_tvec) + 1; it != std::end(m_tvec); it++) {
        if (*(it - 1) > *it) {
            throw CanteraError("Tabulated1::Tabulated1",
                               "time values are not increasing monotonically.");
        }
    }
    m_fvec.resize(n);
    std::copy(fvals, fvals + n, m_fvec.begin());
    setMethod(method);
}

Tabulated1::Tabulated1(const vector<double>& params) : m_isLinear(true)
{
    if (params.size() < 4) {
        throw CanteraError("Tabulated1::Tabulated1",
            "Constructor needs an array with at least 4 entries.");
    }
    if (params.size() % 2 != 0) {
        throw CanteraError("Tabulated1::Tabulated1",
            "Constructor needs an array with an even number of entries.");
    }
    size_t n = params.size() / 2;
    m_tvec.resize(n);
    copy(params.data(), params.data() + n, m_tvec.begin());
    for (auto it = std::begin(m_tvec) + 1; it != std::end(m_tvec); it++) {
        if (*(it - 1) > *it) {
            throw CanteraError("Tabulated1::Tabulated1",
                "Time values are not monotonically increasing.");
        }
    }
    m_fvec.resize(n);
    copy(params.data() + n, params.data() + 2 * n, m_fvec.begin());
}

void Tabulated1::setMethod(const string& method)
{
    if (method == "linear") {
        m_isLinear = true;
    } else if (method == "previous") {
        m_isLinear = false;
    } else {
        throw NotImplementedError("Tabulated1::setMethod",
            "Interpolation method '{}' is not implemented.", method);
    }
}

double Tabulated1::eval(double t) const {
    size_t siz = m_tvec.size();
    // constructor ensures that siz > 0
    if (t <= m_tvec[0]) {
        return m_fvec[0];
    } else if (t >= m_tvec[siz-1]) {
        return m_fvec[siz-1];
    } else {
        size_t ix = 0;
        while (t > m_tvec[ix+1]) {
            ix++;
        }
        if (m_isLinear) {
            double df = m_fvec[ix+1] - m_fvec[ix];
            df /= m_tvec[ix+1] - m_tvec[ix];
            df *= t - m_tvec[ix];
            return m_fvec[ix] + df;
        } else {
            return m_fvec[ix];
        }
    }
}

Func1& Tabulated1::duplicate() const {
    warn_deprecated("Tabulated1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    if (m_isLinear) {
        return *(new Tabulated1(m_tvec.size(), &m_tvec[0], &m_fvec[0],
                                "linear"));
    } else {
        return *(new Tabulated1(m_tvec.size(), &m_tvec[0], &m_fvec[0],
                                "previous"));
    }
}

Func1& Tabulated1::derivative() const {
    warn_deprecated("Tabulated1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    vector<double> tvec;
    vector<double> dvec;
    size_t siz = m_tvec.size();
    if (m_isLinear) {
        // piece-wise continuous derivative
        if (siz>1) {
            for (size_t i=1; i<siz; i++) {
                double d = (m_fvec[i] - m_fvec[i-1]) /
                  (m_tvec[i] - m_tvec[i-1]);
                tvec.push_back(m_tvec[i-1]);
                dvec.push_back(d);
            }
        }
        tvec.push_back(m_tvec[siz-1]);
        dvec.push_back(0.);
    } else {
        // derivative is zero (ignoring discontinuities)
        tvec.push_back(m_tvec[0]);
        tvec.push_back(m_tvec[siz-1]);
        dvec.push_back(0.);
        dvec.push_back(0.);
    }
    return *(new Tabulated1(tvec.size(), &tvec[0], &dvec[0], "previous"));
}

shared_ptr<Func1> Tabulated1::derivative3() const {
    vector<double> tvec;
    vector<double> dvec;
    size_t siz = m_tvec.size();
    if (m_isLinear) {
        // piece-wise continuous derivative
        if (siz>1) {
            for (size_t i=1; i<siz; i++) {
                double d = (m_fvec[i] - m_fvec[i-1]) /
                  (m_tvec[i] - m_tvec[i-1]);
                tvec.push_back(m_tvec[i-1]);
                dvec.push_back(d);
            }
        }
        tvec.push_back(m_tvec[siz-1]);
        dvec.push_back(0.);
    } else {
        // derivative is zero (ignoring discontinuities)
        tvec.push_back(m_tvec[0]);
        tvec.push_back(m_tvec[siz-1]);
        dvec.push_back(0.);
        dvec.push_back(0.);
    }
    return make_shared<Tabulated1>(tvec.size(), &tvec[0], &dvec[0], "previous");
}

/******************************************************************************/

Gaussian::Gaussian(double A, double t0, double fwhm) : Gaussian1(A, t0, fwhm)
{
    warn_deprecated("Gaussian::Gaussian", "To be removed after Cantera 3.0. "
        "Replaced by 'Gaussian1'.");
}

Gaussian::Gaussian(const Gaussian& b) : Gaussian1(b)
{
    warn_deprecated("Gaussian::Gaussian", "To be removed after Cantera 3.0. "
        "Replaced by 'Gaussian1'.");
}

Func1& Gaussian::duplicate() const {
    warn_deprecated("Gaussian::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Gaussian* np = new Gaussian(*this);
    return *((Func1*)np);
}

Func1& Poly1::duplicate() const {
    warn_deprecated("Poly1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Poly1* np = new Poly1(*this);
    return *((Func1*)np);
}

Func1& Fourier1::duplicate() const {
    warn_deprecated("Fourier1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Fourier1* np = new Fourier1(*this);
    return *((Func1*)np);
}

Func1& Arrhenius1::duplicate() const {
    warn_deprecated("Arrhenius1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Arrhenius1* np = new Arrhenius1(*this);
    return *((Func1*)np);
}

Func1& Periodic1::duplicate() const {
    warn_deprecated("Periodic1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Periodic1* np = new Periodic1(*this);
    return *((Func1*)np);
}

string Func1::write(const string& arg) const
{
    return fmt::format("\\mathrm{{{}}}({})", type(), arg);
}

string Pow1::write(const string& arg) const
{
    if (m_c == 0.5) {
        return "\\sqrt{" + arg + "}";
    }
    if (m_c == -0.5) {
        return "\\frac{1}{\\sqrt{" + arg + "}}";
    }
    if (m_c != 1.0) {
        return fmt::format("\\left({}\\right)^{{{}}}", arg, m_c);
    } else {
        return arg;
    }
}

string Tabulated1::write(const string& arg) const
{
    return fmt::format("\\mathrm{{Tabulated}}({})", arg);
}

string Const1::write(const string& arg) const
{
    return fmt::format("{}", m_c);
}

Func1& Const1::duplicate() const {
    warn_deprecated("Const1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    return *(new Const1(m_c));
}

Func1& Const1::derivative() const {
    warn_deprecated("Const1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* z = new Const1(0.0);
    return *z;
}

string Ratio1::write(const string& arg) const
{
    return "\\frac{" + m_f1->write(arg) + "}{"
           + m_f2->write(arg) + "}";
}

Func1& Ratio1::duplicate() const {
    warn_deprecated("Ratio1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1d = m_f1->duplicate();
    Func1& f2d = m_f2->duplicate();
    return newRatioFunction(f1d, f2d);
}

Func1& Ratio1::derivative() const {
    warn_deprecated("Ratio1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1& a1 = newProdFunction(m_f1->derivative(), m_f2->duplicate());
    Func1& a2 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
    Func1& s = newDiffFunction(a1, a2);
    Func1& p = newProdFunction(m_f2->duplicate(), m_f2->duplicate());
    return newRatioFunction(s, p);
}

shared_ptr<Func1> Ratio1::derivative3() const {
    auto a1 = newProdFunction(m_f1_shared->derivative3(), m_f2_shared);
    auto a2 = newProdFunction(m_f1_shared, m_f2_shared->derivative3());
    auto s = newDiffFunction(a1, a2);
    auto p = newProdFunction(m_f2_shared, m_f2_shared);
    return newRatioFunction(s, p);
}

string Product1::write(const string& arg) const
{
    string s = m_f1->write(arg);
    if (m_f1->order() < order()) {
        s = "\\left(" + s + "\\right)";
    }
    string s2 = m_f2->write(arg);
    if (m_f2->order() < order()) {
        s2 = "\\left(" + s2 + "\\right)";
    }
    return s + " " + s2;
}

Func1& Product1::duplicate() const {
    warn_deprecated("Product1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1d = m_f1->duplicate();
    Func1& f2d = m_f2->duplicate();
    return newProdFunction(f1d, f2d);
}

Func1& Product1::derivative() const {
    warn_deprecated("Product1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1& a1 = newProdFunction(m_f1->duplicate(), m_f2->derivative());
    Func1& a2 = newProdFunction(m_f2->duplicate(), m_f1->derivative());
    return newSumFunction(a1, a2);
}

shared_ptr<Func1> Product1::derivative3() const {
    auto a1 = newProdFunction(m_f1_shared, m_f2_shared->derivative3());
    auto a2 = newProdFunction(m_f2_shared, m_f1_shared->derivative3());
    return newSumFunction(a1, a2);
}

string Sum1::write(const string& arg) const
{
    string s1 = m_f1->write(arg);
    string s2 = m_f2->write(arg);
    if (s2[0] == '-') {
        return s1 + " - " + s2.substr(1,s2.size());
    } else {
        return s1 + " + " + s2;
    }
}

Func1& Sum1::duplicate() const {
    warn_deprecated("Sum1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1d = m_f1->duplicate();
    Func1& f2d = m_f2->duplicate();
    return newSumFunction(f1d, f2d);
}

Func1& Sum1::derivative() const {
    warn_deprecated("Sum1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1& d1 = m_f1->derivative();
    Func1& d2 = m_f2->derivative();
    return newSumFunction(d1, d2);
}

string Diff1::write(const string& arg) const
{
    string s1 = m_f1->write(arg);
    string s2 = m_f2->write(arg);
    if (s2[0] == '-') {
        return s1 + " + " + s2.substr(1,s2.size());
    } else {
        return s1 + " - " + s2;
    }
}

Func1& Diff1::duplicate() const {
    warn_deprecated("Diff1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1d = m_f1->duplicate();
    Func1& f2d = m_f2->duplicate();
    return newDiffFunction(f1d, f2d);
}

Func1& Diff1::derivative() const {
    warn_deprecated("Diff1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    return newDiffFunction(m_f1->derivative(), m_f2->derivative());
}

string Composite1::write(const string& arg) const
{
    string g = m_f2->write(arg);
    return m_f1->write(g);
}

Func1& Composite1::duplicate() const {
    warn_deprecated("Composite1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1d = m_f1->duplicate();
    Func1& f2d = m_f2->duplicate();
    return newCompositeFunction(f1d, f2d);
}

Func1& Composite1::derivative() const {
    warn_deprecated("Composite1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1* d1 = &m_f1->derivative();
    Func1* d3 = &newCompositeFunction(*d1, m_f2->duplicate());
    Func1* d2 = &m_f2->derivative();
    Func1* p = &newProdFunction(*d3, *d2);
    return *p;
}

shared_ptr<Func1> Composite1::derivative3() const {
    auto d1 = m_f1_shared->derivative3();
    auto d2 = m_f2_shared->derivative3();
    auto d3 = newCompositeFunction(d1, m_f2_shared);
    return newProdFunction(d3, d2);
}

string TimesConstant1::write(const string& arg) const
{
    string s = m_f1->write(arg);
    if (m_f1->order() < order()) {
        s = "\\left(" + s + "\\right)";
    }
    if (m_c == 1.0) {
        return s;
    }
    if (m_c == -1.0) {
        return "-"+s;
    }
    char n = s[0];
    if (n >= '0' && n <= '9') {
        s = "\\left(" + s + "\\right)";
    }
    return fmt::format("{}{}", m_c, s);
}

Func1& TimesConstant1::duplicate() const {
    warn_deprecated("TimesConstant1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1 = m_f1->duplicate();
    Func1* dup = new TimesConstant1(f1, m_c);
    return *dup;
}

Func1& TimesConstant1::derivative() const {
    warn_deprecated("TimesConstant1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    Func1& f1d = m_f1->derivative();
    Func1* d = &newTimesConstFunction(f1d, m_c);
    return *d;
}

string PlusConstant1::write(const string& arg) const
{
    if (m_c == 0.0) {
        return m_f1->write(arg);
    }
    return fmt::format("{} + {}", m_f1->write(arg), m_c);
}

Func1& PlusConstant1::duplicate() const {
    warn_deprecated("PlusConstant1::duplicate",
        "To be removed after Cantera 3.0; no longer needed.");
    Func1& f1 = m_f1->duplicate();
    Func1* dup = new PlusConstant1(f1, m_c);
    return *dup;
}

Func1& PlusConstant1::derivative() const {
    warn_deprecated("PlusConstant1::derivative",
        "To be changed after Cantera 3.0; for new behavior, see 'derivative3'.");
    return m_f1->derivative();
}


double Func1::isProportional(TimesConstant1& other)
{
    if (isIdentical(other.func1())) {
        return other.c();
    }
    return 0.0;
}
double Func1::isProportional(Func1& other)
{
    if (isIdentical(other)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

namespace { // restrict scope of helper functions to local translation unit

static bool isConstant(Func1& f)
{
    if (f.type() == "constant") {
        return true;
    } else {
        return false;
    }
}

bool isConstant(const shared_ptr<Func1>& f)
{
    return f->type() == "constant";
}

static bool isZero(Func1& f)
{
    if (f.type() == "constant" && f.c() == 0.0) {
        return true;
    } else {
        return false;
    }
}

bool isZero(const shared_ptr<Func1>& f)
{
    return f->type() == "constant" && f->c() == 0.0;
}

static bool isOne(Func1& f)
{
    if (f.type() == "constant" && f.c() == 1.0) {
        return true;
    } else {
        return false;
    }
}

bool isOne(const shared_ptr<Func1>& f)
{
    return f->type() == "constant" && f->c() == 1.0;
}

static bool isTimesConst(Func1& f)
{
    if (f.type() == "times-constant") {
        return true;
    } else {
        return false;
    }
}

bool isTimesConst(const shared_ptr<Func1>& f)
{
    return f->type() == "times-constant";
}

static bool isExp(Func1& f)
{
    if (f.type() == "exp") {
        return true;
    } else {
        return false;
    }
}

bool isExp(const shared_ptr<Func1>& f)
{
    return f->type() == "exp";
}

static bool isPow(Func1& f)
{
    if (f.type() == "pow") {
        return true;
    } else {
        return false;
    }
}

bool isPow(const shared_ptr<Func1>& f)
{
    return f->type() == "pow";
}

} // end unnamed namespace

Func1& newSumFunction(Func1& f1, Func1& f2)
{
    warn_deprecated("newSumFunction(Func1&, Func1&)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (f1.isIdentical(f2)) {
        return newTimesConstFunction(f1, 2.0);
    }
    if (isZero(f1)) {
        delete &f1;
        return f2;
    }
    if (isZero(f2)) {
        delete &f2;
        return f1;
    }
    double c = f1.isProportional(f2);
    if (c != 0) {
        if (c == -1.0) {
            return *(new Const1(0.0));
        } else {
            return newTimesConstFunction(f1, c + 1.0);
        }
    }
    return *(new Sum1(f1, f2));
}

shared_ptr<Func1> newSumFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (f1->isIdentical(*f2)) {
        return newTimesConstFunction(f1, 2.);
    }
    if (isZero(f1)) {
        return f2;
    }
    if (isZero(f2)) {
        return f1;
    }
    double c = f1->isProportional(*f2);
    if (c != 0.) {
        if (c == -1.) {
            return make_shared<Const1>(0.);
        }
        return newTimesConstFunction(f1, c + 1.);
    }
    return make_shared<Sum1>(f1, f2);
}

Func1& newDiffFunction(Func1& f1, Func1& f2)
{
    warn_deprecated("newDiffFunction(Func1&, Func1&)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (isZero(f2)) {
        delete &f2;
        return f1;
    }
    if (f1.isIdentical(f2)) {
        delete &f1;
        delete &f2;
        return *(new Const1(0.0));
    }
    double c = f1.isProportional(f2);
    if (c != 0.0) {
        if (c == 1.0) {
            return *(new Const1(0.0));
        } else {
            return newTimesConstFunction(f1, 1.0 - c);
        }
    }
    return *(new Diff1(f1, f2));
}

shared_ptr<Func1> newDiffFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (isZero(f2)) {
        return f1;
    }
    if (isZero(f1)) {
        return newTimesConstFunction(f2, -1.);
    }

    if (f1->isIdentical(*f2)) {
        return make_shared<Const1>(0.);
    }
    double c = f1->isProportional(*f2);
    if (c != 0.) {
        if (c == 1.0) {
            return make_shared<Const1>(0.);
        } else {
            return newTimesConstFunction(f1, 1. - c);
        }
    }
    return make_shared<Diff1>(f1, f2);
}

Func1& newProdFunction(Func1& f1, Func1& f2)
{
    warn_deprecated("newProdFunction(Func1&, Func1&)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (isOne(f1)) {
        delete &f1;
        return f2;
    }
    if (isOne(f2)) {
        delete &f2;
        return f1;
    }
    if (isZero(f1) || isZero(f2)) {
        delete &f1;
        delete &f2;
        return *(new Const1(0.0));
    }
    if (isConstant(f1) && isConstant(f2)) {
        double c1c2 = f1.c() * f2.c();
        delete &f1;
        delete &f2;
        return *(new Const1(c1c2));
    }
    if (isConstant(f1)) {
        double c = f1.c();
        delete &f1;
        return newTimesConstFunction(f2, c);
    }
    if (isConstant(f2)) {
        double c = f2.c();
        delete &f2;
        return newTimesConstFunction(f1, c);
    }

    if (isPow(f1) && isPow(f2)) {
        Func1& p = *(new Pow1(f1.c() + f2.c()));
        delete &f1;
        delete &f2;
        return p;
    }

    if (isExp(f1) && isExp(f2)) {
        Func1& p = *(new Exp1(f1.c() + f2.c()));
        delete &f1;
        delete &f2;
        return p;
    }

    bool tc1 = isTimesConst(f1);
    bool tc2 = isTimesConst(f2);

    if (tc1 || tc2) {
        double c1 = 1.0, c2 = 1.0;
        Func1* ff1 = 0, *ff2 = 0;
        if (tc1) {
            c1 = f1.c();
            ff1 = &f1.func1_dup();
            delete &f1;
        } else {
            ff1 = &f1;
        }
        if (tc2) {
            c2 = f2.c();
            ff2 = &f2.func1_dup();
            delete &f2;
        } else {
            ff2 = &f2;
        }
        Func1& p = newProdFunction(*ff1, *ff2);

        if (c1* c2 != 1.0) {
            return newTimesConstFunction(p, c1*c2);
        } else {
            return p;
        }
    } else {
        return *(new Product1(f1, f2));
    }
}

shared_ptr<Func1> newProdFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (isOne(f1)) {
        return f2;
    }
    if (isOne(f2)) {
        return f1;
    }
    if (isZero(f1) || isZero(f2)) {
        return make_shared<Const1>(0.);
    }
    if (isConstant(f1) && isConstant(f2)) {
        return make_shared<Const1>(f1->c() * f2->c());
    }
    if (isConstant(f1)) {
        return newTimesConstFunction(f2, f1->c());
    }
    if (isConstant(f2)) {
        return newTimesConstFunction(f1, f2->c());
    }
    if (isPow(f1) && isPow(f2)) {
        return make_shared<Pow1>(f1->c() + f2->c());
    }
    if (isExp(f1) && isExp(f2)) {
        return make_shared<Exp1>(f1->c() + f2->c());
    }

    bool tc1 = isTimesConst(f1);
    bool tc2 = isTimesConst(f2);

    if (tc1 || tc2) {
        double c1 = 1.0;
        auto ff1 = f1;
        if (tc1) {
            c1 = f1->c();
            ff1 = f1->func1_shared();
        }
        double c2 = 1.0;
        auto ff2 = f2;
        if (tc2) {
            c2 = f2->c();
            ff2 = f2->func1_shared();
        }
        auto p = newProdFunction(ff1, ff2);

        if (c1 * c2 != 1.0) {
            return newTimesConstFunction(p, c1 * c2);
        }
        return p;
    }
    return make_shared<Product1>(f1, f2);
}

Func1& newRatioFunction(Func1& f1, Func1& f2)
{
    warn_deprecated("newRatioFunction(Func1&, Func1&)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (isOne(f2)) {
        return f1;
    }
    if (isZero(f1)) {
        return *(new Const1(0.0));
    }
    if (f1.isIdentical(f2)) {
        delete &f1;
        delete &f2;
        return *(new Const1(1.0));
    }
    if (f1.ID() == PowFuncType && f2.ID() == PowFuncType) {
        return *(new Pow1(f1.c() - f2.c()));
    }
    if (f1.ID() == ExpFuncType && f2.ID() == ExpFuncType) {
        return *(new Exp1(f1.c() - f2.c()));
    }
    return *(new Ratio1(f1, f2));
}

shared_ptr<Func1> newRatioFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (isOne(f2)) {
        return f1;
    }
    if (isZero(f1)) {
        return make_shared<Const1>(0.);
    }
    if (f1->isIdentical(*f2)) {
        return make_shared<Const1>(1.);
    }
    if (isPow(f1) && isPow(f2)) {
        return make_shared<Pow1>(f1->c() - f2->c());
    }
    if (isExp(f1) && isExp(f2)) {
        return make_shared<Exp1>(f1->c() - f2->c());
    }
    return make_shared<Ratio1>(f1, f2);
}

Func1& newCompositeFunction(Func1& f1, Func1& f2)
{
    warn_deprecated("newCompositeFunction(Func1&, Func1&)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (isZero(f1)) {
        delete &f1;
        delete &f2;
        return *(new Const1(0.0));
    }
    if (isConstant(f1)) {
        delete &f2;
        return f1;
    }
    if (isPow(f1) && f1.c() == 1.0) {
        delete &f1;
        return f2;
    }
    if (isPow(f1) && f1.c() == 0.0) {
        delete &f1;
        delete &f2;
        return *(new Const1(1.0));
    }
    if (isPow(f1) && isPow(f2)) {
        double c1c2 = f1.c() * f2.c();
        delete &f1;
        delete &f2;
        return *(new Pow1(c1c2));
    }
    return *(new Composite1(f1, f2));
}

shared_ptr<Func1> newCompositeFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (isZero(f1)) {
        return make_shared<Const1>(0.0);
    }
    if (isConstant(f1)) {
        return f1;
    }
    if (isPow(f1) && f1->c() == 1.0) {
        return f2;
    }
    if (isPow(f1) && f1->c() == 0.0) {
        return make_shared<Const1>(1.);
    }
    if (isPow(f1) && isPow(f2)) {
        return make_shared<Pow1>(f1->c() * f2->c());
    }
    return make_shared<Composite1>(f1, f2);
}

Func1& newTimesConstFunction(Func1& f, double c)
{
    warn_deprecated("newTimesConstFunction(Func1&, double)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (c == 0.0) {
        delete &f;
        return *(new Const1(0.0));
    }
    if (c == 1.0) {
        return f;
    }
    if (f.type() == "times-constant") {
        f.setC(f.c() * c);
        return f;
    }
    return *(new TimesConstant1(f, c));
}

shared_ptr<Func1> newTimesConstFunction(shared_ptr<Func1> f, double c)
{
    if (c == 0.0) {
        return make_shared<Const1>(0.0);
    }
    if (c == 1.0) {
        return f;
    }
    if (f->type() == "times-constant") {
        return make_shared<TimesConstant1>(f->func1_shared(), f->c() * c);
    }
    return make_shared<TimesConstant1>(f, c);
}

Func1& newPlusConstFunction(Func1& f, double c)
{
    warn_deprecated("newPlusConstFunction(Func1&, double)",
        "To be removed after Cantera 3.0; replaced by shared pointer version.");

    if (c == 0.0) {
        return f;
    }
    if (isConstant(f)) {
        double cc = f.c() + c;
        delete &f;
        return *(new Const1(cc));
    }
    if (f.type() == "plus-constant") {
        f.setC(f.c() + c);
        return f;
    }
    return *(new PlusConstant1(f, c));
}

shared_ptr<Func1> newPlusConstFunction(shared_ptr<Func1> f, double c)
{
    if (c == 0.0) {
        return f;
    }
    if (isConstant(f)) {
        return make_shared<Const1>(f->c() + c);
    }
    if (f->type() == "plus-constant") {
        return make_shared<PlusConstant1>(f->func1_shared(), f->c() + c);
    }
    return make_shared<PlusConstant1>(f, c);
}

}
