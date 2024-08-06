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

shared_ptr<Func1> Func1::derivative() const
{
    throw CanteraError("Func1::derivative",
        "Needs to be overloaded by Func1 specialization.");
}

bool Func1::isIdentical(Func1& other) const
{
    if (type() == "functor" || type() != other.type() || m_c != other.m_c) {
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

int Func1::order() const
{
    return 3;
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

shared_ptr<Func1> Sin1::derivative() const
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

shared_ptr<Func1> Cos1::derivative() const
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

shared_ptr<Func1> Exp1::derivative() const
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

shared_ptr<Func1> Log1::derivative() const
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

shared_ptr<Func1> Pow1::derivative() const
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

Poly13::Poly13(const vector<double>& params)
{
    if (params.size() == 0) {
        throw CanteraError("Poly13::Poly13",
            "Constructor needs an array that is not empty.");
    }
    size_t n = params.size() - 1;
    m_cpoly.resize(n + 1);
    copy(params.data(), params.data() + m_cpoly.size(), m_cpoly.begin());
}

string Poly13::write(const string& arg) const
{
    // write terms in reverse order
    string out = "";
    if (m_cpoly[m_cpoly.size()-1] != 0.) {
        out = fmt::format("{}", m_cpoly[m_cpoly.size()-1]);
    }
    for (size_t n=1; n<m_cpoly.size(); n++) {
        if (m_cpoly[m_cpoly.size()-1-n] == 0.) {
            continue;
        }
        string term;
        if (m_cpoly[m_cpoly.size()-1-n] == 1.) {
            term = fmt::format("{}", arg);
        } else if (m_cpoly[m_cpoly.size()-1-n] == -1.) {
            term = fmt::format("-{}", arg);
        } else {
            term = fmt::format("{}{}", m_cpoly[m_cpoly.size()-1-n], arg);
        }
        if (n > 9) {
            term = fmt::format("{}^{{{}}}", term, n);
        } else if (n > 1) {
            term = fmt::format("{}^{}", term, n);
        }
        if (!out.size()) {
            out = term;
        } else if (out[0] == '-') {
            out = fmt::format("{} - {}", term, out.substr(1));
        } else {
            out = fmt::format("{} + {}", term, out);
        }
    }
    return out;
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

shared_ptr<Func1> Tabulated1::derivative() const {
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

string Ratio1::write(const string& arg) const
{
    return "\\frac{" + m_f1->write(arg) + "}{" + m_f2->write(arg) + "}";
}

shared_ptr<Func1> Ratio1::derivative() const {
    auto a1 = newProdFunction(m_f1->derivative(), m_f2);
    auto a2 = newProdFunction(m_f1, m_f2->derivative());
    auto s = newDiffFunction(a1, a2);
    auto p = newProdFunction(m_f2, m_f2);
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

shared_ptr<Func1> Product1::derivative() const {
    auto a1 = newProdFunction(m_f1, m_f2->derivative());
    auto a2 = newProdFunction(m_f2, m_f1->derivative());
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

string Composite1::write(const string& arg) const
{
    string g = m_f2->write(arg);
    return m_f1->write(g);
}

shared_ptr<Func1> Composite1::derivative() const {
    auto d1 = m_f1->derivative();
    auto d2 = m_f2->derivative();
    auto d3 = newCompositeFunction(d1, m_f2);
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

string PlusConstant1::write(const string& arg) const
{
    if (m_c == 0.0) {
        return m_f1->write(arg);
    }
    return fmt::format("{} + {}", m_f1->write(arg), m_c);
}

double Func1::isProportional(TimesConstant1& other)
{
    warn_deprecated(
        "Func1::isProportional",
        "Deprecated in Cantera 3.1; replaced by internal function.");
    if (isIdentical(*other.func1_shared())) {
        return other.c();
    }
    return 0.0;
}

double Func1::isProportional(Func1& other)
{
    warn_deprecated(
        "Func1::isProportional",
        "Deprecated in Cantera 3.1; replaced by internal function.");
    if (isIdentical(other)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

namespace { // restrict scope of helper functions to local translation unit

bool isConstant(const shared_ptr<Func1>& f)
{
    return f->type() == "constant";
}

bool isZero(const shared_ptr<Func1>& f)
{
    return f->type() == "constant" && f->c() == 0.0;
}

bool isOne(const shared_ptr<Func1>& f)
{
    return f->type() == "constant" && f->c() == 1.0;
}

bool isTimesConst(const shared_ptr<Func1>& f)
{
    return f->type() == "times-constant";
}

bool isExp(const shared_ptr<Func1>& f)
{
    return f->type() == "exp";
}

bool isPow(const shared_ptr<Func1>& f)
{
    return f->type() == "pow";
}

pair<bool, double> isProportional(
    const shared_ptr<Func1>&f1, const shared_ptr<Func1>&f2)
{
    bool tc1 = isTimesConst(f1);
    bool tc2 = isTimesConst(f2);
    if (!tc1 && !tc2) {
        if (f1->isIdentical(*f2)) {
            return std::make_pair(true, 1.);
        }
        return std::make_pair(false, 0.);
    }
    if (!tc1 && tc2) {
        if (f1->isIdentical(*(f2->func1_shared()))) {
            return std::make_pair(true, f2->c());
        }
        return std::make_pair(false, 0.);
    }
    if (tc1 && !tc2) {
        if (f2->isIdentical(*(f1->func1_shared()))) {
            return std::make_pair(true, 1. / f1->c());
        }
        return std::make_pair(false, 0.);
    }
    if (f2->func1_shared()->isIdentical(*(f1->func1_shared()))) {
        return std::make_pair(true, f2->c() / f1->c());
    }
    return std::make_pair(false, 0.);
}

} // end unnamed namespace

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
    if (isConstant(f2)) {
        return newPlusConstFunction(f1, f2->c());
    }
    if (isConstant(f1)) {
        return newPlusConstFunction(f2, f1->c());
    }
    auto prop = isProportional(f1, f2);
    if (prop.first) {
        double c = prop.second;
        if (c == -1.) {
            return make_shared<Const1>(0.);
        }
        return newTimesConstFunction(f1, c + 1.);
    }
    return make_shared<Sum1>(f1, f2);
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
    if (isConstant(f2)) {
        return newPlusConstFunction(f1, -f2->c());
    }
    auto prop = isProportional(f1, f2);
    if (prop.first) {
        double c = prop.second;
        if (c == 1.0) {
            return make_shared<Const1>(0.);
        }
        return newTimesConstFunction(f1, 1. - c);
    }
    return make_shared<Diff1>(f1, f2);
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

shared_ptr<Func1> newRatioFunction(shared_ptr<Func1> f1, shared_ptr<Func1> f2)
{
    if (isOne(f2)) {
        return f1;
    }
    if (isZero(f1)) {
        return make_shared<Const1>(0.);
    }
    if (isZero(f2)) {
        throw CanteraError("newRatioFunction", "Division by zero.");
    }
    if (f1->isIdentical(*f2)) {
        return make_shared<Const1>(1.);
    }
    if (isConstant(f2)) {
        return newTimesConstFunction(f1, 1. / f2->c());
    }
    if (isPow(f1) && isPow(f2)) {
        return make_shared<Pow1>(f1->c() - f2->c());
    }
    if (isExp(f1) && isExp(f2)) {
        return make_shared<Exp1>(f1->c() - f2->c());
    }
    return make_shared<Ratio1>(f1, f2);
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
