//! @file Func1Factory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/Func1Factory.h"

namespace Cantera
{

Func1Factory* Func1Factory::s_factory = 0;
std::mutex Func1Factory::domain_mutex;

Func1Factory::Func1Factory()
{
    reg("functor", [](size_t n, const vector<double>& params) {
        return new Func1();
    });
    reg("sin", [](size_t n, const vector<double>& params) {
        return new Sin1(n, params);
    });
    reg("cos", [](size_t n, const vector<double>& params) {
        return new Cos1(n, params);
    });
    reg("exp", [](size_t n, const vector<double>& params) {
        return new Exp1(n, params);
    });
    reg("log", [](size_t n, const vector<double>& params) {
        return new Log1(n, params);
    });
    reg("pow", [](size_t n, const vector<double>& params) {
        return new Pow1(n, params);
    });
    reg("constant", [](size_t n, const vector<double>& params) {
        return new Const1(n, params);
    });
    reg("polynomial", [](size_t n, const vector<double>& params) {
        return new Poly1(n, params);
    });
    reg("fourier", [](size_t n, const vector<double>& params) {
        return new Fourier1(n, params);
    });
    reg("gaussian", [](size_t n, const vector<double>& params) {
        return new Gaussian1(n, params);
    });
    reg("arrhenius", [](size_t n, const vector<double>& params) {
        return new Arrhenius1(n, params);
    });
    reg("tabulated-linear", [](size_t n, const vector<double>& params) {
        return new Tabulated1(n, params);
    });
    reg("tabulated-previous", [](size_t n, const vector<double>& params) {
        auto fcn = new Tabulated1(n, params);
        fcn->setMethod("previous");
        return fcn;
    });
}

Func1Factory* Func1Factory::factory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    if (!s_factory) {
        s_factory = new Func1Factory;
    }
    return s_factory;
}

void Func1Factory::deleteFactory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    delete s_factory;
    s_factory = 0;
}

Math1FactoryA* Math1FactoryA::s_factory = 0;
std::mutex Math1FactoryA::domain_mutex;

Math1FactoryA::Math1FactoryA()
{
    reg("sum", [](const shared_ptr<Func1> f1, const shared_ptr<Func1> f2) {
        return new Sum1(f1, f2);
    });
    reg("diff", [](const shared_ptr<Func1> f1, const shared_ptr<Func1> f2) {
        return new Diff1(f1, f2);
    });
    reg("product", [](const shared_ptr<Func1> f1, const shared_ptr<Func1> f2) {
        return new Product1(f1, f2);
    });
    reg("ratio", [](const shared_ptr<Func1> f1, const shared_ptr<Func1> f2) {
        return new Ratio1(f1, f2);
    });
    reg("composite", [](const shared_ptr<Func1> f1, const shared_ptr<Func1> f2) {
        return new Composite1(f1, f2);
    });
}

Math1FactoryA* Math1FactoryA::factory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    if (!s_factory) {
        s_factory = new Math1FactoryA;
    }
    return s_factory;
}

void Math1FactoryA::deleteFactory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    delete s_factory;
    s_factory = 0;
}

Math1FactoryB* Math1FactoryB::s_factory = 0;
std::mutex Math1FactoryB::domain_mutex;

Math1FactoryB::Math1FactoryB()
{
    reg("times-constant", [](const shared_ptr<Func1> f, double c) {
        return new TimesConstant1(f, c);
    });
    reg("plus-constant", [](const shared_ptr<Func1> f, double c) {
        return new PlusConstant1(f, c);
    });
    reg("periodic", [](const shared_ptr<Func1> f, double c) {
        return new Periodic1(f, c);
    });
}

Math1FactoryB* Math1FactoryB::factory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    if (!s_factory) {
        s_factory = new Math1FactoryB;
    }
    return s_factory;
}

void Math1FactoryB::deleteFactory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    delete s_factory;
    s_factory = 0;
}

shared_ptr<Func1> newFunc1(const string& func1Type, double c)
{
    return shared_ptr<Func1>(
        Func1Factory::factory()->create(func1Type, 0, {c}));
}

shared_ptr<Func1> newFunc1(const string& func1Type, size_t n,
                           const vector<double>& params)
{
    return shared_ptr<Func1>(
        Func1Factory::factory()->create(func1Type, n, params));
}

shared_ptr<Func1> newMath1(const string& mathType, const shared_ptr<Func1> f1,
                           const shared_ptr<Func1> f2)
{
    return shared_ptr<Func1>(
        Math1FactoryA::factory()->create(mathType, f1, f2));
}

shared_ptr<Func1> newMath1(const string& mathType, const shared_ptr<Func1> f, double c)
{
    return shared_ptr<Func1>(
        Math1FactoryB::factory()->create(mathType, f, c));
}

}
