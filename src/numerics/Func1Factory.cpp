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
    reg("functor", [](const vector<double>& params, size_t n) {
        return new Func1();
    });
    reg("sin", [](const vector<double>& params, size_t n) {
        return new Sin1(params);
    });
    reg("cos", [](const vector<double>& params, size_t n) {
        return new Cos1(params);
    });
    reg("exp", [](const vector<double>& params, size_t n) {
        return new Exp1(params);
    });
    reg("log", [](const vector<double>& params, size_t n) {
        return new Log1(params);
    });
    reg("pow", [](const vector<double>& params, size_t n) {
        return new Pow1(params);
    });
    reg("constant", [](const vector<double>& params, size_t n) {
        return new Const1(params);
    });
    reg("polynomial", [](const vector<double>& params, size_t n) {
        return new Poly1(params, n);
    });
    reg("fourier", [](const vector<double>& params, size_t n) {
        return new Fourier1(params, n);
    });
    reg("gaussian", [](const vector<double>& params, size_t n) {
        return new Gaussian1(params, n);
    });
    reg("arrhenius", [](const vector<double>& params, size_t n) {
        return new Arrhenius1(params, n);
    });
    reg("tabulated-linear", [](const vector<double>& params, size_t n) {
        return new Tabulated1(params, n);
    });
    reg("tabulated-previous", [](const vector<double>& params, size_t n) {
        auto fcn = new Tabulated1(params, n);
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
        Func1Factory::factory()->create(func1Type, {c}, npos));
}

shared_ptr<Func1> newFunc1(const string& func1Type,
                           const vector<double>& params, size_t n)
{
    return shared_ptr<Func1>(
        Func1Factory::factory()->create(func1Type, params, n));
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
