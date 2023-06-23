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
    reg("sin", [](size_t n, const vector<double>& params) {
        return new Sin1(n, params);
    });
    reg("cos", [](size_t n, const vector<double>& params) {
        return new Cos1(n, params);
    });
    reg("exp", [](size_t n, const vector<double>& params) {
        return new Exp1(n, params);
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
    reg("times-constant", [](const shared_ptr<Func1> f1, double c) {
        return new TimesConstant1(f1, c);
    });
    reg("plus-constant", [](const shared_ptr<Func1> f1, double c) {
        return new PlusConstant1(f1, c);
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

}
