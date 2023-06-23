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

}
