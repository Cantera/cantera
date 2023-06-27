/**
 * @file ctfunc.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctfunc.h"

#include "cantera/numerics/Func1Factory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "clib_utils.h"

using namespace Cantera;

typedef SharedCabinet<Func1> FuncCabinet;
// Assign storage to the SharedCabinet<Func1> static member
template<> FuncCabinet* FuncCabinet::s_storage = 0;

extern "C" {

    // functions

    int func_new(int type, size_t n, size_t lenp, const double* params)
    {
        try {
            shared_ptr<Func1> r;
            size_t m = lenp;
            if (type == SinFuncType) {
                r = newFunc1("sin", params[0]);
            } else if (type == CosFuncType) {
                r = newFunc1("cos", params[0]);
            } else if (type == ExpFuncType) {
                r = newFunc1("exp", params[0]);
            } else if (type == PowFuncType) {
                r = newFunc1("pow", params[0]);
            } else if (type == ConstFuncType) {
                r = newFunc1("constant", params[0]);
            } else if (type == FourierFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = newFunc1("fourier", par, n);
            } else if (type == GaussianFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = newFunc1("gaussian", par, n);
            } else if (type == PolyFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = newFunc1("polynomial", par, n);
            } else if (type == ArrheniusFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = newFunc1("arrhenius", par, n);
            } else if (type == PeriodicFuncType) {
                r = newMath1("periodic", FuncCabinet::at(n), params[0]);
            } else if (type == SumFuncType) {
                r = newMath1("sum", FuncCabinet::at(n), FuncCabinet::at(m));
            } else if (type == DiffFuncType) {
                r = newMath1("diff", FuncCabinet::at(n), FuncCabinet::at(m));
            } else if (type == ProdFuncType) {
                r = newMath1("product", FuncCabinet::at(n), FuncCabinet::at(m));
            } else if (type == RatioFuncType) {
                r = newMath1("ratio", FuncCabinet::at(n), FuncCabinet::at(m));
            } else if (type == CompositeFuncType) {
                r = newMath1("composite", FuncCabinet::at(n), FuncCabinet::at(m));
            } else if (type == TimesConstantFuncType) {
                r = newMath1("times-constant", FuncCabinet::at(n), params[0]);
            } else if (type == PlusConstantFuncType) {
                r = newMath1("plus-constant", FuncCabinet::at(n), params[0]);
            } else {
                throw CanteraError("func_new", "unknown function type");
            }
            return FuncCabinet::add(r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_basic(const char* type, double c)
    {
        try {
            return FuncCabinet::add(newFunc1(type, c));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_advanced(const char* type, size_t lenp, const double* params, size_t n)
    {
        try {
            vector<double> par(lenp);
            std::copy(params, params + lenp, par.data());
            return FuncCabinet::add(newFunc1(type, par, n));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_math(const char* type, size_t a, size_t b)
    {
        try {
            return FuncCabinet::add(
                newMath1(type, FuncCabinet::at(a), FuncCabinet::at(b)));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_misc(const char* type, size_t a, double c)
    {
        try {
            return FuncCabinet::add(newMath1(type, FuncCabinet::at(a), c));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_del(int i)
    {
        try {
            FuncCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearFunc()
    {
        try {
            FuncCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_type(int i, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(FuncCabinet::item(i).type(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double func_value(int i, double t)
    {
        try {
            return FuncCabinet::item(i).eval(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int func_derivative(int i)
    {
        try {
            return FuncCabinet::add(FuncCabinet::at(i)->derivative3());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_duplicate(int i)
    {
        try {
            return FuncCabinet::add(FuncCabinet::at(i));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_write(int i, size_t lennm, const char* arg, char* nm)
    {
        try {
            copyString(FuncCabinet::item(i).write(arg), nm, lennm);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
