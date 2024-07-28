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
            int m = static_cast<int>(lenp);
            int nn = static_cast<int>(n);
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
                vector<double> par(params, params + lenp);
                r = newFunc1("Fourier", par);
            } else if (type == GaussianFuncType) {
                vector<double> par(params, params + lenp);
                r = newFunc1("Gaussian", par);
            } else if (type == PolyFuncType) {
                vector<double> par(params, params + lenp);
                r = newFunc1("polynomial", par);
            } else if (type == ArrheniusFuncType) {
                vector<double> par(params, params + lenp);
                r = newFunc1("Arrhenius", par);
            } else if (type == PeriodicFuncType) {
                r = newFunc1("periodic", FuncCabinet::at(nn), params[0]);
            } else if (type == SumFuncType) {
                r = newFunc1("sum", FuncCabinet::at(nn), FuncCabinet::at(m));
            } else if (type == DiffFuncType) {
                r = newFunc1("diff", FuncCabinet::at(nn), FuncCabinet::at(m));
            } else if (type == ProdFuncType) {
                r = newFunc1("product", FuncCabinet::at(nn), FuncCabinet::at(m));
            } else if (type == RatioFuncType) {
                r = newFunc1("ratio", FuncCabinet::at(nn), FuncCabinet::at(m));
            } else if (type == CompositeFuncType) {
                r = newFunc1("composite", FuncCabinet::at(nn), FuncCabinet::at(m));
            } else if (type == TimesConstantFuncType) {
                r = newFunc1("times-constant", FuncCabinet::at(nn), params[0]);
            } else if (type == PlusConstantFuncType) {
                r = newFunc1("plus-constant", FuncCabinet::at(nn), params[0]);
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

    int func_new_advanced(const char* type, size_t lenp, const double* params)
    {
        try {
            vector<double> par(params, params + lenp);
            return FuncCabinet::add(newFunc1(type, par));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_compound(const char* type, int a, int b)
    {
        try {
            return FuncCabinet::add(
                newFunc1(type, FuncCabinet::at(a), FuncCabinet::at(b)));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_modified(const char* type, int a, double c)
    {
        try {
            return FuncCabinet::add(newFunc1(type, FuncCabinet::at(a), c));
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
            return FuncCabinet::add(FuncCabinet::at(i)->derivative());
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

    int func_write3(int i, const char* arg, size_t len, char* buf)
    {
        try {
            return static_cast<int>(copyString(FuncCabinet::item(i).write(arg), buf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
