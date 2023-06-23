/**
 * @file ctfunc.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctfunc.h"

#include "cantera/numerics/Func1.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "clib_utils.h"

using namespace Cantera;

typedef Func1 func_t;

typedef Cabinet<Func1> FuncCabinet;
// Assign storage to the Cabinet<Func1> static member
template<> FuncCabinet* FuncCabinet::s_storage = 0;

extern "C" {

    // functions

    int func_new(int type, size_t n, size_t lenp, const double* params)
    {
        try {
            func_t* r=0;
            size_t m = lenp;
            if (type == SinFuncType) {
                vector<double> par = {params[0]};
                r = new Sin1(0, par);
            } else if (type == CosFuncType) {
                vector<double> par = {params[0]};
                r = new Cos1(0, par);
            } else if (type == ExpFuncType) {
                vector<double> par = {params[0]};
                r = new Exp1(0, par);
            } else if (type == PowFuncType) {
                vector<double> par = {params[0]};
                r = new Pow1(0, par);
            } else if (type == ConstFuncType) {
                vector<double> par = {params[0]};
                r = new Const1(0, par);
            } else if (type == FourierFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = new Fourier1(n, par);
            } else if (type == GaussianFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = new Gaussian1(n, par);
            } else if (type == PolyFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = new Poly1(n, par);
            } else if (type == ArrheniusFuncType) {
                vector<double> par(lenp);
                std::copy(params, params + lenp, par.data());
                r = new Arrhenius1(n, par);
            } else if (type == PeriodicFuncType) {
                r = new Periodic1(FuncCabinet::item(n), params[0]);
            } else if (type == SumFuncType) {
                r = new Sum1(FuncCabinet::item(n).duplicate(),
                             FuncCabinet::item(m).duplicate());
            } else if (type == DiffFuncType) {
                r = new Diff1(FuncCabinet::item(n).duplicate(),
                              FuncCabinet::item(m).duplicate());
            } else if (type == ProdFuncType) {
                r = new Product1(FuncCabinet::item(n).duplicate(),
                                 FuncCabinet::item(m).duplicate());
            } else if (type == RatioFuncType) {
                r = new Ratio1(FuncCabinet::item(n).duplicate(),
                               FuncCabinet::item(m).duplicate());
            } else if (type == CompositeFuncType) {
                r = new Composite1(FuncCabinet::item(n).duplicate(),
                                   FuncCabinet::item(m).duplicate());
            } else if (type == TimesConstantFuncType) {
                r = new TimesConstant1(FuncCabinet::item(n).duplicate(), params[0]);
            } else if (type == PlusConstantFuncType) {
                r = new PlusConstant1(FuncCabinet::item(n).duplicate(), params[0]);
            } else {
                throw CanteraError("func_new","unknown function type");
                r = new Func1();
            }
            return FuncCabinet::add(r);
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
            func_t* r = 0;
            r = &FuncCabinet::item(i).derivative();
            return FuncCabinet::add(r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_duplicate(int i)
    {
        try {
            func_t* r = 0;
            r = &FuncCabinet::item(i).duplicate();
            return FuncCabinet::add(r);
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
