/**
 * @file ctfunc.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CANTERA_USE_INTERNAL
#include "cantera/clib/ctfunc.h"

#include "cantera/numerics/Func1.h"
#include "cantera/base/ctexceptions.h"

#include "Cabinet.h"

using namespace Cantera;
using namespace std;

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
                r = new Sin1(params[0]);
            } else if (type == CosFuncType) {
                r = new Cos1(params[0]);
            } else if (type == ExpFuncType) {
                r = new Exp1(params[0]);
            } else if (type == PowFuncType) {
                if (lenp < 1) {
                    throw CanteraError("func_new",
                                       "exponent for pow must be supplied");
                }
                r = new Pow1(params[0]);
            } else if (type == ConstFuncType) {
                r = new Const1(params[0]);
            } else if (type == FourierFuncType) {
                if (lenp < 2*n + 2) {
                    throw CanteraError("func_new",
                                       "not enough Fourier coefficients");
                }
                r = new Fourier1(n, params[n+1], params[0], params + 1,
                                 params + n + 2);
            } else if (type == GaussianFuncType) {
                if (lenp < 3) {
                    throw CanteraError("func_new",
                                       "not enough Gaussian coefficients");
                }
                r = new Gaussian(params[0], params[1], params[2]);
            } else if (type == PolyFuncType) {
                if (lenp < n + 1) {
                    throw CanteraError("func_new",
                                       "not enough polynomial coefficients");
                }
                r = new Poly1(n, params);
            } else if (type == ArrheniusFuncType) {
                if (lenp < 3*n) {
                    throw CanteraError("func_new",
                                       "not enough Arrhenius coefficients");
                }
                r = new Arrhenius1(n, params);
            } else if (type == PeriodicFuncType) {
                r = new Periodic1(FuncCabinet::item(n), params[0]);
            } else if (type == SumFuncType) {
                r = &newSumFunction(FuncCabinet::item(n).duplicate(),
                                    FuncCabinet::item(m).duplicate());
            } else if (type == DiffFuncType) {
                r = &newDiffFunction(FuncCabinet::item(n).duplicate(),
                                     FuncCabinet::item(m).duplicate());
            } else if (type == ProdFuncType) {
                r = &newProdFunction(FuncCabinet::item(n).duplicate(),
                                     FuncCabinet::item(m).duplicate());
            } else if (type == RatioFuncType) {
                r = &newRatioFunction(FuncCabinet::item(n).duplicate(),
                                      FuncCabinet::item(m).duplicate());
            } else if (type == CompositeFuncType) {
                r = &newCompositeFunction(FuncCabinet::item(n).duplicate(),
                                          FuncCabinet::item(m).duplicate());
            } else if (type == TimesConstantFuncType) {
                r = &newTimesConstFunction(FuncCabinet::item(n).duplicate(), params[0]);
            } else if (type == PlusConstantFuncType) {
                r = &newPlusConstFunction(FuncCabinet::item(n).duplicate(), params[0]);
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
