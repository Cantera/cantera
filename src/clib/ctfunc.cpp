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

    int func_check(const char* type, int len, char* buf)
    {
        try {
            return static_cast<int>(copyString(checkFunc1(type), buf, len));
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

    int func_new_advanced(const char* type, int lenp, const double* params)
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

    int func_new_sum(int a, int b)
    {
        try {
            return FuncCabinet::add(
                newSumFunction(FuncCabinet::at(a), FuncCabinet::at(b)));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_diff(int a, int b)
    {
        try {
            return FuncCabinet::add(
                newDiffFunction(FuncCabinet::at(a), FuncCabinet::at(b)));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_prod(int a, int b)
    {
        try {
            return FuncCabinet::add(
                newProdFunction(FuncCabinet::at(a), FuncCabinet::at(b)));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int func_new_ratio(int a, int b)
    {
        try {
            return FuncCabinet::add(
                newRatioFunction(FuncCabinet::at(a), FuncCabinet::at(b)));
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

    int func_type(int i, int lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(FuncCabinet::at(i)->type(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double func_value(int i, double t)
    {
        try {
            return FuncCabinet::at(i)->eval(t);
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

    int func_write(int i, const char* arg, int len, char* buf)
    {
        try {
            return static_cast<int>(
                copyString(FuncCabinet::at(i)->write(arg), buf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
