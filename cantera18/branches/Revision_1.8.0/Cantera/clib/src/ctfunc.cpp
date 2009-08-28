/**
 * @file ctfunc.cpp
 */
/*
 *      $Id: ctfunc.cpp,v 1.10 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#pragma warning(disable:4996)
#endif

#define CANTERA_USE_INTERNAL
#include "ctfunc.h"

#include "Func1.h"
#include "ctexceptions.h"

using namespace Cantera;

#include "Cabinet.h"


typedef Func1 func_t;

// Assign storage to the Cabinet<Func1> static member
template<> Cabinet<func_t>*       Cabinet<func_t>::__storage = 0;

inline func_t* _func(int i) {
    return Cabinet<func_t>::cabinet()->item(i);
}
 
extern "C" {

    // functions

    int DLL_EXPORT func_new(int type, int n, int lenp, double* params) {
        func_t* r=0;
        int m = lenp;
        try {
            if (type == SinFuncType) {
                r = new Sin1(params[0]);
            }
            else if (type == CosFuncType) {
                r = new Cos1(params[0]);
            }
            else if (type == ExpFuncType) {
                r = new Exp1(params[0]);
            }
            else if (type == PowFuncType) {
                if (lenp < 1) 
                    throw CanteraError("func_new", 
                        "exponent for pow must be supplied");                
                r = new Pow1(params[0]);
            }
            else if (type == ConstFuncType) {
                r = new Const1(params[0]);
            }
            else if (type == FourierFuncType) {
                if (lenp < 2*n + 2) 
                    throw CanteraError("func_new", 
                        "not enough Fourier coefficients");
                r = new Fourier1(n, params[n+1], params[0], params + 1, 
                    params + n + 2);
            }
            else if (type == GaussianFuncType) {
                if (lenp < 3) 
                    throw CanteraError("func_new", 
                        "not enough Gaussian coefficients");
                r = new Gaussian(params[0], params[1], params[2]);
            }
            else if (type == PolyFuncType) {
                if (lenp < n + 1) 
                    throw CanteraError("func_new", 
                        "not enough polynomial coefficients");
                r = new Poly1(n, params);
            }
            else if (type == ArrheniusFuncType) {
                if (lenp < 3*n) 
                    throw CanteraError("func_new", 
                        "not enough Arrhenius coefficients");
                r = new Arrhenius1(n, params);
            }
            else if (type == PeriodicFuncType) {
                r = new Periodic1(*_func(n), params[0]);
            }
            else if (type == SumFuncType) {
                r = &newSumFunction(_func(n)->duplicate(),
                    _func(m)->duplicate());
            }
            else if (type == DiffFuncType) {
                r = &newDiffFunction(_func(n)->duplicate(), 
                    _func(m)->duplicate());
            }
            else if (type == ProdFuncType) {
                r = &newProdFunction(_func(n)->duplicate(), 
                    _func(m)->duplicate());
            }
            else if (type == RatioFuncType) {
                r = &newRatioFunction(_func(n)->duplicate(), 
                    _func(m)->duplicate());
            }
            else if (type == CompositeFuncType) {
                r = &newCompositeFunction(_func(n)->duplicate(), 
                    _func(m)->duplicate());
            }
            else if (type == TimesConstantFuncType) {
                r = &newTimesConstFunction(_func(n)->duplicate(), params[0]);
            }
            else if (type == PlusConstantFuncType) {
                r = &newPlusConstFunction(_func(n)->duplicate(), params[0]);
            }
            else {
                throw CanteraError("func_new","unknown function type");
                r = new Func1();
            }
            return Cabinet<func_t>::cabinet()->add(r);
        }
        catch (CanteraError) {return -1;}
    }


    int DLL_EXPORT func_del(int i) {
        Cabinet<func_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT func_copy(int i) {
        return Cabinet<func_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT func_assign(int i, int j) {
        return Cabinet<func_t>::cabinet()->assign(i,j);
    }

    double DLL_EXPORT func_value(int i, double t) {
        return _func(i)->eval(t);
    }

    int DLL_EXPORT func_derivative(int i) {
        func_t* r = 0;
        r = &_func(i)->derivative();
        return Cabinet<func_t>::cabinet()->add(r);
    }

    int DLL_EXPORT func_duplicate(int i) {
        func_t* r = 0;
        r = &_func(i)->duplicate();
        return Cabinet<func_t>::cabinet()->add(r);
    }

    int DLL_EXPORT func_write(int i, int lennm, const char* arg, char* nm) {
        try {
            string a = string(arg);
            string w = _func(i)->write(a);
            int ws = w.size();
            int lout = (lennm > ws ? ws : lennm);
			std::copy(w.c_str(), w.c_str() + lout, nm);
            nm[lout] = '\0';
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

}
