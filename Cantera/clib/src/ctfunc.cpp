
#include "Func1.h"
#include "ctexceptions.h"
using namespace Cantera;

#include "Cabinet.h"

// Build as a DLL under Windows
#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#else
#define DLL_EXPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999

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
                r = new Sum1(*_func(n), *_func(m));
            }
            else if (type == DiffFuncType) {
                r = new Diff1(*_func(n), *_func(m));
            }
            else if (type == ProdFuncType) {
                r = new Product1(*_func(n), *_func(m));
            }
            else if (type == RatioFuncType) {
                r = new Ratio1(*_func(n), *_func(m));
            }
            else if (type == CompositeFuncType) {
                r = new Composite1(*_func(n), *_func(m));
            }
            else if (type == TimesConstantFuncType) {
                r = new TimesConstant1(*_func(n), params[0]);
            }
            else if (type == PlusConstantFuncType) {
                r = new PlusConstant1(*_func(n), params[0]);
            }
            else 
                r = new Func1();
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


}
