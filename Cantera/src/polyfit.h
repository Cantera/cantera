
/* C interface for Fortran DPOLFT subroutine */

#ifndef CT_POLYFIT_H
#define CT_POLYFIT_H

#include <vector>
using namespace std;

#include "ct_defs.h"

namespace Cantera {

    doublereal polyfit(int n, doublereal* x, doublereal* y, doublereal* w, 
        int maxdeg, int& ndeg, doublereal eps, doublereal* r);

    template<class D, class R>
    R poly6(D x, R* c) {
        return ((((((c[6]*x + c[5])*x + c[4])*x + c[3])*x + 
                     c[2])*x + c[1])*x + c[0]);
    }

    template<class D, class R>
    R poly8(D x, R* c) {
        return ((((((((c[8]*x + c[7])*x + c[6])*x + c[5])*x + c[4])*x + c[3])*x + 
                     c[2])*x + c[1])*x + c[0]);
    }

    template<class D, class R>
    R poly10(D x, R* c) {
        return ((((((((((c[10]*x + c[9])*x + c[8])*x + c[7])*x 
                         + c[6])*x + c[5])*x + c[4])*x + c[3])*x 
                     + c[2])*x + c[1])*x + c[0]);
    }
    
    template<class D, class R>
    R poly5(D x, R* c) {
        return (((((c[5]*x + c[4])*x + c[3])*x + 
                     c[2])*x + c[1])*x + c[0]);
    }
    
    template<class D, class R>
    R poly4(D x, R* c) {
        return ((((c[4]*x + c[3])*x + 
                     c[2])*x + c[1])*x + c[0]);
    }
    
    template<class D, class R>
    R poly3(D x, R* c) {
        return (((c[3]*x + c[2])*x + c[1])*x + c[0]);
    }
}    
#endif

 
