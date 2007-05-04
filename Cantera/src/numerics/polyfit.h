
/* C interface for Fortran DPOLFT subroutine */

#ifndef CT_POLYFIT_H
#define CT_POLYFIT_H

#include <vector>
//using namespace std;

#include "ct_defs.h"

namespace Cantera {

    doublereal polyfit(int n, doublereal* x, doublereal* y, doublereal* w, 
        int maxdeg, int& ndeg, doublereal eps, doublereal* r);

}    
#endif

 
