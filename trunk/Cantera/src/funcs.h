#ifndef CT_FUNCS_H
#define CT_FUNCS_H

#include "ct_defs.h"
 
namespace Cantera {
    doublereal linearInterp(doublereal x, const vector_fp& xpts, 
        const vector_fp& fpts);
}

#endif
