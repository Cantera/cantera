/**
 *  @file sort.h
 */

#ifndef CT_SORT_H
#define CT_SORT_H

#include "ct_defs.h"

namespace Cantera {

  // sort (x,y) pairs by x

    void heapsort(vector_fp& x, vector_int& y);
    void heapsort(vector_fp& x, vector_fp& y);
}

#endif
