/**
 *  @file sort.h
 *
 *  $Id: sort.h,v 1.3 2006/10/20 21:42:37 hkmoffa Exp $
 */

#ifndef CT_SORT_H
#define CT_SORT_H

#include "ct_defs.h"

namespace Cantera {

    /// Given two arrays x and y, sort the (x,y) pairs by the x
    /// values. This version is for floating-point x, and integer y.
    void heapsort(vector_fp& x, vector_int& y);

    /// Given two arrays x and y, sort the (x,y) pairs by the x
    /// values. This version is for floating-point x, and
    /// floating-point y.
    void heapsort(vector_fp& x, vector_fp& y);
}

#endif
