/**
 * @file ctnum.h
 */
/*
 *      $Id: ctnum.h,v 1.2 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_CTNUM_H
#define CTC_CTNUM_H

#include "clib_defs.h"

extern "C" {

    EEXXTT int DLL_CPREFIX newMatrix(int m, int n);
    EEXXTT int DLL_CPREFIX delMatrix(int i);
    EEXXTT int DLL_CPREFIX matrix_copy(int i);
    EEXXTT int DLL_CPREFIX matrix_assign(int i, int j);
    EEXXTT int DLL_CPREFIX matrix_nRows(int i);
    EEXXTT int DLL_CPREFIX matrix_nColumns(int i);
    EEXXTT int DLL_CPREFIX matrix_resize(int i, int m, int n, double v);
    EEXXTT int DLL_CPREFIX matrix_appendColumn(int i, double* c);
    EEXXTT double DLL_CPREFIX matrix_value(int i, int m, int n);
    EEXXTT double DLL_CPREFIX matrix_setvalue(int i, int m, int n, double v);
    EEXXTT int DLL_CPREFIX matrix_solve(int i1, int i2);
    EEXXTT int DLL_CPREFIX matrix_multiply(int ma, int mb, int mp);
    EEXXTT int DLL_CPREFIX matrix_invert(int ma);

    EEXXTT int DLL_CPREFIX bmatrix_new(int n, int kl, int ku);
    EEXXTT int DLL_CPREFIX bmatrix_del(int i);
    EEXXTT int DLL_CPREFIX bmatrix_copy(int i);
    EEXXTT int DLL_CPREFIX bmatrix_assign(int i, int j);
    EEXXTT int DLL_CPREFIX bmatrix_nRows(int i);
    EEXXTT int DLL_CPREFIX bmatrix_nColumns(int i);
    EEXXTT int DLL_CPREFIX bmatrix_resize(int i, int m, int n, double v);
    EEXXTT double DLL_CPREFIX bmatrix_value(int i, int m, int n);
    EEXXTT double DLL_CPREFIX bmatrix_setvalue(int i, int m, int n, double v);
    EEXXTT int DLL_CPREFIX bmatrix_solve(int ma, int mb);
    EEXXTT int DLL_CPREFIX bmatrix_multiply(int ma, int mb, int mp);


}

#endif
