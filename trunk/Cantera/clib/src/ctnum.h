#ifndef CTC_CTNUM_H
#define CTC_CTNUM_H

#include "clib_defs.h"

extern "C" {

    int DLL_IMPORT newMatrix(int m, int n);
    int DLL_IMPORT delMatrix(int i);
    int DLL_IMPORT matrix_copy(int i);
    int DLL_IMPORT matrix_assign(int i, int j);
    int DLL_IMPORT matrix_nRows(int i);
    int DLL_IMPORT matrix_nColumns(int i);
    int DLL_IMPORT matrix_resize(int i, int m, int n, double v);
    int DLL_IMPORT matrix_appendColumn(int i, double* c);
    double DLL_IMPORT matrix_value(int i, int m, int n);
    double DLL_IMPORT matrix_setvalue(int i, int m, int n, double v);
    int DLL_IMPORT matrix_solve(int i1, int i2);
    int DLL_IMPORT matrix_multiply(int ma, int mb, int mp);
    int DLL_IMPORT matrix_invert(int ma);

    int DLL_IMPORT bmatrix_new(int n, int kl, int ku);
    int DLL_IMPORT bmatrix_del(int i);
    int DLL_IMPORT bmatrix_copy(int i);
    int DLL_IMPORT bmatrix_assign(int i, int j);
    int DLL_IMPORT bmatrix_nRows(int i);
    int DLL_IMPORT bmatrix_nColumns(int i);
    int DLL_IMPORT bmatrix_resize(int i, int m, int n, double v);
    double DLL_IMPORT bmatrix_value(int i, int m, int n);
    double DLL_IMPORT bmatrix_setvalue(int i, int m, int n, double v);
    int DLL_IMPORT bmatrix_solve(int ma, int mb);
    int DLL_IMPORT bmatrix_multiply(int ma, int mb, int mp);


}

#endif
