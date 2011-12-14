
// Cantera includes
#include "numerics.h"
#include "Cabinet.h"
inline DenseMatrix* _matrix(int i) {
    return Cabinet<DenseMatrix>::cabinet()->item(i);
}

inline BandMatrix* _bmatrix(int i) {
    return Cabinet<BandMatrix>::cabinet()->item(i);
}

// Build as a DLL under Windows
#ifdef WIN32
#ifdef NO_DLL_BUILD
#define DLL_EXPORT
#else
#define DLL_EXPORT __declspec(dllexport)
#endif
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#else
#define DLL_EXPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999

Cabinet<DenseMatrix>* Cabinet<DenseMatrix>::__storage = 0;
Cabinet<BandMatrix>* Cabinet<BandMatrix>::__storage = 0;

extern "C" {


    ///// Matrix //////

    int DLL_EXPORT newMatrix(int m, int n) {
        DenseMatrix* x = new DenseMatrix(m,n);
        return Cabinet<DenseMatrix>::cabinet()->add(x);
    }

    int DLL_EXPORT delMatrix(int i) {
        Cabinet<DenseMatrix>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT matrix_copy(int i) {
        return Cabinet<DenseMatrix>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT matrix_assign(int i, int j) {
        return Cabinet<DenseMatrix>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT matrix_nRows(int i) {
        return _matrix(i)->nRows();
    }

    int DLL_EXPORT matrix_nColumns(int i) {
        return _matrix(i)->nColumns();
    }

    int DLL_EXPORT matrix_resize(int i, int m, int n, double v) {
        _matrix(i)->resize(m,n,v);
        return 0;
    }

    int DLL_EXPORT matrix_appendColumn(int i, double* c) {
        _matrix(i)->appendColumn(c);
        return 0;
    }

    double DLL_EXPORT matrix_value(int i, int m, int n) {
        return _matrix(i)->value(m,n);
    }

    double DLL_EXPORT matrix_setvalue(int i, int m, int n, double v) {
        _matrix(i)->value(m,n) = v;
        return v;
    }

    int DLL_EXPORT matrix_solve(int i1, int i2) {
        try {
            int info =  solve(*_matrix(i1), *_matrix(i2));
            return info;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }

    int DLL_EXPORT matrix_multiply(int ma, int mb, int mp) {
        try {
            DenseMatrix* a = _matrix(ma);
            DenseMatrix* b = _matrix(mb);
            DenseMatrix* p = _matrix(mp);
            multiply(*a, b->begin(), p->begin());
            return 0;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }

    int DLL_EXPORT matrix_invert(int ma) {
        try {
            invert(*_matrix(ma));
            return 0;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }


    /////////////////  BandMatrix  //////////////////////


    int DLL_EXPORT bmatrix_new(int n, int kl, int ku) {
        BandMatrix* x = new BandMatrix(n, kl, ku);
        return Cabinet<BandMatrix>::cabinet()->add(x);
    }

    int DLL_EXPORT bmatrix_del(int i) {
        Cabinet<BandMatrix>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT bmatrix_copy(int i) {
        return Cabinet<BandMatrix>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT bmatrix_assign(int i, int j) {
        return Cabinet<BandMatrix>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT bmatrix_nRows(int i) {
        return _bmatrix(i)->rows();
    }

    int DLL_EXPORT bmatrix_nColumns(int i) {
        return _bmatrix(i)->columns();
    }

    int DLL_EXPORT bmatrix_resize(int i, int m, int n, double v) {
        _bmatrix(i)->resize(m,n,v);
        return 0;
    }

    double DLL_EXPORT bmatrix_value(int i, int m, int n) {
        return _bmatrix(i)->value(m,n);
    }

    double DLL_EXPORT bmatrix_setvalue(int i, int m, int n, double v) {
        try {
            _bmatrix(i)->value(m,n) = v;
            return v;
        }
        catch (...) { return ERR; }
    }

    int DLL_EXPORT bmatrix_solve(int ma, int mb) {
        try {
            int n = _bmatrix(ma)->nColumns();
            _bmatrix(ma)->solve(n, 
                _matrix(mb)->begin());
            return 0;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }

    int DLL_EXPORT bmatrix_multiply(int ma, int mb, int mp) {
        try {
            BandMatrix* a = _bmatrix(ma);
            DenseMatrix* b = _matrix(mb);
            DenseMatrix* p = _matrix(mp);
            a->mult(b->begin(), p->begin());
            return 0;
        }
        catch (CanteraError) { return -1; }
        catch (...) { return ERR; }
    }

}
