/**
 *  @file BandMatrix.cpp
 *
 *  Banded matrices.
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "BandMatrix.h"
#include "ctlapack.h"
#include "utilities.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "global.h"

using namespace std;

namespace Cantera {


    /// Default constructor.
    BandMatrix::BandMatrix() : m_factored(false), m_n(0), 
        m_kl(0), m_ku(0), m_zero(0.0) { 
        data.clear(); ludata.clear(); 
    }


    /**
     * Constructor. Create an n by n banded matrix.
     * @param n number of rows and columns
     * @param kl number of subdiagonals
     * @param ku number of superdiagonals
     * @param v initial value (default = 0.0)
     */
    BandMatrix::BandMatrix(int n, int kl, int ku, doublereal v)  
        : m_factored(false), m_n(n), m_kl(kl), m_ku(ku) {
        data.resize(n*(2*kl + ku + 1));
        ludata.resize(n*(2*kl + ku + 1));
        fill(data.begin(), data.end(), v);
        fill(ludata.begin(), ludata.end(), 0.0);
        m_ipiv.resize(m_n);
    }

    /// copy constructor
    BandMatrix::BandMatrix(const BandMatrix& y) {
        m_n = y.m_n;
        m_kl = y.m_kl;
        m_ku = y.m_ku;
        data = y.data;
        ludata = y.ludata;
        m_factored = y.m_factored;
        m_ipiv = y.m_ipiv;
    }

    BandMatrix& BandMatrix::operator=(const BandMatrix& y) {
        if (&y == this) return *this;
        m_n = y.m_n;
        m_kl = y.m_kl;
        m_ku = y.m_ku;
        m_ipiv = y.m_ipiv;
        data = y.data;
        ludata = y.ludata;
        m_factored = y.m_factored;
        return *this;
    }

    void BandMatrix::resize(int n, int kl, int ku, doublereal v) {
        m_n = n;
        m_kl = kl;
        m_ku = ku;
        data.resize(n*(2*kl + ku + 1));
        ludata.resize(n*(2*kl + ku + 1));
        m_ipiv.resize(m_n);
        fill(data.begin(), data.end(), v);
        fill(data.begin(), data.end(), 0.0);
        m_factored = false;
    }

                
    /**
     * Multiply A*b and write result to \c prod.
     */
    void BandMatrix::mult(const double* b, double* prod) const {
        int nr = rows();
        int m, j;
        double sum = 0.0;
        for (m = 0; m < nr; m++) {
            sum = 0.0;
            for (j = m - m_kl; j <= m + m_ku; j++) {
                if (j >= 0 && j < m_n)
                    sum += _value(m,j)*b[j];
            }
            prod[m] = sum;
        }
    }


    /**
     * Multiply b*A and write result to \c prod.
     */
    void BandMatrix::leftMult(const double* b, double* prod) const {
        int nc = columns();
        int n, i;
        double sum = 0.0;
        for (n = 0; n < nc; n++) {
            sum = 0.0;
            for (i = n - m_ku; i <= n + m_kl; i++) {
                if (i >= 0 && i < m_n)
                    sum += _value(i,n)*b[i];
            }
            prod[n] = sum;
        }
    }


    /**
     * Perform an LU decomposition. LAPACK routine DGBTRF is used.
     * The factorization is saved in ludata.
     */
    int BandMatrix::factor() {
        int info=0;
        copy(data.begin(), data.end(), ludata.begin());
        ct_dgbtrf(rows(), columns(), nSubDiagonals(), nSuperDiagonals(), 
            DATA_PTR(ludata), ldim(), DATA_PTR(ipiv()), info);

        // if info = 0, LU decomp succeeded. 
        if (info == 0) {
            m_factored = true;
        }
        else {
	  m_factored = false;
          ofstream fout("bandmatrix.csv");
          fout << *this << endl;
          fout.close();
        }
	return info;
    }



    int BandMatrix::solve(int n, const doublereal* b, doublereal* x) {
        copy(b, b+n, x);
        return solve(n, x);
    }

    int BandMatrix::solve(int n, doublereal* b) {
        int info = 0;
        if (!m_factored) info = factor();
	if (info == 0)
	  ct_dgbtrs(ctlapack::NoTranspose, columns(), nSubDiagonals(), 
              nSuperDiagonals(), 1, DATA_PTR(ludata), ldim(), 
              DATA_PTR(ipiv()), b, columns(), info);

        // error handling
        if (info != 0) {
            ofstream fout("bandmatrix.csv");
            fout << *this << endl;
            fout.close();
        }
	return info;
    }

    ostream& operator<<(ostream& s, const BandMatrix& m) {
        int nr = m.rows();
        int nc = m.columns();
        int i,j;
        for (i = 0; i < nr; i++) {
            for (j = 0; j < nc; j++) {
                s << m(i,j) << ", ";
            }
            s << endl;
        }
        return s;
    }
}

