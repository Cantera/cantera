/**
 *  @file BandMatrix.h 
 *
 *  Banded matrices.
 *
 *  $Author$
 *  $Revision$
 *  $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_BANDMATRIX_H
#define CT_BANDMATRIX_H

#include <iostream>
#include <vector>
using namespace std;

#include "ct_defs.h"
#include "ctlapack.h"
#include "utilities.h"
#include "ctexceptions.h"

namespace Cantera { 

    /**
     *  A class for banded matrices. 
     */
    class BandMatrix {

    public:

        BandMatrix();
        BandMatrix(int n, int kl, int ku, doublereal v = 0.0);

        /// copy constructor
        BandMatrix(const BandMatrix& y);

        /// Destructor. Does nothing.
        virtual ~BandMatrix(){}

        /// assignment.
        BandMatrix& operator=(const BandMatrix& y);

        void resize(int n, int kl, int ku, doublereal v = 0.0);

        void bfill(doublereal v) {
            fill(data.begin(), data.end(), v);
            m_factored = false;
        }

        doublereal& operator()( int i, int j) {
            return value(i,j);
        }

        doublereal operator() ( int i, int j) const {
            return value(i,j);
        }

        doublereal& value( int i, int j) {
            m_factored = false;
            if (i < j - m_ku || i > j + m_kl) {
                m_zero = 0.0;
                return m_zero;
            }
            return data[index(i,j)];
        }

        doublereal value( int i, int j) const {
            if (i < j - m_ku || i > j + m_kl) return 0.0;
            return data[index(i,j)];
        }

        int index(int i, int j) const {
            int rw = m_kl + m_ku + i - j;
            return (2*m_kl + m_ku + 1)*j + rw;
        }
            
        doublereal _value(int i, int j) const {
            return data[index(i,j)];
        }

        /// Number of rows
        int nRows() const { return m_n; }
        int rows() const { return m_n; }

        /// Number of columns
        int nColumns() const { return m_n; }
        int columns() const { return m_n; }

        /// Number of subdiagonals
        int nSubDiagonals() const { return m_kl; }

        /// Number of superdiagonals
        int nSuperDiagonals() const { return m_ku; }

        int ldim() const { return 2*m_kl + m_ku + 1; }
        vector_int& ipiv() { return m_ipiv; }
        
        /// Multiply A*b and write result to prod.
        void mult(const double* b, double* prod) const;

        /// Multiply b*A and write result to prod.
        void leftMult(const double* b, double* prod) const;

        void factor();

        //void solve(const vector_fp& b, vector_fp& x);

        void solve(int n, const doublereal* b, doublereal* x);
        void solve(int n, doublereal* b);
            
        vector_fp::iterator begin() { 
            m_factored = false;
            return data.begin(); 
        }
        vector_fp::iterator end() { 
            m_factored = false;
            return data.end(); 
        }
        vector_fp::const_iterator begin() const { return data.begin(); }
        vector_fp::const_iterator end() const { return data.end(); }

        vector_fp data;
        vector_fp ludata;
        bool m_factored;

    protected:

        int m_n, m_kl, m_ku;
        doublereal m_zero;
        vector_int     m_ipiv;

    };

    ostream& operator<<(ostream& s, const BandMatrix& m);

    /**
     * Solve Ax = b. Array b is overwritten on exit with x.
     */
    //    int bsolve(BandMatrix& A, double* b);

}

#endif



