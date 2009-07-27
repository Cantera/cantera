/**
 *  @file BandMatrix.h
 *
 *  Banded matrices.
 */

/*
 *  $Author: dggoodwin $
 *  $Revision: 1.1 $
 *  $Date: 2007/05/04 14:40:26 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_BANDMATRIX_H
#define CT_BANDMATRIX_H

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
          std::fill(data.begin(), data.end(), v);
            m_factored = false;
        }

        doublereal& operator()( int i, int j) {
            return value(i,j);
        }

        doublereal operator() ( int i, int j) const {
            return value(i,j);
        }

        /// Return a reference to element (i,j). Since this method may
        /// alter the element value, it may need to be refactored, so
        /// the flag m_factored is set to false.
        doublereal& value( int i, int j) {
            m_factored = false;
            if (i < j - m_ku || i > j + m_kl) {
                m_zero = 0.0;
                return m_zero;
            }
            return data[index(i,j)];
        }

        /// Return the value of element (i,j). This method does not
        /// alter the array.
        doublereal value( int i, int j) const {
            if (i < j - m_ku || i > j + m_kl) return 0.0;
            return data[index(i,j)];
        }

        /// Return the location in the internal 1D array corresponding to
        /// the (i,j) element in the banded array.
        int index(int i, int j) const {
            int rw = m_kl + m_ku + i - j;
            return (2*m_kl + m_ku + 1)*j + rw;
        }

        /// Return the value of the (i,j) element for (i,j) within the
        /// bandwidth. For efficiency, this method does not check that
        /// (i,j) are within the bandwidth; it is up to the calling
        /// program to insure that this is true.
        doublereal _value(int i, int j) const {
            return data[index(i,j)];
        }

        /// Number of rows
        int nRows() const { return m_n; }
        /// @deprecated Redundant.
        int rows() const { return m_n; }

        /// Number of columns
        int nColumns() const { return m_n; }
        /// @deprecated Redundant.
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

        int factor();

        //void solve(const vector_fp& b, vector_fp& x);

        int solve(int n, const doublereal* b, doublereal* x);
        int solve(int n, doublereal* b);

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

    std::ostream& operator<<(std::ostream& s, const BandMatrix& m);

}

#endif
