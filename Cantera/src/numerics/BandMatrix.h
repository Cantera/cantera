/**
 *  @file BandMatrix.h
 *
 *  Banded matrices.
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
        BandMatrix(size_t n, size_t kl, size_t ku, doublereal v = 0.0);

        /// copy constructor
        BandMatrix(const BandMatrix& y);

        /// Destructor. Does nothing.
        virtual ~BandMatrix(){}

        /// assignment.
        BandMatrix& operator=(const BandMatrix& y);

        void resize(size_t n, size_t kl, size_t ku, doublereal v = 0.0);

        void bfill(doublereal v) {
          std::fill(data.begin(), data.end(), v);
            m_factored = false;
        }

        doublereal& operator()(size_t i, size_t j) {
            return value(i,j);
        }

        doublereal operator() (size_t i, size_t j) const {
            return value(i,j);
        }

        /// Return a reference to element (i,j). Since this method may
        /// alter the element value, it may need to be refactored, so
        /// the flag m_factored is set to false.
        doublereal& value(size_t i, size_t j) {
            m_factored = false;
            if (i + m_ku < j || i > j + m_kl) {
                m_zero = 0.0;
                return m_zero;
            }
            return data[index(i,j)];
        }

        /// Return the value of element (i,j). This method does not
        /// alter the array.
        doublereal value(size_t i, size_t j) const {
            if (i + m_ku < j || i > j + m_kl) return 0.0;
            return data[index(i,j)];
        }

        /// Return the location in the internal 1D array corresponding to
        /// the (i,j) element in the banded array.
        size_t index(size_t i, size_t j) const {
            size_t rw = m_kl + m_ku + i - j;
            return (2*m_kl + m_ku + 1)*j + rw;
        }

        /// Return the value of the (i,j) element for (i,j) within the
        /// bandwidth. For efficiency, this method does not check that
        /// (i,j) are within the bandwidth; it is up to the calling
        /// program to insure that this is true.
        doublereal _value(size_t i, size_t j) const {
            return data[index(i,j)];
        }

        /// Number of rows
        size_t nRows() const { return m_n; }
        /// @deprecated Redundant.
        size_t rows() const { return m_n; }

        /// Number of columns
        size_t nColumns() const { return m_n; }
        /// @deprecated Redundant.
        size_t columns() const { return m_n; }

        /// Number of subdiagonals
        size_t nSubDiagonals() const { return m_kl; }

        /// Number of superdiagonals
        size_t nSuperDiagonals() const { return m_ku; }

        size_t ldim() const { return 2*m_kl + m_ku + 1; }
        vector_int& ipiv() { return m_ipiv; }

        /// Multiply A*b and write result to prod.
        void mult(const double* b, double* prod) const;

        /// Multiply b*A and write result to prod.
        void leftMult(const double* b, double* prod) const;

        int factor();

        //void solve(const vector_fp& b, vector_fp& x);

        int solve(size_t n, const doublereal* b, doublereal* x);
        int solve(size_t n, doublereal* b);

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

    protected:
        vector_fp data;
        vector_fp ludata;
        bool m_factored;


        size_t m_n, m_kl, m_ku;
        doublereal m_zero;
        vector_int     m_ipiv;

    };

    std::ostream& operator<<(std::ostream& s, const BandMatrix& m);

}

#endif
