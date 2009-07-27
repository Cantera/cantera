/**
 *  @file DenseMatrix.h
 *
 * Dense (not sparse) matrices.
 */

/*
 *  $Author: dggoodwin $
 *  $Date: 2008/02/05 23:36:12 $
 *  $Revision: 1.2 $
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_DENSEMATRIX_H
#define CT_DENSEMATRIX_H

#include "ct_defs.h"
#include "Array.h"

namespace Cantera { 

    /**
     *  A class for full (non-sparse) matrices with Fortran-compatible
     *  data storage. Adds matrix operations to class Array2D.
     */
    class DenseMatrix : public Array2D {

    public:

        DenseMatrix(){}

        /** 
         * Constructor. Create an \c n by \c m matrix, and initialize
         * all elements to \c v.
         */
        DenseMatrix(int n, int m, doublereal v = 0.0) : Array2D(n,m,v) {
            m_ipiv.resize( max(n, m) );
        }

        /// copy constructor
        DenseMatrix(const DenseMatrix& y) : Array2D(y) {
            m_ipiv = y.ipiv();
        }

        /// assignment.
        DenseMatrix& operator=(const DenseMatrix& y);

        void resize(int n, int m, doublereal v = 0.0);

        /// Destructor. Does nothing.
        virtual ~DenseMatrix(){}


        /**
         * Multiply A*b and write result to \c prod.
         */
        virtual void mult(const double* b, double* prod) const;

        /**
         * Left-multiply the matrix by transpose(b), and write the
         * result to prod.
         */
        virtual void leftMult(const double* b, double* prod) const;

        vector_int& ipiv() { return m_ipiv; }
        const vector_int& ipiv() const { return m_ipiv; }

    protected:

        vector_int     m_ipiv;
    };


    /**
     * Solve Ax = b. Array b is overwritten on exit with x.
     */
     int solve(DenseMatrix& A, double* b);

    /**  Solve Ax = b for multiple right-hand-side vectors. */
     int solve(DenseMatrix& A, DenseMatrix& b);

#ifdef INCL_LEAST_SQUARES
    /** @todo fix lwork */
     int leastSquares(DenseMatrix& A, double* b);
#endif
    /**
     * Multiply \c A*b and return the result in \c prod. Uses BLAS
     * routine DGEMV.
     */
     void multiply(const DenseMatrix& A, const double* b, double* prod);

     void increment(const DenseMatrix& A,
        const double* b, double* prod);

    /**
     * invert A. A is overwritten with A^-1.
     */
     int invert(DenseMatrix& A, int nn=-1);

}

#endif



