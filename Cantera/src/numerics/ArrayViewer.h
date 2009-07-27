/**
 *  @file ArrayViewer.h
 *
 *  Header file for class ArrayViewer
 */

/*  $Author: dggoodwin $
 *  $Revision: 1.1 $
 *  $Date: 2007/05/04 14:40:26 $
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_ARRAYVIEWER_H
#define CT_ARRAYVIEWER_H

#include "ct_defs.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "utilities.h"

namespace Cantera { 


    /**
     *  An interface for 2D arrays stored in column-major
     *  (Fortran-compatible) form. This class is designed for
     *  situations when you have a Fortran-compatible 2D array that
     *  you want to be able to easily get/set individual elements or
     *  entire rows or columns.  Instances of ArrayViewer store only a
     *  pointer to the first element in the array; no copy is made of
     *  the array.  @todo This class should probably be renamed, since
     *  it not only views, but also can modify, array elements. This
     *  class is hardly used; candidate for removal.
     */
    class ArrayViewer {

    public:

        typedef doublereal* iterator;
        typedef const doublereal* const_iterator;


        /**
         * Default constructor. Create an empty array viewer.
         */
        ArrayViewer() : m_nrows(0), m_ncols(0) { data = 0; }


        /** 
         *  Constructor. Create an \c m by \c n array viewer for array v.
         */
        ArrayViewer(int m, int n, doublereal* v) 
            : m_nrows(m), m_ncols(n) {
            data = v;
        }

        /// resize the array viewer
        void resize(int n, int m) {
            m_nrows = n;
            m_ncols = m;
        }

        /// set the nth row to array rw
        void setRow(int n, doublereal* rw) {
            for (int j = 0; j < m_ncols; j++) {
                data[m_nrows*j + n] = rw[j];
            }
        }

        /// get the nth row
        void getRow(int n, doublereal* rw) {
            for (int j = 0; j < m_ncols; j++) {
                rw[j] = data[m_nrows*j + n];
            }
        }

        /// set the values in column m to those in array col
        void setColumn(int m, doublereal* col) {
            for (int i = 0; i < m_nrows; i++) {
                data[m_nrows*m + i] = col[i];
            }
        }

        /// get the values in column m
        void getColumn(int m, doublereal* col) {
            for (int i = 0; i < m_nrows; i++) {
                col[i] = data[m_nrows*m + i];
            }
        }
                
        /// Destructor. Does nothing.
        virtual ~ArrayViewer(){}

        doublereal& operator()( int i, int j) {return value(i,j);}
        doublereal operator() ( int i, int j) const {return value(i,j);}

        /// Return a reference to the (i,j) array element.
        doublereal& value( int i, int j) {return data[m_nrows*j + i];}

        /// Return the value of the (i,j) array element.
        doublereal value( int i, int j) const {return data[m_nrows*j + i];}

        /// Number of rows
        size_t nRows() const { return m_nrows; }

        /// Number of columns
        size_t nColumns() const { return m_ncols; }

        iterator begin() { return data; }
        iterator end() { return data + m_nrows*m_ncols; }
        const_iterator begin() const { return data; }
        const_iterator end() const { return data + m_nrows*m_ncols; }

        doublereal* data;

    protected:

        int m_nrows, m_ncols;
    };

    /// output the array
    inline std::ostream& operator<<(std::ostream& s, const ArrayViewer& m) {
        int nr = static_cast<int>(m.nRows());
        int nc = static_cast<int>(m.nColumns());
        int i,j;
        for (i = 0; i < nr; i++) {
            for (j = 0; j < nc; j++) {
                s << m(i,j) << ", ";
            }
            s << std::endl;
        }
        return s;
    }

    /// Multiply the array by a constant.
    inline void operator*=(ArrayViewer& m, doublereal a) {
        scale(m.begin(), m.end(), m.begin(), a);
    }

    /// Increment the entire array by a constant.
    inline void operator+=(ArrayViewer& x, const ArrayViewer& y) {
        sum_each(x.begin(), x.end(), y.begin());
    }

}

#endif



