/**
 * @file ctvector.h
 *
 * These classes are defined to represent floating-point and integer
 * arrays. They act as simplified replacements for std::vector<double>
 * and std::vector<int>, respectively. These classes are designed for
 * use where it is necessary to access the array of data values
 * directly, which not all implementations of std::vector allow. This
 * is required, for example, to pass the arrays to Fortran or C
 * routines.
 *
 * They are not implemented as templates since only two data types are
 * usually needed, and to allow for different implementations for
 * floating-point and integer arrays (i.e., using BLAS).  Note that
 * the allocated memory is doubled each time new memory is required,
 * as is done also in vector<T>.
 */

#ifndef CT_CTVECTOR_H
#define CT_CTVECTOR_H

#include <iostream>
#include "config.h"

#ifdef DEBUG_MODE
#include "ctexceptions.h"
#endif

namespace ct {

    /**
     *  A class for floating-point arrays
     */
    class ctvector_fp {
    public:

        //typedef unsigned int size_t;
        //typedef double doublereal;
        typedef doublereal* iterator;
        typedef const doublereal* const_iterator;
        typedef doublereal* pointer;
        typedef size_t difference_type;
        typedef doublereal value_type;

        ctvector_fp(size_t n=0);
        ctvector_fp(size_t n, value_type v0);
        ctvector_fp(const ctvector_fp& x);
        ctvector_fp operator=(const ctvector_fp& x);
        virtual ~ctvector_fp();

        value_type operator[](size_t n) const {
#ifdef DEBUG_MODE
            if (n < 0 || n >= _size)
                throw CanteraError("ctvector_fp","index out of range");
#endif
            return _data[n]; 
        }
        value_type& operator[](size_t n) { return _data[n]; }

        void resize(size_t n);
        void resize(size_t n, value_type v0);

        value_type back() const { return _data[_size-1]; }
        const_iterator begin() const { return _data; }
        iterator begin() { return _data; }
        const_iterator end() const { return _data + _size; }
        iterator end() { return _data + _size; }
        void push_back(value_type x);
        size_t size() const { return _size; }
        void clear();
        bool empty() const { return (_size == 0); }

    protected:
        size_t _size, _alloc;
        iterator _data;
    private:
    };


//     /**
//      *  A class for single-precision floating-point arrays
//      */
//     class ctvector_float {
//     public:

//         //typedef unsigned int size_t;
//         //typedef float doublereal;
//         typedef real* iterator;
//         typedef const real* const_iterator;
//         typedef real* pointer;
//         typedef size_t difference_type;
//         typedef real value_type;

//         ctvector_float(size_t n=0);
//         ctvector_float(size_t n, value_type v0);
//         ctvector_float(const ctvector_float& x);
//         ctvector_float operator=(const ctvector_float& x);
//         virtual ~ctvector_float();

//         value_type operator[](size_t n) const { return _data[n]; }
//         value_type& operator[](size_t n) { return _data[n]; }

//         void resize(size_t n);
//         void resize(size_t n, value_type v0);

//         value_type back() const { return _data[_size-1]; }
//         const_iterator begin() const { return _data; }
//         iterator begin() { return _data; }
//         const_iterator end() const { return _data + _size; }
//         iterator end() { return _data + _size; }
//         void push_back(value_type x);
//         size_t size() const { return _size; }
//         void clear();
//         bool empty() const { return (_size == 0); }

//     protected:
//         size_t _size, _alloc;
//         iterator _data;
//     private:
//     };


    // integer arrays
    class ctvector_int {
    public:

        //typedef unsigned int size_t;
        //typedef long int integer;
        typedef integer* iterator;
        typedef const integer* const_iterator;
        typedef integer* pointer;
        typedef size_t difference_type;
        typedef integer value_type;

        ctvector_int(size_t n=0);
        ctvector_int(size_t n, value_type v0);
        ctvector_int(const ctvector_int& x);
        ctvector_int operator=(const ctvector_int& x);
        virtual ~ctvector_int();

        value_type operator[](size_t n) const { 
            return _data[n];
        }
        value_type& operator[](size_t n) { 
            return _data[n];
        }

        void resize(size_t n);
        void resize(size_t n, value_type v0);

        value_type back() const { return _data[_size-1]; }
        const_iterator begin() const { return _data; }
        iterator begin() { return _data; }
        const_iterator end() const { return _data + _size; }
        iterator end() { return _data + _size; }
        void push_back(value_type x);
        size_t size() const { return _size; }
        void clear();
        bool empty() const { return (_size == 0); }
    protected:
        size_t _size, _alloc;
        iterator _data;
    private:
    };

std::ostream& operator<<(std::ostream& s, const ct::ctvector_fp& v);
    //std::ostream& operator<<(std::ostream& s, ct::ctvector_fp& v);
//std::ostream& operator<<(std::ostream& s, const ct::ctvector_float& v);
//std::ostream& operator<<(std::ostream& s, ct::ctvector_float& v);
std::ostream& operator<<(std::ostream& s, const ct::ctvector_int& v);
}

#endif
