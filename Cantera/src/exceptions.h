/**
 *  @file exceptions.h
 *
 *      @deprecated
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_DEBUG_EXC_H
#define CT_DEBUG_EXC_H

#include <vector>

#ifdef WIN32 
#define _TYPENAME_ 
#else
#define _TYPENAME_ typename
#endif

__BEGIN_DEBUG_NAMESPACE

template<class T>
void show_values(const T& x, size_t max_show = 20) {
        
        cerr << "values: <";
        if (x.size() < max_show) {
                copy(x.begin(), x.begin() + x.size(), 
                        ostream_iterator<_TYPENAME_ T::value_type>(cerr, " "));
        }
        else {
                size_t m2 = max_show/2;
                copy(x.begin(), x.begin() + m2, 
                        ostream_iterator<_TYPENAME_ T::value_type>(cerr, " "));
                cerr << "...(skipping " << x.size() - max_show << ")... ";   
                copy(x.end() - max_show + m2, x.end(), 
                        ostream_iterator<_TYPENAME_ T::value_type>(cerr, " "));
        }   
        cerr << ">" << endl << endl;
}
            
    
template<class T>
class RangeError {
        
public:

        RangeError(const T& vec, typename T::size_type n) {

                cerr << "\n###RANGE ERROR###\n"  
                     << "attempt to access element " 
                     << n << " outside valid range [0,"
                     << vec.size()-1 << "]" << endl;
                cerr << "object at " << &vec << endl;
                show_values(vec);
        }
};
        
    
template<class T>
class SizeError {

public:

        SizeError(const T& vec, typename T::size_type n) {

                cerr << "\n###SIZE ERROR###\n"  
                     << "size = " << vec.size() << ", but should be " << n <<endl;
                cerr << "object at " << &vec << endl;
                show_values(vec);
        }
};
    
__END_DEBUG_NAMESPACE

#endif

