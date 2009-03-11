
#ifndef CT_DEFINES_H 
#define CT_DEFINES_H

namespace Cantera {
  typedef double Real;
  typedef int   Integer;
}


#include <valarray>
#include <numeric>
#include <iostream>

namespace Cantera {

  typedef std::valarray<Real>    array_fp;
  typedef std::valarray<Integer> array_int;

}

#endif
