#ifndef CT_TYPES_H
#define CT_TYPES_H

#include <valarray>
#include <vector>
#include <algorithm>

using namespace std;


namespace Cantera {
  typedef double Real;

  class array_fp : public valarray<Real> {
  public:
    array_fp(size_t n);
  private:
    size_t sz_;
  };
}

#endif
