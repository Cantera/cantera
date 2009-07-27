#ifndef CT_TYPES_H
#define CT_TYPES_H
/**
 *
 */


#include <valarray>
#include <algorithm>

namespace Cantera {

  typedef double  Real;
  typedef long    Integer;
  typedef std::valarray<Real>  Array_FP;
  typedef std::valarray<Integer>  Array_Int;

  namespace units {  
    const Real cm = 0.01;
    const Real microns = 1.0e-6;
    const Real nm = 1.0e-9;
    const Real km = 1.0E3;

    const Real kg = 1.0;
    const Real g  = 1.0E-3;

    const Real s  = 1.0;
    const Real min  = 60.0;
    const Real hr   = 3600.0;
    const Real days = 24*hr;

    const Real kJ   = 1.0E3;
    const Real J    = 1.0;

    const Real kmol = 1.0;
    const Real mol  = 1.0E-3;


  }
}
#endif
