#include "cttypes.h"

namespace Cantera {

  array_fp::array_fp(size_t n) : valarray<Real>(n), sz_(2)  {  }
    
}
