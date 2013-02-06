/*
 * $Id: surfdemo.cpp 255 2009-11-09 23:36:49Z hkmoffa $
 *
 *  Sample program that solves an implicit problem for surface
 *  site fractions.
 */
#include "cantera/base/ct_defs.h"
#include "cantera/Cantera.h"
#include <iostream>

#include "cantera/IdealGasMix.h"
//#include "cantera/Interface.h"

using namespace Cantera;
using namespace Cantera_CXX;
using namespace std;
using namespace Cantera_CXX;

int main() {

  try {

    cout << CANTERA_VERSION << std::endl;

  }
  catch (CanteraError) {
    showErrors(cout);
  }
  return 0;
}
      
