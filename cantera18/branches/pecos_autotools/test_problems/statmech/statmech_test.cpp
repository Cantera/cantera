/**
 *  @file statmech
 *       test problem for statistical mechanics in cantera
 */

//  Example 
//
// Test case for the statistical mechanics in cantera
//

#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

/*****************************************************************/
/*****************************************************************/

#include "Cantera.h"
#include "transport.h"
#include "IdealGasMix.h"

#include "TransportFactory.h"

using namespace Cantera;
using namespace Cantera_CXX;

int main(int argc, char** argv) 
{
  int k;
  IdealGasMix g("gri30.xml", "gri30_mix");
  int nsp = g.nSpecies();
  double pres = 1.0E5;
  
  // transport
  int log_level = 0;
  Transport * tran = newTransportMgr("Mix", &g, log_level=0);
  
  
}
