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
  IdealGasMix g("test.xml");
  int nsp = g.nSpecies();
  double pres = 1.0E5;

  //vector_fp Xset(nsp, 0.0);
  //Xset[0] =  0.5 ;
  //Xset[1] =  0.5;
  
  //g.setState_TPX(1500.0, pres, DATA_PTR(Xset));


  vector_fp cp_R(nsp, 0.0);
  g.getCp_R(DATA_PTR(cp_R));
  
  //vector_fp S_R(nsp, 0.0);
  //g.getEntropy_R(DATA_PTR(S_R));
  //printf("  S/R: %11.6g %11.6g %11.6g\n", S_R[0], S_R[1], S_R[2]);


}
