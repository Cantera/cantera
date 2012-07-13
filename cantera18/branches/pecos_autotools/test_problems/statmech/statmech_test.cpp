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
#include "equil.h"

#include "TransportFactory.h"

using namespace Cantera;
using namespace Cantera_CXX;

int main(int argc, char** argv) 
{

  try
    {
      int k;
      IdealGasMix g("test_stat.xml");
      int nsp = g.nSpecies();
      double pres = 1.0E5;

      vector_fp Xset(nsp, 0.0);
      Xset[0] =  0.5 ;
      Xset[1] =  0.5;
  
      g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
      equilibrate(g, "TP", -1);

      vector_fp cp_R(nsp, 0.0);
      g.getCp_R(DATA_PTR(cp_R));

      //for(int i=0;i<nsp;i++)
      //{
      //  std::cout.precision(10);
      //  std::cout << cp_R[i] << std::endl;
      //	}  

      // error check-- exactly 2.5 for atoms
      if(cp_R[0] != 2.5)
	{
	  std::cout << "Error for monotomic Species!\n";
	  return 1;
	}

      // error check: analytical result is more complicated for 
      // molecules. One species should suffice, lets try NO2, with
      // three vibrational modes:
      /// theta[0]: 1.07900e3
      /// theta[1]: 1.9000003
      /// theta[2]: 2.32700e3
      // at T = 2000
      //
      // This is precisely: 6.655804161 (e.g. 5/2 + 1 + 3.1558..)
      //
      double sol = 6.655804161; 
      double tol = 1e-10;

      if(cp_R[3] - sol >= tol )
	{
	  double diff = cp_R[3]-6.655804161;
	  std::cout << "Error for Species NO2!\n";
	  std::cout << "Diff was: " << diff << "\n";
	  return 1;
	}

    }
    catch (CanteraError) 
      {
	showErrors(cout);
	return 1;
      }

  // Mark it zero!
  return 0;

}
