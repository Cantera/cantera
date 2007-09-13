/**
 *  @file NASA9poly_test
 *       test problem for NASA 9 coefficient formulation
 */




#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;


#define MAX(x,y) (( (x) > (y) ) ? (x) : (y))    

/*****************************************************************/
/*****************************************************************/

#include "Cantera.h"
#include "transport.h"
#include "IdealGasMix.h"

#include "kernel/TransportFactory.h"

using namespace Cantera;

void printDbl(double val) {
  if (fabs(val) < 5.0E-17) {
    cout << " nil";
  } else {
    cout << val;
  }
}

int main(int argc, char** argv) {

  try {

   
    IdealGasMix g("gasNASA9.xml", "gri30_mix");
    int nsp = g.nSpecies();
    double pres = 1.0E5;
    vector_fp Xset(nsp, 0.0);
    Xset[0] =  0.5 ;
    Xset[1] =  0.5;


    g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
 
         
    vector_fp cp_R(nsp, 0.0);
    g.getCp_R(DATA_PTR(cp_R));
    
    printf("Comparisons of H2 calculated via several equivalent classes:\n");
    printf("1500 K and 1 atm:\n");
    printf("         NasaThermo   Nasa9   Nasa9_4reg \n");
    printf("  cp/R: %11.6g %11.6g %11.6g\n", cp_R[0], cp_R[1], cp_R[2]);

    vector_fp H_RT(nsp, 0.0);
    g.getEnthalpy_RT(DATA_PTR(H_RT));

    printf("  H/RT: %11.6g %11.6g %11.6g\n", H_RT[0], H_RT[1], H_RT[2]);

 
    vector_fp S_R(nsp, 0.0);
    g.getEntropy_R(DATA_PTR(S_R));
    printf("  S/R: %11.6g %11.6g %11.6g\n", S_R[0], S_R[1], S_R[2]);


  }
  catch (CanteraError) {
    showErrors(cout);
  }

  return 0;
}
/***********************************************************/
