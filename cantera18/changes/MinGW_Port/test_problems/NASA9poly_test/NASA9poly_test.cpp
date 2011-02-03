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
using namespace Cantera_CXX;

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
    
    vector_fp Gvalues(nsp, 0.0);

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

   Transport * tran = newTransportMgr("Mix", &g);

    // MultiTransport * tranMix = dynamic_cast<MultiTransport *>(tran);
    printf("Viscoscity and thermal Cond vs. T\n");
    for ( int k = 0; k < 40; k++) {
      double T1 = 400. + 200. * k;
      g.setState_TPX(T1, pres, DATA_PTR(Xset));
      g.getPureGibbs(DATA_PTR(Gvalues));
      //printf("     --  %13g %13.5g %13.5g %13.5g %13.5g \n",
      //	     Gvalues[0],  Gvalues[1],  Gvalues[2], 
      //	     Gvalues[3],  Gvalues[4]); 
      double visc = tran->viscosity();
      double cond = tran->thermalConductivity();
      printf("    %13g %13.5g %13.5g\n", T1, visc, cond);
    }

  }
  catch (CanteraError) {
    showErrors(cout);
  }

  return 0;
}
/***********************************************************/
