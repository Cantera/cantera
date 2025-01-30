/**
 *  @file CustomKinetics.cpp 
 *
 * @ingroup chemkinetics
 */

//
// Author: Q. Cazeres, A. Felden, P. Pepiot 
//
//
#include "cantera/kinetics/CustomKinetics.h"
// #include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"

#include <iostream>
#include <dlfcn.h>
 
using namespace std;

namespace Cantera 
{
  CustomKinetics::CustomKinetics() : GasKinetics() 
  {
    printf("WARNING: Using customized kinetics from f90 file.\n");
    handle = dlopen("customkinetics.so", RTLD_LAZY);

    // load symbol 
    ck = (ck_t) dlsym(handle, "customkinetics_");

  }

  void CustomKinetics::get_wdot_reduced(double* wdot)
  {

    double P = thermo().pressure();
    double T = thermo().temperature();
    // New yarc2 format
    //double rho = thermo().density();
    const double* m_y = thermo().massFractions();

    ck(&P,&T,&m_y[0],&wdot[0]);

    // New yarc2 format
    //ck(&P,&T,&rho,&m_y[0],&wdot[0]);

    // Old yarc format
    //mol/kmol conversion cantera is in kmol
    for (size_t i=0;i<thermo().nSpecies();i++) {
    	wdot[i]=wdot[i]/1000.0;
    }
  }

  void CustomKinetics::close_dl()
  {
   dlclose(handle);
  }

}
