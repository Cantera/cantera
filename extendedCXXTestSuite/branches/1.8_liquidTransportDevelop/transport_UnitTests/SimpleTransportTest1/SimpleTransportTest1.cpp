/* ======================================================================= */
/* $RCSfile: SimpleTransportTest1.cpp,v $ */
/* $Author$ */
/* $Date$ */
/* $Revision$ */
/* ======================================================================= */

#include <stdio.h>


#include "cantera/Cantera.h"
#include "cantera/thermo.h"

#include "cantera/equilibrium.h"

#include "cantera/thermo.h"

#include "cantera/kernel/HMWSoln.h"
#include "cantera/kernel/vcs_internal.h"
#include "cantera/kernel/logger.h"

#include "cantera/transport.h"
#include "cantera/kernel/SimpleTransport.h"



using namespace std;
using namespace Cantera;

int CHECK_DEBUG_MODE = 0;

void printUsage() {
    cout << "usage: Margules_test_4 " <<  endl;
    cout <<"                " << endl;
}


int main(int argc, char **argv)
{
  int k;
  int retn = 0;

  try {
    char iFile[80];
    strcpy(iFile, "HMW_NaCl_pdss.xml");
    if (argc > 1) {
      strcpy(iFile, argv[1]);
    }


    double x[20];
    for (k = 0; k < 20; k++) {
      x[k] = 0.0;
    }
    //double um[20];
    //double pres = OneAtm;
    int log_level = 3;

    HMWSoln *HMW = new HMWSoln(iFile, "NaCl_electrolyte");

    Transport * tran = newDefaultTransportMgr(HMW, log_level);

    SimpleTransport * tranSimple = dynamic_cast<SimpleTransport *>(tran);
    int nsp = HMW->nSpecies();

    HMW->setState_TP(30+273.13, OneAtm);

    double visc = tranSimple->viscosity();
    printf("visc = %g\n", visc);



    tranSimple->getSpeciesViscosities(x);
    for (k = 0; k < nsp; k++) {
      string spName = HMW->speciesName(k);
      printf("sp visc (%s) = %g\n", spName.c_str(), x[k]);
    }


    double cond = tranSimple->thermalConductivity();
    printf("cond = %g\n", cond);



    tranSimple->getMixDiffCoeffs(x);
    for (k = 0; k < nsp; k++) {
      string spName = HMW->speciesName(k);
      printf("sp diff (%s) = %g\n", spName.c_str(), x[k]);
    }

    tranSimple->getMobilities(x);
    for (k = 0; k < nsp; k++) {
      string spName = HMW->speciesName(k);
      printf("Mobility (%s) = %g\n", spName.c_str(), x[k]);
    }

    double gradX[20];
    for (k = 0; k < 20; k++) {
      gradX[k] = 0.0;
    }
    gradX[1] = 1.0;
    double gradT = 0.0;

    tranSimple->getSpeciesFluxes(1, &gradT, 5, &gradX[0], 5, x);
    for (k = 0; k < nsp; k++) {
      string spName = HMW->speciesName(k);
      printf("SpeciesFlux (%s) = %g\n", spName.c_str(), x[k]);
    }

   
    gradX[1] = 0.0;
    double gradV = 1.0;

    tranSimple->set_Grad_T(&gradT);
    tranSimple->set_Grad_V(&gradV);
    tranSimple->set_Grad_X(gradX);

    tranSimple->getSpeciesFluxesExt(5, x);
    for (k = 0; k < nsp; k++) {
      string spName = HMW->speciesName(k);
      printf("SpeciesFlux (%s) = %g\n", spName.c_str(), x[k]);
    }

    delete HMW;
    HMW=0;

    delete tranSimple;
    tranSimple = 0;

    Cantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
