/* ======================================================================= */
/* $RCSfile: SimpleTransportTest1.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2009/10/13 01:49:33 $ */
/* $Revision: 1.1 $ */
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
#include "cantera/kernel/DustyGasTransport.h"
#include "cantera/kernel/TransportFactory.h"


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
    strcpy(iFile, "h2o2.xml");
    if (argc > 1) {
      strcpy(iFile, argv[1]);
    }


    double x[20];
    for (k = 0; k < 20; k++) {
      x[k] = 0.0;
    }
    double um[20];
    double pres = OneAtm;
    int log_level = 3;


#ifdef CANTERA_HAS_TEMPLATES
    thermo_t *g = newPhase<doublereal>(iFile, "ohmech");
#else
    thermo_t *g = newPhase(iFile, "ohmech");
#endif

    double Temp = 500;
    Transport* tran = newTransportMgr("DustyGas", g, log_level);

    DustyGasTransport * tranDusty = dynamic_cast<DustyGasTransport *>(tran);
    int nsp = g->nSpecies();
    std::vector<double> multiD(nsp*nsp);

    g->setState_TP(Temp, OneAtm);

    g->setMoleFractionsByName("OH:1, H:2, O2:3, O:1.0E-8, H2:1.0E-8, H2O:1.0E-8, H2O2:1.0E-8, HO2:1.0E-8, AR:1.0E-8");

    /*
     *  Again, the density is a primative variable. So we need to reset the pressure to OneAtm again after
     *  changing the mole fractions.
     */
    g->setState_TP(Temp, OneAtm);

    //! Set the Parameters in the model
    /*!
     *    @param type     Type of the parameter to set
     *                     0 - porosity
     *                     1 - tortuosity
     *                     2 - mean pore radius
     *                     3 - mean particle radius
     *                     4 - permeability
     *    @param k        Unused int
     *    @param p         pointer to double for the input list of parameters
     *
     */
    double por = 0.2;
    tran->setParameters(0, 0, &por);

    double tort = 4.0;
    tran->setParameters(1, 0, &tort);

   
    //! Set the mean pore radius (m)
    /*!
     *     @param   rbar    Value of the pore radius ( m)
     */
    double meanPoreRadius = 1.5E-7;
    tranDusty->setMeanPoreRadius(meanPoreRadius);

    //! Set the mean particle diameter
    /*!
     *   @param dbar   Set the mean particle diameter (m)
     */
    double meanParticleDiameter = 1.5E-6;
    tranDusty->setMeanParticleDiameter(meanParticleDiameter);



    tranDusty->getMultiDiffCoeffs(nsp, DATA_PTR(multiD));

    printf("MultiDiffusion coefficients: \n");
    for (int i = 0; i < nsp; i++) {
      for (int j = 0; j < nsp; j++) {
	printf(" %15.8E,", multiD[nsp*j + i]);
      }
      printf("\n");
    }

    //double visc = tranDusty->viscosity();
    //printf("visc = %g\n", visc);

    // tranDusty->getSpeciesViscosities(x);
    //for (k = 0; k < nsp; k++) {
    // string spName = g->speciesName(k);
    // printf("sp visc (%s) = %g\n", spName.c_str(), x[k]);
    //}


    // double cond = tranDusty->thermalConductivity();
    //printf("cond = %g\n", cond);



    // tranDusty->getMixDiffCoeffs(x);
    //for (k = 0; k < nsp; k++) {
    // string spName = g->speciesName(k);
    // printf("sp diff (%s) = %g\n", spName.c_str(), x[k]);
    //}

    //  tranDusty->getMobilities(x);
    // for (k = 0; k < nsp; k++) {
    // string spName = g->speciesName(k);
    // printf("Mobility (%s) = %g\n", spName.c_str(), x[k]);
    //}

    //   double gradX[20];
    // for (k = 0; k < 20; k++) {
    //  gradX[k] = 0.0;
    // }
    //gradX[1] = 1.0;
    //double gradT = 0.0;

    // tranDusty->getSpeciesFluxes(1, &gradT, 5, &gradX[0], 5, x);
    //for (k = 0; k < nsp; k++) {
    // string spName = g->speciesName(k);
    // printf("SpeciesFlux (%s) = %g\n", spName.c_str(), x[k]);
    //}

   
    //   gradX[1] = 0.0;
    //double gradV = 1.0;

    //  tranDusty->set_Grad_T(&gradT);
    // tranDusty->set_Grad_V(&gradV);
    // tranDusty->set_Grad_X(gradX);

    //   tranDusty->getSpeciesFluxesExt(5, x);
    //  for (k = 0; k < nsp; k++) {
    // string spName = g->speciesName(k);
    // printf("SpeciesFlux (%s) = %g\n", spName.c_str(), x[k]);
    // }



    vector_fp state1;
    g->saveState(state1);
    g->setState_TP(Temp, 1.2 * OneAtm);
    vector_fp state2;
    g->saveState(state2);
    double delta = 0.001;
    vector_fp fluxes;
    fluxes.resize(nsp);

  //! Get the molar fluxes [kmol/m^2/s], given the thermodynamic state at two nearby points. 
    /*!
     *
     *   \f[
     *       J_k = - \sum_{j = 1, N} \left[D^{multi}_{kj}\right]^{-1} \left( \nabla C_j  + \frac{C_j}{\mathcal{D}^{knud}_j} \frac{\kappa}{\mu} \nabla p \right) 
     *   \f]
     *
     * @param  state1  Array of temperature, density, and mass fractions for state 1.
     * @param  state2  Array of temperature, density, and mass  fractions for state 2.  
     * @param  delta   Distance from state 1 to state 2 (m).
     *
     * @param fluxes   Vector of species molar fluxes due to diffusional driving force 
     */
    tranDusty->getMolarFluxes(DATA_PTR(state1), DATA_PTR(state1), delta, DATA_PTR(fluxes));
    for (int i = 0; i < nsp; i++) {
      printf(" flux[%d] = %13.8E\n", i, fluxes[i]);
    }
    tranDusty->getMolarFluxes(DATA_PTR(state1), DATA_PTR(state2), delta, DATA_PTR(fluxes));
    for (int i = 0; i < nsp; i++) {
      printf(" flux[%d] = %13.8E\n", i, fluxes[i]);
    }


    delete g;
    g=0;
    delete tran;
    tran = 0;

    Cantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
