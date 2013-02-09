/* ======================================================================= */
/* $RCSfile: latsol.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2011/03/10 23:58:53 $ */
/* $Revision: 1.1 $ */
/* ======================================================================= */

#include <stdio.h>


#include "cantera/Cantera.h"
#include "cantera/thermo.h"
#include "cantera/kernel/MargulesVPSSTP.h"
#include "cantera/kernel/FixedChemPotSSTP.h"
#include "cantera/kernel/LatticeSolidPhase.h"


#include "cantera/equilibrium.h"

#include "cantera/thermo.h"
#include "cantera/kernel/vcs_internal.h"
#include "cantera/kernel/logger.h"




using namespace std;
using namespace Cantera;

int CHECK_DEBUG_MODE = 0;

void printUsage() {
    cout << "usage: HMW_test_1 " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}


int main(int argc, char **argv)
{

   int retn = 0;

   try {
     double x[20];
     for (int k = 0; k < 20; k++) {
       x[k] = 0.0;
     }
     double um[20];
     VCSnonideal::vcs_timing_print_lvl = 0;

     string f_li7si3 = "Li7Si3_ls.xml";
     string id = "Li7Si3_and_Interstitials(S)";
#ifdef CANTERA_HAS_TEMPLATES
     Cantera::ThermoPhase<doublereal> *LiSi_solid = Cantera::newPhase<doublereal>(f_li7si3, id);
#else
     Cantera::ThermoPhase *LiSi_solid = Cantera::newPhase(f_li7si3, id);
#endif

     MargulesVPSSTP *salt = new MargulesVPSSTP(1);

     //int iKCl_l = salt->speciesIndex("KCl(L)");
     //int iLiCl_l = salt->speciesIndex("LiCl(L)");

#ifdef CANTERA_HAS_TEMPLATES
     ThermoPhase<doublereal> *Li_liq = newPhase<doublereal>("Li_Liquid.xml","Li(L)");
#else
     ThermoPhase *Li_liq = newPhase("Li_Liquid.xml","Li(L)");
#endif

     /*
      * set states
      */
     x[0] = 0.7;
     x[1] = 1.0 - x[0];
     double T = 273.15 + 352.;
     salt->setState_TPX(T, OneAtm, x);
     LiSi_solid->setState_TP(T, OneAtm);

     salt->getChemPotentials(um);

     LiSi_solid->getChemPotentials(um);

     int ee =  LiSi_solid->nElements();
     printf("Number of elements = %d\n", ee);




    FixedChemPotSSTP *LiFixed = new FixedChemPotSSTP("Li", -2.3E7);


    LiFixed->setState_TP(T, OneAtm);

    LiFixed->getChemPotentials(um);
    printf(" chem pot = %g\n", um[0]);


    double  volts = 1.635;     // has some Fe in it // test suite

    double dg_corr =  - volts * Faraday;
    printf("dg_corr = %g\n", dg_corr);

    Li_liq->getChemPotentials(um);
    double um_Li_liq = um[0];

    double um_li_chempot = um_Li_liq + dg_corr;

    printf("um_li_chempot = %g\n", um_li_chempot);
    LiFixed->setChemicalPotential(um_li_chempot);




    Cantera::MultiPhase mmm;

    mmm.addPhase(salt, 10.);
    mmm.addPhase(LiSi_solid, 1.);
    mmm.addPhase(LiFixed, 100.);

    int   solver = 2;
    int printLvl = 10;
    int estimateEquil = 0;

    vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    cout << mmm << endl;


    delete LiSi_solid;
    LiSi_solid = 0;
    delete Li_liq;
    Li_liq = 0;

    delete LiFixed;
    LiFixed = 0;

     delete salt;
     salt = 0;
     Cantera::appdelete();

     return retn;

   } catch (CanteraError) {

     showErrors();
     return -1;
   }
} 
