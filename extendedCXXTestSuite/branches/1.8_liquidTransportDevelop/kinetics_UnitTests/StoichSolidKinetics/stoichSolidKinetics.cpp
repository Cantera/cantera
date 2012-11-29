/**
 *  @file 
 *
 */

//  Example 
//
// Note that this example needs updating. It works fine, but is
// written in a way that is less than transparent or
// user-friendly. This could be rewritten using class Interface to
// make things simpler.


#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

#ifdef DEBUG_HKM
int iDebug_HKM = 0;
#endif

/*****************************************************************/
/*****************************************************************/

#include "cantera/Cantera.h"
#include "cantera/electrolyteThermo.h"
#include "cantera/kinetics.h"
#include "cantera/kernel/MargulesVPSSTP.h"
#include "cantera/kernel/IonsFromNeutralVPSSTP.h"
#include "cantera/IdealGasMix.h"
#include "cantera/kernel/PrintCtrl.h"


using namespace Cantera;

void printDbl(double val) {
  if (fabs(val) < 5.0E-17) {
    cout << " nil";
  } else {
    cout << val;
  }
}

bool ctversionCanHandle(std::string ctv) {
  if (ctv.substr(0,2) == "1.7") {
    return false;
  } 
  if (ctv == "1.8.0") {
    return false;
  }
  if (ctv == "1.8_liquidTransportDevelop") {
    return true;
  }
  if (ctv == "1.8.x") {
    return true;
  }
  if (ctv.substr(0,2) == "1.8") {
    return false;
  }
  printf("ctversionCanHandle: error: unknown version: %s\n", ctv.c_str());
  exit(-1);
  return false;
}


int main(int argc, char** argv) {
  int i;
  //double ca[10];

  std::string ctv = CANTERA_VERSION;
  bool canHandle = ctversionCanHandle(ctv);

  try {
    string id;
    Cantera::PrintCtrl pc;

    XML_Node *xc = new XML_Node();
    string path = findInputFile("ReactionSurf.xml");
    ctml::get_CTML_Tree(xc, path);

    XML_Node * const xg = xc->findNameID("phase", "reaction_surface");
    if (!xg) {
      throw CanteraError("couldn't find file", "");
    }
    ThermoPhase *surfTP = newPhase(*xg);

  
    ThermoPhase *gasTP =  newPhase("gas.xml", "");
    //IdealGasPhase *gasTP = new IdealGasPhase("gas.xml", "");

    string s_cao_s = "CaO_Solid.xml";
    id = "CaO(S)";
    Cantera::ThermoPhase *cao_s = Cantera::newPhase(s_cao_s, id);

    string  s_caco3_s = "CaCO3_Solid.xml";
    id = "CaCO3(S)";
    Cantera::ThermoPhase *caco3_s = Cantera::newPhase(s_caco3_s, id);

    string  s_c_s = "C_Solid.xml";
    id = "C(S)";
    Cantera::ThermoPhase *c_s = Cantera::newPhase(s_c_s, id);

    string  s_fe3o4_s = "Fe3O4_Solid.xml";
    id = "Fe3O4(S)";
    Cantera::ThermoPhase *fe3o4_s = Cantera::newPhase(s_fe3o4_s, id);

    string  s_feo_s = "FeO_Solid.xml";
    id = "FeO(S)";
    Cantera::ThermoPhase *feo_s = Cantera::newPhase(s_feo_s, id);


    string  s_fe_s = "Fe_Solid.xml";
    id = "Fe(S)";
    Cantera::ThermoPhase *fe_s = Cantera::newPhase(s_fe_s, id);

    //   string s_reactSurf = "ReactionSurf.xml";
    // id = "reaction_surface";
    //Cantera::ThermoPhase *surfTP = Cantera::newPhase(s_reactSurf, id);


    vector<ThermoPhase *> phaseList;
     
    phaseList.push_back(gasTP);
    phaseList.push_back(cao_s);
    phaseList.push_back(caco3_s);
    phaseList.push_back(c_s);
    phaseList.push_back(fe3o4_s);
    phaseList.push_back(feo_s);
    phaseList.push_back(fe_s);
    phaseList.push_back(surfTP);
    InterfaceKinetics *iKin_ptr = new InterfaceKinetics();
    importKinetics(*xg, phaseList, iKin_ptr);
    //    int nr = iKin_ptr->nReactions();
   

    double x[20];
    for (i = 0; i < 20; i++) x[i] = 0.0;
    double mll[20];
    for (i = 0; i < 20; i++) mll[i] = 0.0;
    //int ikca0 = iKin_ptr->kineticsSpeciesIndex("CaO(S)");
    //int ikcaco3 = iKin_ptr->kineticsSpeciesIndex("CaCO3(S)");
    //int ikco2 = iKin_ptr->kineticsSpeciesIndex("CO2");


    int igco2 = gasTP->speciesIndex("CO2");
    int igo2  = gasTP->speciesIndex("O2");
    int ign2  = gasTP->speciesIndex("N2");

    mll[igco2] = 0.2;
    mll[igo2] = 0.1;
    mll[ign2] = 0.7;
    double Temp = 1000.;
    gasTP->setState_TPX(Temp, OneAtm, mll);
   
    double p = 20.0*OneAtm/760.0;
    
    for (i = 0; i < 20; i++) x[i] = 0.0;
   

    for (i = 0; i < 20; i++) x[i] = 0.0;
    x[0] = 1.0;


    // Set the bath gas of TempK and 1 atm
    cao_s->setState_TP(Temp, OneAtm);
    caco3_s->setState_TP(Temp, OneAtm);
    c_s->setState_TP(Temp, OneAtm);
    fe3o4_s->setState_TP(Temp, OneAtm);
    surfTP->setState_TP(Temp, OneAtm);

    // Throw some asserts in here to test that they compile
    AssertTrace(p == p); 
    AssertThrow(p == p, "main"); 
    AssertThrowMsg(i == 20, "main", "are you kidding"); 

    double src[20];
    for (i = 0; i < 20; i++) src[i] = 0.0;
    //iKin_ptr->getNetProductionRates(src);
    double deltaG[20], deltaGSS[20];

    // Preamble
    printf("==========================================================================\n");
    printf("   StoichSolidKinetics Test \n");
    printf("==========================================================================\n");

    printf("  Tests for the proper behavior of heterogeneous reactions \n");
    printf("  when phases may or may not exist:\n");
    printf("        CaCO3(s) =   CO2(g) +  CaO(s)  \n");
    printf("\n\n");

    for (int ktrials = 0; ktrials < 2; ktrials++) {

    iKin_ptr->getDeltaSSGibbs(deltaGSS);

    cout << "   deltaGSS      = ";
    pc.pr_de(deltaGSS[0], 4, 13);
    cout << "\n";

    iKin_ptr->getDeltaGibbs(deltaG);
    cout << "   deltaG        = ";
    pc.pr_de(deltaG[0], 4, 13);
    cout << "\n";

    double mu[10];
    gasTP->getChemPotentials(mu);
    double mu_CO2 = mu[igco2];
    cout << "   mu_CO2(g)     = ";
    pc.pr_de(mu_CO2, 4, 13);
    cout << "\n";

 
    cao_s->getGibbs_RT(mu);
    double mu_cao = mu[0] * GasConstant * Temp;
    cout << "   mu_cao(s)     = ";
    pc.pr_de(mu_cao, 4, 13);
    cout << "\n";

    caco3_s->getChemPotentials(mu);
    double mu_caco3  = mu[0];
    cout << "   mu_caco3      = ";
    pc.pr_de(mu_caco3, 4, 13);
    cout << "\n";

    double deltaG_calc =  mu_CO2 +  mu_cao - mu_caco3;
    cout << "   deltaG_calc   = ";
    pc.pr_de(deltaG_calc, 4, 13);
    cout << "\n";
 
    double act[20];
    gasTP->getActivities(act);
    double act_CO2 = act[igco2];
    cout << "   act_CO2       = ";
    pc.pr_de(act_CO2, 4, 13);
    cout << "\n";

    cao_s->getActivities(act);
    cout << "   act_cao(s)    = ";
    pc.pr_de(act[0], 4, 13);
    cout << "\n";

    caco3_s->getActivities(act);
    cout << "   act_caco3(s)  = ";
    pc.pr_de(act[0], 4, 13);
    cout << "\n";

 
    vector<double> Rfwd(5);
    vector<double> Rrev(5);
    vector<double> Rnet(5);
    printf("==========================================================================\n");
    printf("  OK base problem assuming that all phases exist:\n");
    printf("==========================================================================\n");

  
    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;

    double kfwd[5];
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

    double krev[5];
    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;


#if  defined (CANTERA_VERSION_18_LTD) || defined (CANTERA_VERSION_18_XXX)

    printf("==========================================================================\n");
    printf("  OK Changing the problem setting CaO(S) phase to nonexistent:\n");
    printf("==========================================================================\n");

    int ip_cao = iKin_ptr->phaseIndex("CaO(S)");
    iKin_ptr->setPhaseExistence(ip_cao, false);
    iKin_ptr->setPhaseStability(ip_cao, true);

    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;


   
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;

    printf("==========================================================================\n");
    printf("  OK Changing the problem setting CaCO3(S) phase to nonexistent:\n");
    printf("==========================================================================\n");

    int ip_caco3 = iKin_ptr->phaseIndex("CaCO3(S)");
    iKin_ptr->setPhaseExistence(ip_caco3, false);
    iKin_ptr->setPhaseStability(ip_caco3, true);


    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;

   
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

 
    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;


    printf("==========================================================================\n");
    printf("  OK Changing the problem setting CaO(S) phase to existent, CaCO3 nonexistent:\n");
    printf("==========================================================================\n");


    iKin_ptr->setPhaseExistence(ip_cao, true);

    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;

   
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

 
    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;


    printf("==========================================================================\n");
    printf("  OK Changing the problem setting Gas phase to nonexistent, CaCO3 nonexistent:\n");
    printf("==========================================================================\n");

    int ip_gas = iKin_ptr->phaseIndex("air");
    iKin_ptr->setPhaseExistence(ip_gas, false);
    iKin_ptr->setPhaseStability(ip_gas, true);

    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;

   
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

 
    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;


    printf("==========================================================================\n");
    printf("  OK Changing the problem setting to all phases existing:\n");
    printf("==========================================================================\n");
    iKin_ptr->setPhaseExistence(ip_gas, true);
    iKin_ptr->setPhaseExistence(ip_cao, true);
    iKin_ptr->setPhaseExistence(ip_caco3, true);


    iKin_ptr->getNetRatesOfProgress(DATA_PTR(Rnet));
    cout << "ROP_net       = ";
    pc.pr_de(Rnet[0], 4, 13);
    cout << endl;

  
    iKin_ptr->getFwdRatesOfProgress(DATA_PTR(Rfwd));
    cout << "ROP_forward   = ";
    pc.pr_de(Rfwd[0], 4, 13);
    cout << endl;

    iKin_ptr->getRevRatesOfProgress(DATA_PTR(Rrev));
    cout << "ROP_reverse   = ";
    pc.pr_de(Rrev[0], 4, 13);
    cout << endl;

   
    iKin_ptr->getFwdRateConstants(DATA_PTR(kfwd));
    cout << "    kfwd      = ";
    pc.pr_de(kfwd[0], 4, 13);
    cout << endl;

 
    iKin_ptr->getRevRateConstants(DATA_PTR(krev));
    cout << "    krev      = ";
    pc.pr_de(krev[0], 4, 13);
    cout << endl;

#endif

    if (ktrials == 0) {
      printf("==========================================================================\n");
      printf("  OK Changing the problem setting so that forward rate if faster:\n");
      printf("==========================================================================\n");
      mll[igco2] = 0.002;
      mll[igo2] =  0.1;
      mll[ign2] =  0.898;
      gasTP->setState_TPX(Temp, OneAtm, mll);
    }

    }


    printf("==========================================================================\n");

  }
  catch (CanteraError) {
    showErrors(cout);
  }

  return 0;
}
/***********************************************************/
