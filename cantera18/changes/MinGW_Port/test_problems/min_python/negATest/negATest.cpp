
/*
 *  $Author: dggoodwin $
 *  $Date: 2007/05/04 15:18:45 $
 *  $Revision: 1.3 $
 *

 */
#ifdef SRCDIRTREE
#include "ct_defs.h"
#include "ctml.h"
#include "GasKinetics.h"
#include "importCTML.h"
#include "ThermoPhase.h"
#include "InterfaceKinetics.h"
#else
#include "Cantera.h"
#include "kinetics.h"
#endif

using namespace Cantera;
using namespace std;

#
int main(int argc, char **argv) {
  try {   
    int i;
    string infile = "noxNeg.cti";
    double x[20];
    double cdot[20], ddot[20];
    
    XML_Node *xc = new XML_Node();
    string path = findInputFile(infile);
    ctml::get_CTML_Tree(xc, path);

    XML_Node * const xg = xc->findNameID("phase", "air");
    ThermoPhase *gasTP = newPhase(*xg);
    int nsp = gasTP->nSpecies();
    cout << "Number of species = " << nsp << endl;


    vector<ThermoPhase *> phaseList;
    phaseList.push_back(gasTP);
    GasKinetics *iKin_ptr = new GasKinetics();
    importKinetics(*xg, phaseList, iKin_ptr);
    int nr = iKin_ptr->nReactions();
    cout << "Number of reactions = " << nr << endl;

    int iH2 = gasTP->speciesIndex("H2");
    int iH = gasTP->speciesIndex("H");
    int iO2 = gasTP->speciesIndex("O2");
    int iOH = gasTP->speciesIndex("OH");
    int iH2O = gasTP->speciesIndex("H2O");
    int iNH = gasTP->speciesIndex("NH");
    int iNO = gasTP->speciesIndex("NO");
    int iN2O = gasTP->speciesIndex("N2O");

    for (i = 0; i < nsp; i++) {
      x[i] = 0.0;
    }
    x[iH2O] = 1.0 /2.0;
    x[iOH]  = 0.1 /2.0;
    x[iH]   = 0.2 /2.0;
    x[iO2]  = 0.3 /2.0;
    x[iH2]  = 0.25 /2.0;
    x[iNH]  = 0.05/2.0;
    x[iNO]  = 0.05/2.0;
    x[iN2O]  = 0.05/2.0;

    double p = OneAtm;

    gasTP->setState_TPX(2000., p, x);


    double src[20];
    for (i = 0; i < 20; i++) src[i] = 0.0;
    iKin_ptr->getNetProductionRates(src);

    for (i = 0; i < nsp; i++) {
      string sSt = gasTP->speciesName(i);
      printf("rop [ %d:%s ] = %g \n", i, sSt.c_str(), src[i]);
    }

    int nReactions = iKin_ptr->nReactions();
    cout << "number of reactions = " << nReactions << endl;

    double fwd_rop[20];
    double rev_rop[20];
    iKin_ptr->getFwdRatesOfProgress(fwd_rop);
    iKin_ptr->getRevRatesOfProgress(rev_rop);
    for (i = 0; i < nReactions; i++) {
      printf("fwd_rop[%3d] = %13g    rev_rop[%3d] = %13g\n", i, fwd_rop[i],
	     i, rev_rop[i]);
    }



    iKin_ptr->getCreationRates(cdot);
    iKin_ptr->getDestructionRates(ddot);
    

    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
