/**
 *  @file example2.cpp
 *
 */

//  Example 
//
//  Read a mechanism and a thermodynamics file for the 
//  class IdealSolidSolnPhase in order to test that it's
//  working correctly
//

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
/*****************************************************************/
static void printUsage()
{

}


#include "kernel/ct_defs.h"
#include "kernel/ctml.h"
#include "kernel/GasKinetics.h"
#include "kernel/importCTML.h"
#include "kernel/ThermoPhase.h"
#include "kernel/InterfaceKinetics.h"

using namespace Cantera;

int main(int argc, char** argv) {
    int i, k;
    string infile = "diamond.cti";
 
    try {
      XML_Node *xc = new XML_Node();
      string path = findInputFile(infile);
      ctml::get_CTML_Tree(xc, path);

      XML_Node *xg = findXMLPhase(xc, "gas");
      ThermoPhase *gasTP = newPhase(*xg);
      int nsp = gasTP->nSpecies();
      cout << "Number of species = " << nsp << endl;

      XML_Node *xd = findXMLPhase(xc, "diamond");
      ThermoPhase *diamondTP = newPhase(*xd);
      int nsp_diamond = diamondTP->nSpecies();
      cout << "Number of species in diamond = " << nsp_diamond << endl;


      XML_Node *xs = findXMLPhase(xc, "diamond_100");
      ThermoPhase *diamond100TP = newPhase(*xs);
      int nsp_d100 = diamond100TP->nSpecies();
      cout << "Number of species in diamond_100 = " << nsp_d100 << endl;

      vector<ThermoPhase *> phaseList;
      phaseList.push_back(gasTP);     
      phaseList.push_back(diamondTP);
      phaseList.push_back(diamond100TP);
      InterfaceKinetics *iKin_ptr = new InterfaceKinetics();
      importKinetics(*xs, phaseList, iKin_ptr);
      int nr = iKin_ptr->nReactions();
      cout << "Number of reactions = " << nr << endl;

      double x[20];
      for (i = 0; i < 20; i++) x[i] = 0.0;
      x[0] = 0.0010;
      x[1] = 0.9888;
      x[2] = 0.0002;
      x[3] = 0.0100;
      double p = 20.0*OneAtm/760.0;

      gasTP->setState_TPX(1200., p, x);
      
      for (i = 0; i < 20; i++) x[i] = 0.0;
      int i0 = diamond100TP->speciesIndex("c6H*");
      x[i0] = 0.1;
      int i1 = diamond100TP->speciesIndex("c6HH");
      x[i1] = 0.9;
      diamond100TP->setState_TX(1200., x);

      for (i = 0; i < 20; i++) x[i] = 0.0;
      x[0] = 1.0;
      diamondTP->setState_TPX(1200., p, x);

      iKin_ptr->advanceCoverages(100.);

   
      double src[20];
      for (i = 0; i < 20; i++) src[i] = 0.0;
      iKin_ptr->getNetProductionRates(src);
      double sum = 0.0;
      double naH;
      for (k = 0; k < 13; k++) {
	if (k < 4) {
	  naH = gasTP->nAtoms(k, 0);
	} else if (k == 4) {
	  naH = 0;
	} else if (k > 4) {
	  int itp = k - 5;
	  naH = diamond100TP->nAtoms(itp, 0);
	}
	cout << k << "  " << naH << "  " ;
	if (fabs(src[k]) < 2.0E-17) {
	  cout << " nil" << endl;
	} else {
	  cout << src[k] << endl;
	}
	sum += naH * src[k];
	
      }
  
      cout << "sum = " << sum << endl;
      double mwd = diamondTP->molecularWeight(0);
      double dens = diamondTP->density();
      double gr = src[4] * mwd / dens;
      gr *= 1.0E6 * 3600.;
      cout << "growth rate = " << gr << " microns per hour" << endl;

      
      diamond100TP->getMoleFractions(x);
      cout << "Coverages:" << endl;
      for (k = 0; k < 8; k++) {
	int iii = 3;
	cout << k << "   " << diamond100TP->speciesName(k) 
	     << "   " 
	     << x[k] << endl;
      }    
    }
    catch (CanteraError) {
      showErrors(cout);
    }

    return 0;
}
/***********************************************************/
