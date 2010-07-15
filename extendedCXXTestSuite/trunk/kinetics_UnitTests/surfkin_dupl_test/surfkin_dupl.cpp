/*
 * $Id: surfdemo.cpp 255 2009-11-09 23:36:49Z hkmoffa $
 *
 *  Sample program that solves an implicit problem for surface
 *  site fractions.
 */
#include "cantera/Cantera.h"
#include <iostream>

#include "cantera/IdealGasMix.h"
#include "cantera/Interface.h"

using namespace Cantera;
using namespace Cantera_CXX;
using namespace std;
using namespace Cantera_CXX;

int main() {

  try {
    IdealGasMix gas("gri30.xml", "gri30");
    gas.setState_TPX(1200.0, OneAtm, 
		     "H2:2, O2:1, OH:0.01, H:0.01, O:0.01");

    vector<ThermoPhase*> phases;
    phases.push_back(&gas);
    Interface surf("surface.xml", "surface", phases);
    vector_fp cov;
    cov.push_back(0.8);
    cov.push_back(0.2);
    surf.setCoverages(DATA_PTR(cov));
    vector_fp wdot(gas.nSpecies() + surf.nSpecies());
    surf.getNetProductionRates(DATA_PTR(wdot));
    int k;
    for (k = 0; k < gas.nSpecies(); k++)
      cout << gas.speciesName(k) << "  " << wdot[k] << endl;

    for (k = 0; k < surf.nSpecies(); k++)
      cout << surf.speciesName(k) << "   " 
	   << wdot[k+gas.nSpecies()] << endl;


    int nph = surf.nPhases();
    std::vector<ThermoPhase *> thermoVec;
    for (int iph = 0; iph < nph; iph++) {
      ThermoPhase * tp_ph = & surf.thermo(iph);
      ThermoPhase *tp2_ph = tp_ph->duplMyselfAsThermoPhase();
      thermoVec.push_back(tp2_ph);  
    }
    
    InterfaceKinetics *ik = dynamic_cast<InterfaceKinetics *>(&surf);

    Kinetics *ksurf2 = ik->duplMyselfAsKinetics(thermoVec);
    Cantera::InterfaceKinetics *surf2 = dynamic_cast<InterfaceKinetics *>(ksurf2);
    
    int isurf = surf2->surfacePhaseIndex();
    ThermoPhase *tsurp = thermoVec[isurf];
    SurfPhase *surp = dynamic_cast<SurfPhase *>(tsurp);

    surp->setCoverages(DATA_PTR(cov));
    surf2->getNetProductionRates(DATA_PTR(wdot));
 
    for (k = 0; k < gas.nSpecies(); k++)
      cout << gas.speciesName(k) << "  " << wdot[k] << endl;

    for (k = 0; k < surf.nSpecies(); k++) {
      cout << surp->speciesName(k) << "   " 
	   << wdot[k+gas.nSpecies()] << endl;
    }

    std::vector<ThermoPhase *> thermoVec3;
    for (int iph = 0; iph < nph; iph++) {
      ThermoPhase * tp_ph = thermoVec[iph];
      ThermoPhase *tp3_ph = tp_ph->duplMyselfAsThermoPhase();
      thermoVec3.push_back(tp3_ph);  
    }



    Kinetics *ksurf3 = ksurf2->duplMyselfAsKinetics(thermoVec3);
    Cantera::InterfaceKinetics *surf3 = dynamic_cast<InterfaceKinetics *>(ksurf3);

    delete ksurf2;
    for (int iph = 0; iph < nph; iph++) {
      delete  thermoVec[iph];
    }

    isurf = surf3->surfacePhaseIndex();
    tsurp = thermoVec3[isurf];
    surp = dynamic_cast<SurfPhase *>(tsurp);

    surp->setCoverages(DATA_PTR(cov));
    surf3->getNetProductionRates(DATA_PTR(wdot));
 

    for (k = 0; k < gas.nSpecies(); k++)
      cout << gas.speciesName(k) << "  " << wdot[k] << endl;

    for (k = 0; k < surf.nSpecies(); k++) {
      cout << surp->speciesName(k) << "   " 
	   << wdot[k+gas.nSpecies()] << endl;
    }
    delete ksurf3;
 
    for (int iph = 0; iph < nph; iph++) {
      delete  thermoVec3[iph];
    }
    Cantera::appdelete();

  }
  catch (CanteraError) {
    showErrors(cout);
  }
  return 0;
}
      
