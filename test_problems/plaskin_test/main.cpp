
#include <cmath>
#include <vector>
#include <iostream>
#include "cantera/IdealGasMix.h"  
#include "cantera/equilibrium.h"   
#include "cantera/transport.h"    
#include <cantera/Cantera.h>
#include "cantera/zerodim.h"

using namespace std;
using namespace Cantera;

int main()
{

  /* declarations */
  cout.precision(6);
  cout.setf(ios::scientific);

  IdealGasMix gas("mech.cti","mymech");
  string  file = "source";
  gas.writeMech(file);
  
  int     nr   = gas.nReactions();
  int     nrte = gas.nElectronImpactReactions();
  int     nrvr = gas.nVibRelaxationReactions();
  int     ns   = gas.nSpecies();
  int     ne   = gas.nElements();

  double  T      = 1200.0;
  double  Te     = OneEV;
  double  p      = OneAtm;
  double* x      = new double[ns];
  double* mw     = new double[ns];
  double* deltaE = new double[nrte];
  double* kfwd   = new double[nr];
  double* krev   = new double[nr];
  double* omey   = new double[ns];
  int*    teindx = new int[nrte];

  cout << endl;
  cout << "nr   = " << nr << endl;
  cout << "nrte = " << nrte << endl;
  cout << "nrvr = " << nrvr << endl;
  cout << "ns   = " << ns << endl;
  cout << "ne   = " << ne << endl;

  for(int i = 0; i < ns; ++i){
    x[i] = 1.0;
  }
  
  gas.setState_TeTPX(Te, T, p, x);
  gas.getMolecularWeights(mw);
  gas.getDeltaEPlasma(deltaE);
  gas.getTeDepIndices(teindx);

  cout << endl;
  cout << "r  = " << gas.density() << endl;
  cout << "T  = " << gas.temperature() << endl;
  cout << "Te = " << gas.elec_temperature() << endl;
  cout << endl;

  for(int i = 0; i < ns; ++i) {
    cout << "m_speciesIndex[" << "\"" << gas.speciesName(i).c_str() << "\"" << "] = " << i << ";" << endl;
  }
  cout << endl;

  for(int i = 0; i < ns; ++i) {
    cout << "m_mw[" << i << "] = " << mw[i] << ";" << endl;
  }
  cout << endl;

  for(int i = 0; i < ns; ++i) {
    cout << "m_charge[" << i << "] = " << gas.charge(i) << ";" << endl;
  }
  cout << endl;

  for(int i = 0; i < nrte; ++i) {
    cout << "m_eirIndex[" << i << "] = " << teindx[i] << ";" << endl;
  }
  cout << endl;

  for(int i = 0; i < nrte; ++i) {
    cout << "m_dHinel[" << i << "] = " << deltaE[i] << ";" << endl;
  }
  cout << endl;

  // gas.writeMech(file);
   
}

