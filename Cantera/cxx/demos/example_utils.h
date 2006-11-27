#ifndef CT_EXAMPLE_UTILS_H
#define CT_EXAMPLE_UTILS_H

#include <cantera/kernel/Array.h>
#include <cantera/kernel/plots.h>

namespace Cantera{}
using namespace Cantera;
namespace std{}
using namespace std;
namespace CanteraZeroD{}
using namespace CanteraZeroD;

// Save the temperature, density, pressure, and mole fractions at one
// time
template<class G, class A>
void saveSoln(int i, double time, const G& gas, A& soln) {  
    soln(0,i) = time;
    soln(1,i) = gas.temperature();
    soln(2,i) = gas.density();
    soln(3,i) = gas.pressure();
    gas.getMoleFractions(&soln(4,i));
}

template<class G, class A>
void saveSoln(double time, const G& gas, A& soln) {  
    soln.resize(soln.nRows(), soln.nColumns() + 1);
    int back = soln.nColumns() - 1;
    soln(0,back) = time;
    soln(1,back) = gas.temperature();
    soln(2,back) = gas.density();
    soln(3,back) = gas.pressure();
    int nsp = gas.nSpecies();
    for (int k = 0; k < nsp; k++) 
        soln(4+k,back) = gas.moleFraction(k);
}

template<class G, class V>
void makeDataLabels(const G& gas, V& names) {
    int nsp = gas.nSpecies();
    names.resize(nsp + 4);
    names[0] = "time (s)";
    names[1]  = "Temperature (K)";
    names[2]  = "Density (kg/m3)";
    names[3]  = "Pressure (Pa)";
    int k;
    for (k = 0; k < nsp; k++) names[4+k] = gas.speciesName(k);
}

template<class G, class A>
void plotSoln(string fname, string fmt, string title, const G& gas, const A& soln) {
    vector<string> names;
    makeDataLabels(gas, names);
    writePlotFile(fname, fmt, title, names, soln);
}

inline void writeCanteraHeader(ostream& s) {
    s << endl;
    s << "     Cantera version " << "CANTERA_VERSION" << endl;         
    s << "     Copyright California Institute of Technology, 2002." << endl;
    s << "     http://www.cantera.org" << endl;
    s << endl;
}

#endif
