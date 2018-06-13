#ifndef CT_EXAMPLE_UTILS_H
#define CT_EXAMPLE_UTILS_H

#include "cantera/base/Array.h"
#include "cantera/base/plots.h"

// Save the temperature, density, pressure, and mole fractions at one
// time
template<class G, class A>
void saveSoln(int i, double time, const G& gas, A& soln)
{
    soln(0,i) = time;
    soln(1,i) = gas.temperature();
    soln(2,i) = gas.density();
    soln(3,i) = gas.pressure();
    gas.getMoleFractions(&soln(4,i));
}

template<class G, class A>
void saveSoln(double time, const G& gas, A& soln)
{
    soln.resize(soln.nRows(), soln.nColumns() + 1);
    int back = soln.nColumns() - 1;
    soln(0,back) = time;
    soln(1,back) = gas.temperature();
    soln(2,back) = gas.density();
    soln(3,back) = gas.pressure();
    int nsp = gas.nSpecies();
    for (int k = 0; k < nsp; k++) {
        soln(4+k,back) = gas.moleFraction(k);
    }
}

template<class G, class V>
void makeDataLabels(const G& gas, V& names)
{
    int nsp = gas.nSpecies();
    names.resize(nsp + 4);
    names[0] = "time (s)";
    names[1] = "Temperature (K)";
    names[2] = "Density (kg/m3)";
    names[3] = "Pressure (Pa)";
    int k;
    for (k = 0; k < nsp; k++) {
        names[4+k] = gas.speciesName(k);
    }
}

template<class G, class A>
void plotSoln(const std::string& fname, const std::string& fmt,
              const std::string& title, const G& gas, const A& soln)
{
    std::vector<std::string> names;
    makeDataLabels(gas, names);
    writePlotFile(fname, fmt, title, names, soln);
}

#endif
