// Example utilities
// =================
//
// Defines some functions used by :doc:`kinetics1.cpp <kinetics1>`.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EXAMPLE_UTILS_H
#define CT_EXAMPLE_UTILS_H

#include "cantera/base/Array.h"
#include <fstream>

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

void writeCsv(const std::string& fname, const Cantera::ThermoPhase& gas,
              const Cantera::Array2D& data)
{
    std::ofstream s(fname);
    // Write labels
    s << "time (s),Temperature (K),Density (kg/m3),Pressure (Pa),";

    for (size_t k = 0; k < gas.nSpecies(); k++) {
        s << gas.speciesName(k);
        if (k != gas.nSpecies()-1) {
            s << ",";
        }
    }
    s << std::endl;
    for (size_t i = 0; i < data.nColumns(); i++) {
        for (size_t j = 0; j < data.nRows(); j++) {
            s << data(j,i);
            if (j != data.nRows()-1) {
                s << ",";
            }
        }
        s << std::endl;
    }
}

#endif
