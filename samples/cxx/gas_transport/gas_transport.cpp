/*
 * Gas phase transport properties
 * ==============================
 *
 * Construct a gas phase Solution object and use it to compute viscosity,
 * thermal conductivity, mixture-averaged diffusion coefficients, and thermal
 * diffusivities for a range of temperatures.
 *
 * .. tags:: C++, tutorial, transport, saving output
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/core.h"
#include "cantera/base/Array.h"

#include <iostream>
#include <fstream>

using namespace Cantera;
using std::cout;
using std::endl;

void write_csv(const string& name, const vector<string>& names, const Array2D& data)
{
    std::ofstream s(name);
    for (size_t i = 0; i < data.nRows(); i++) {
        s << names[i];
        if (i != data.nRows()-1) {
            s << ",";
        }
    }
    s << endl;
    for (size_t i = 0; i < data.nColumns(); i++) {
        for (size_t j = 0; j < data.nRows(); j++) {
            s << data(j,i);
            if (j != data.nRows()-1) {
                s << ",";
            }
        }
        s << endl;
    }
}

void transport_example()
{
    // create a gas mixture, and set its state
    auto sol = newSolution("gri30.yaml", "gri30", "mixture-averaged");
    auto gas = sol->thermo();
    double temp = 500.0;
    double pres = 2.0*OneAtm;
    gas->setState_TPX(temp, pres, "H2:1.0, CH4:0.1");

    size_t nsp = gas->nSpecies();

    // create a 2D array to hold the outputs
    int ntemps = 20;
    Array2D output(nsp+3, ntemps);

    // Get mixture-averaged properties at several temperatures
    for (int i = 0; i < ntemps; i++) {
        temp = 500.0 + 100.0*i;
        gas->setState_TP(temp, pres);
        output(0,i) = temp;
        output(1,i) = sol->transport()->viscosity();
        output(2,i) = sol->transport()->thermalConductivity();
        sol->transport()->getMixDiffCoeffs(span<double>(&output(3,i), nsp));
    }

    // Create a list of labels for the CSV output file
    vector<string> labels {
        "Temperature (K)",
        "Viscosity (Pa*s)",
        "Thermal Conductivity (W/m*K)"
    };
    for (size_t k = 0; k < gas->nSpecies(); k++) {
        labels.push_back(gas->speciesName(k));
    }

    // Save transport properties to a file
    write_csv("transport_mix.csv", labels, output);

    // Switch transport manager to multicomponent properties
    sol->setTransportModel("multicomponent");

    // Get multicomponent properties at several temperatures
    for (int i = 0; i < ntemps; i++) {
        temp = 500.0 + 100.0*i;
        gas->setState_TP(temp, pres);
        output(0,i) = temp;
        output(1,i) = sol->transport()->viscosity();
        output(2,i) = sol->transport()->thermalConductivity();
        sol->transport()->getThermalDiffCoeffs(span<double>(&output(3,i), nsp));
    }

    // Save transport properties to a file
    write_csv("transport_multi.csv", labels, output);

    cout << "Output files:" << endl
            << "  transport_mix.csv" << endl
            << "  transport_multi.csv" << endl;
}

int main(int argc, char** argv) {
    try {
        transport_example();
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        return -1;
    }
    return 0;
}
