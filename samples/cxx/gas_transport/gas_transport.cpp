/*!
 * @file gas_transport.cpp
 *
 * Gas phase transport properties
 *
 * Construct a gas phase Solution object and use it to compute viscosity,
 * thermal conductivity, mixture-averaged diffusion coefficients, and thermal
 * diffusivities for a range of temperatures.
 *
 * Keywords: tutorial, transport, saving output
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Array.h"
#include "cantera/base/plots.h"

#include <iostream>

using namespace Cantera;
using std::cout;
using std::endl;

void write_csv(const std::string& name, const std::vector<std::string>& names,
               const Array2D& data)
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
    auto sol = newSolution("gri30.yaml", "gri30", "Mix");
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
        sol->transport()->getMixDiffCoeffs(&output(3,i));
    }

    // Create a list of labels for the CSV output file
    std::vector<std::string> labels {
        "Temperature (K)",
        "Viscosity (Pa*s)",
        "Thermal Conductivity (W/m*K)"
    };
    for (size_t k = 0; k < gas->nSpecies(); k++) {
        labels.push_back(gas->speciesName(k));
    }

    // Save transport properties to a file
    write_csv("transport_mix.csv", labels, output);

    // Create a new transport manager for multicomponent properties
    unique_ptr<Transport> multi(
        newTransportMgr("multicomponent", sol->thermo().get()));

    // Get multicomponent properties at several temperatures
    for (int i = 0; i < ntemps; i++) {
        temp = 500.0 + 100.0*i;
        gas->setState_TP(temp, pres);
        output(0,i) = temp;
        output(1,i) = multi->viscosity();
        output(2,i) = multi->thermalConductivity();
        multi->getThermalDiffCoeffs(&output(3,i));
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
