/*
 * LiC6 electrode
 * ==============
 *
 * Calculate the open-circuit potential of an LiC6 electrode and activity
 * coefficients of each species as a function of the mole fraction of
 * intercalated lithium.
 *
 * .. tags:: C++, thermodynamics, electrochemistry, battery
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/core.h"
#include <iostream>
#include <fstream>

using namespace Cantera;

void calc_potentials()
{
    double Tk = 273.15 + 25.0;

    string filename = "LiC6_electrodebulk.yaml";
    string phasename = "LiC6_and_Vacancies";
    auto sol = newSolution(filename, phasename);
    auto electrodebulk = sol->thermo();
    string intercalatingSpeciesName("Li(C6)");
    size_t intercalatingSpeciesIdx = electrodebulk->speciesIndex(intercalatingSpeciesName);
    size_t nsp_tot = electrodebulk->nSpecies();

    std::ofstream fout("LiC6_electrode_output.csv", std::ofstream::out);
    fout << "x[LiC6], ChemPotential[LiC6], ChemPotential[C6], Uref, ActCoeff[LiC6], ActCoeff[C6], dlnActCoeffdx[LiC6], dlnActCoeffdx[C6]" << std::endl;

    vector<double> spvals(nsp_tot);
    vector<double> actCoeff(nsp_tot);
    vector<double> dlnActCoeffdlnX_diag(nsp_tot);
    double xmin = 0.6;
    double xmax = 0.9;

    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    size_t nsp_electrodeBulk = electrodebulk->nSpecies();
    vector<double> xv(nsp_electrodeBulk, 0.0);

    for (int i = 0; i < numSteps; ++i) {
        double x = xmin + i*dx;

        vector<double> xv(nsp_electrodeBulk, 0.0);
        // Set the fraction of intercalated lithium
        xv[intercalatingSpeciesIdx] = x;

        //Set so that mole fractions sum to 1
        for (size_t j = 0; j < nsp_electrodeBulk; ++j) {
            if (j != intercalatingSpeciesIdx) {
                xv[j] = (1.0 - xv[intercalatingSpeciesIdx]);
            }
        }

        electrodebulk->setMoleFractions(xv);
        electrodebulk->setTemperature(Tk);
        electrodebulk->getChemPotentials(spvals);

        // Calculate the open circuit potential
        double Uref = (spvals[1] - spvals[0])/Faraday;

        electrodebulk->getdlnActCoeffdlnX_diag(dlnActCoeffdlnX_diag);
        electrodebulk->getActivityCoefficients(actCoeff);

        fout << fmt::format("{}, {}, {}, {}, {}, {}, {}, {}\n",
            xv[0], spvals[0], spvals[1], Uref, actCoeff[0],
            actCoeff[1], dlnActCoeffdlnX_diag[0], dlnActCoeffdlnX_diag[1]);
    }
}

int main(int argc, char** argv)
{
    try {
        calc_potentials();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
