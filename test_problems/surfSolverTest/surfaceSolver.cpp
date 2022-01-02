/**
 *  @file surfaceSolver.cpp
 *
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/solveSP.h"
#include "cantera/base/fmt.h"
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Cantera;

#define MSSIZE 200

void printGas(ostream& oooo, ThermoPhase* gasTP, InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    double C[MSSIZE];
    oooo.precision(3);
    string gasPhaseName = gasTP->name();
    gasTP->getMoleFractions(x);
    gasTP->getConcentrations(C);
    double Temp = gasTP->temperature();
    double p = gasTP->pressure();
    oooo << "Gas Temperature = " << Temp << endl;
    oooo << "Gas Pressure    = " << p << endl;
    size_t iPhase = iKin_ptr->phaseIndex(gasPhaseName);
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, iPhase);
    oooo << "Gas Phase:  " << gasPhaseName << "   "
         << "(" << kstart << ")" << endl;
    oooo << "                       Name      "
         << "     Conc              MoleF       SrcRate " << endl;
    oooo << "                                 "
         << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
    double sum = 0.0;
    size_t nspGas = gasTP->nSpecies();
    for (size_t k = 0; k < nspGas; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(k, iPhase);
        fmt::print(oooo, "{:4d} {:>24s}   {:14.3g} {:14.3g}  {:14.3e}\n",
                   k, gasTP->speciesName(k), C[k], x[k], src[kstart]);
        sum += x[k];
    }
    oooo << "Sum of gas mole fractions= " << sum << endl;
    oooo << endl;
}

void printBulk(ostream& oooo,
               ThermoPhase* bulkPhaseTP, InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    double C[MSSIZE];
    oooo.precision(3);
    string bulkParticlePhaseName = bulkPhaseTP->name();
    bulkPhaseTP->getMoleFractions(x);
    bulkPhaseTP->getConcentrations(C);
    size_t iPhase = iKin_ptr->phaseIndex(bulkParticlePhaseName);
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, iPhase);
    double dens = bulkPhaseTP->density();
    oooo << "Bulk Phase:  " << bulkParticlePhaseName << "   "
         << "(" << kstart << ")" << endl;
    double Temp = bulkPhaseTP->temperature();
    double p = bulkPhaseTP->pressure();
    oooo << "Bulk Temperature = " << Temp << endl;
    oooo << "Bulk Pressure    = " << p << endl;
    oooo << "                       Name      "
         << "     Conc              MoleF       SrcRate " << endl;
    oooo << "                                 "
         << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
    double sum = 0.0;
    double Wsum = 0.0;
    const vector_fp& molecW = bulkPhaseTP->molecularWeights();
    size_t nspBulk = bulkPhaseTP->nSpecies();
    for (size_t k = 0; k < nspBulk; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(k, iPhase);
        fmt::print(oooo, "{:4d} {:>24s}   {:14.3g} {:14.3g}  {:14.3e}\n",
                   k, bulkPhaseTP->speciesName(k), C[k], x[k], src[kstart]);
        sum += x[k];
        Wsum += src[kstart] * molecW[k];
    }
    oooo.precision(3);
    oooo << "Bulk Weight Growth Rate = " << Wsum << " kg/m^2/s" << endl;
    double gr = Wsum / dens;
    oooo << "Bulk Growth Rate = " << gr << " m/s" << endl;
    oooo << "Bulk Growth Rate = " << gr * 1.0E6 * 3600.
         << " microns / hour" << endl;
    oooo << "Density of bulk phase = " << dens << " kg / m^3 "<< endl;
    oooo << "                      = " << dens / 1.0E3
         <<" gm / cm^3 " << endl;
    oooo << "Sum of bulk mole fractions= " << sum << endl;
    oooo << endl;
}

void printSurf(ostream& oooo,
               ThermoPhase* surfPhaseTP, InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    string surfParticlePhaseName = surfPhaseTP->name();
    surfPhaseTP->getMoleFractions(x);
    size_t iPhase = iKin_ptr->phaseIndex(surfParticlePhaseName);
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, iPhase);
    oooo << "Surface Phase:  " << surfParticlePhaseName
         << " (" << kstart << ")" << endl;
    double Temp = surfPhaseTP->temperature();
    double p = surfPhaseTP->pressure();
    oooo << "Surface Temperature = " << Temp << endl;
    oooo << "Surface Pressure    = " << p << endl;
    oooo << "                       Name      "
         << "   Coverage         SrcRate " << endl;
    double sum = 0.0;
    size_t nspSurf = surfPhaseTP->nSpecies();
    for (size_t k = 0; k < nspSurf; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(0, iPhase);
        double srcK = src[kstart];
        if (fabs(srcK) < 1.0E-7) {
            srcK = 0.0;
        }
        fmt::print(oooo, "{:4d} {:>24s}   {:14.3g}   {:14.3e}\n",
                   k, surfPhaseTP->speciesName(k), x[k], srcK);
        sum += x[k];
    }
    oooo << "Sum of coverages = " << sum << endl;
}

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    string infile = "haca2.yaml";
    string gasPhaseName = "gas";
    string bulkParticlePhaseName = "soot";
    string surfParticlePhaseName = "soot_interface";
    int ioflag = 1;

    try {
        auto iface = newInterface(infile, surfParticlePhaseName);
        ThermoPhase* gasTP = iface->adjacent(gasPhaseName)->thermo().get();
        size_t nspGas = gasTP->nSpecies();
        cout << "Number of species = " << nspGas << endl;

        ThermoPhase* bulkPhaseTP = iface->adjacent(bulkParticlePhaseName)->thermo().get();
        size_t nspBulk = bulkPhaseTP->nSpecies();
        cout << "Number of species in bulk phase named " <<
             bulkParticlePhaseName << " = " << nspBulk << endl;

        ThermoPhase* surfPhaseTP = iface->thermo().get();
        size_t nsp_d100 = surfPhaseTP->nSpecies();
        cout << "Number of species in surface phase, " << surfParticlePhaseName
             << " = " << nsp_d100 << endl;

        InterfaceKinetics* iKin_ptr = iface->kinetics().get();
        size_t nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        double x[MSSIZE];

        ofstream ofile("results.txt");

        iKin_ptr->setIOFlag(ioflag);
        /*
         *  Solve the Equation system
         */
        iKin_ptr->solvePseudoSteadyStateProblem();

        /*
         * Download the source terms for the rate equations
         */
        double src[MSSIZE];
        iKin_ptr->getNetProductionRates(src);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;

        printGas(ofile, gasTP, iKin_ptr, src);
        printBulk(ofile, bulkPhaseTP, iKin_ptr, src);
        printSurf(ofile, surfPhaseTP, iKin_ptr, src) ;
        /*****************************************************************************/
        /*  Now Tweak the inputs and do a quick calculation */
        /****************************************************************************/

        /*
         * Set the Gas State:
         * -> note that the states are set in the input file too
         */
        double pres = gasTP->pressure();
        gasTP->getMoleFractions(x);
        double tmp = 0.3 * std::min(x[0], x[1]);
        x[0] += tmp;
        x[1] -= tmp;
        gasTP->setState_PX(pres, x);

        iKin_ptr->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;

        printGas(ofile, gasTP, iKin_ptr, src);
        printBulk(ofile, bulkPhaseTP, iKin_ptr, src);
        printSurf(ofile, surfPhaseTP, iKin_ptr, src) ;
        /*****************************************************************************/
        /*  Now Tweak the inputs and do a quick calculation */
        /****************************************************************************/

        pres = surfPhaseTP->pressure();
        double temp = surfPhaseTP->temperature();
        temp += 95;
        surfPhaseTP->setState_TP(temp, pres);

        iKin_ptr->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;

        printGas(ofile, gasTP, iKin_ptr, src);
        printBulk(ofile, bulkPhaseTP, iKin_ptr, src);
        printSurf(ofile, surfPhaseTP, iKin_ptr, src) ;

        /*****************************************************************************/
        /*  Now Don't Tweak the inputs at all */
        /****************************************************************************/
        surfPhaseTP->setState_TP(temp, pres);

        iKin_ptr->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;

        printGas(ofile, gasTP, iKin_ptr, src);
        printBulk(ofile, bulkPhaseTP, iKin_ptr, src);
        printSurf(ofile, surfPhaseTP, iKin_ptr, src) ;
        appdelete();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }

    return 0;
}
/***********************************************************/
