/*!
 * @file derivative_speed.cpp
 *
 * Benchmark derivative evaluations
 *
 * Calculate derivatives of reaction rates of progress and species production and
 * destruction rates with respect to temperature, pressure, molar concentration,
 * and mole fractions. Time evaluation for different chemical mechanisms.
 *
 * Keywords: kinetics, benchmarking
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <chrono>
#include <iostream>
#include <iomanip>
#include <numeric>
#include "cantera/base/Solution.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics.h"
#include "cantera/numerics/eigen_sparse.h"

using namespace Cantera;

void statistics(vector_fp times, size_t loops, size_t runs)
{
    double average = accumulate(times.begin(), times.end(), 0.0) / times.size();
    for (auto& v : times) {
        v = (v - average) * (v - average);
    }
    double std = accumulate(times.begin(), times.end(), 0.0) / times.size();
    std = pow(std, 0.5);

    // output statistics
    std::cout << std::setprecision(5) << average / 1000. << " μs ± "
        << std::setprecision(3) << std / 1000. << " μs "
        << "per loop (" << runs << " runs, " << loops << " loops each)\n";
}

//! timer for standard getters
void timeit_array(void (Kinetics::*function)(double*),
                  Kinetics* kin,
                  ThermoPhase& gas,
                  size_t siz,
                  size_t loops=10000,
                  size_t runs=7)
{
    vector_fp out(siz);

    double T = gas.temperature();
    double pressure = gas.pressure();
    double deltaT = 1e-5;

    vector_fp times;
    for (size_t run = 0; run < runs; ++run) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0.; i < loops; ++i) {
            gas.setState_TP(T + i * deltaT, pressure);
            (kin->*function)(out.data());
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        times.push_back(
            std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /
            loops);
    }
    gas.setState_TP(T, pressure);
    statistics(times, loops, runs);
}

//! timer for Eigen::SparseMatrix<double> getters
void timeit_matrix(Eigen::SparseMatrix<double> (Kinetics::*function)(),
                   Kinetics* kin,
                   ThermoPhase& gas,
                   size_t loops=1000,
                   size_t runs=7)
{
    Eigen::SparseMatrix<double> out;

    double T = gas.temperature();
    double pressure = gas.pressure();
    double deltaT = 1e-5;

    vector_fp times;
    for (size_t run = 0; run < runs; ++run) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0.; i < loops; ++i) {
            gas.setState_TP(T + i * deltaT, pressure);
            out = (kin->*function)();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        times.push_back(
            std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /
            loops);
    }
    gas.setState_TP(T, pressure);
    statistics(times, loops, runs);
}

void benchmark(const std::string& mech, const std::string& phase,
    const std::string& fuel)
{
    auto sol = newSolution(mech, phase, "None");
    auto& gas = *(sol->thermo());
    auto& kin = *(sol->kinetics());

    auto nSpecies = gas.nSpecies();
    auto nReactions = kin.nReactions();
    std::cout << mech << ": "
        << nSpecies << " species, " << nReactions << " reactions." << std::endl;

    double Tu = 300.0;
    double pressure = 101325.0;
    gas.setState_TP(Tu, pressure);
    gas.setEquivalenceRatio(1.0, fuel, "O2:1, N2:3.76");
    gas.equilibrate("HP");

    std::cout << std::endl;
    std::cout << "getFwdRatesOfProgress:     ";
    timeit_array(&Kinetics::getFwdRatesOfProgress, &kin, gas, nReactions);

    std::cout << "getFwdRatesOfProgress_ddT: ";
    timeit_array(&Kinetics::getFwdRatesOfProgress_ddT, &kin, gas, nReactions);

    std::cout << "getFwdRatesOfProgress_ddP: ";
    timeit_array(&Kinetics::getFwdRatesOfProgress_ddP, &kin, gas, nReactions);

    std::cout << "getFwdRatesOfProgress_ddC: ";
    timeit_array(&Kinetics::getFwdRatesOfProgress_ddC, &kin, gas, nReactions);

    std::cout << "fwdRatesOfProgress_ddX:    ";
    timeit_matrix(&Kinetics::fwdRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getRevRatesOfProgress:     ";
    timeit_array(&Kinetics::getRevRatesOfProgress, &kin, gas, nReactions);

    std::cout << "getRevRatesOfProgress_ddT: ";
    timeit_array(&Kinetics::getRevRatesOfProgress_ddT, &kin, gas, nReactions);

    std::cout << "getRevRatesOfProgress_ddP: ";
    timeit_array(&Kinetics::getRevRatesOfProgress_ddP, &kin, gas, nReactions);

    std::cout << "getRevRatesOfProgress_ddC: ";
    timeit_array(&Kinetics::getRevRatesOfProgress_ddC, &kin, gas, nReactions);

    std::cout << "revRatesOfProgress_ddX:    ";
    timeit_matrix(&Kinetics::revRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getNetRatesOfProgress:     ";
    timeit_array(&Kinetics::getNetRatesOfProgress, &kin, gas, nReactions);

    std::cout << "getNetRatesOfProgress_ddT: ";
    timeit_array(&Kinetics::getNetRatesOfProgress_ddT, &kin, gas, nReactions);

    std::cout << "getNetRatesOfProgress_ddP: ";
    timeit_array(&Kinetics::getNetRatesOfProgress_ddP, &kin, gas, nReactions);

    std::cout << "getNetRatesOfProgress_ddC: ";
    timeit_array(&Kinetics::getNetRatesOfProgress_ddC, &kin, gas, nReactions);

    std::cout << "netRatesOfProgress_ddX:    ";
    timeit_matrix(&Kinetics::netRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getCreationRates:          ";
    timeit_array(&Kinetics::getCreationRates, &kin, gas, nSpecies);

    std::cout << "getCreationRates_ddT:      ";
    timeit_array(&Kinetics::getCreationRates_ddT, &kin, gas, nSpecies);

    std::cout << "getCreationRates_ddP:      ";
    timeit_array(&Kinetics::getCreationRates_ddP, &kin, gas, nSpecies);

    std::cout << "getCreationRates_ddC:      ";
    timeit_array(&Kinetics::getCreationRates_ddC, &kin, gas, nSpecies);

    std::cout << "creationRates_ddX:         ";
    timeit_matrix(&Kinetics::creationRates_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getDestructionRates:       ";
    timeit_array(&Kinetics::getDestructionRates, &kin, gas, nSpecies);

    std::cout << "getDestructionRates_ddT:   ";
    timeit_array(&Kinetics::getDestructionRates_ddT, &kin, gas, nSpecies);

    std::cout << "getDestructionRates_ddP:   ";
    timeit_array(&Kinetics::getDestructionRates_ddP, &kin, gas, nSpecies);

    std::cout << "getDestructionRates_ddC:   ";
    timeit_array(&Kinetics::getDestructionRates_ddC, &kin, gas, nSpecies);

    std::cout << "destructionRates_ddX:      ";
    timeit_matrix(&Kinetics::destructionRates_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getNetProductionRates:     ";
    timeit_array(&Kinetics::getNetProductionRates, &kin, gas, nSpecies);

    std::cout << "getNetProductionRates_ddT: ";
    timeit_array(&Kinetics::getNetProductionRates_ddT, &kin, gas, nSpecies);

    std::cout << "getNetProductionRates_ddP: ";
    timeit_array(&Kinetics::getNetProductionRates_ddP, &kin, gas, nSpecies);

    std::cout << "getNetProductionRates_ddC: ";
    timeit_array(&Kinetics::getNetProductionRates_ddC, &kin, gas, nSpecies);

    std::cout << "netProductionRates_ddX:    ";
    timeit_matrix(&Kinetics::netProductionRates_ddX, &kin, gas);
}

int main()
{
    std::cout << "Benchmark tests for derivative evaluations." << std::endl;
    std::cout << std::endl;
    benchmark("h2o2.yaml", "ohmech", "H2");
    std::cout << std::endl;
    benchmark("gri30.yaml", "gri30", "CH4");
    std::cout << std::endl;
    benchmark("nDodecane_Reitz.yaml", "nDodecane_IG", "c12h26");
    return 0;
}
