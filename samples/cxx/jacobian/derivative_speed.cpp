/*!
 * @file derivative_speed.cpp
 *
 * Benchmark tests for derivative evaluations
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
void timeit_base(void (Kinetics::*function)(double*),
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

//! timer for Eigen::VectorXd getters
void timeit_vector(Eigen::VectorXd (Kinetics::*function)(),
                   Kinetics* kin,
                   ThermoPhase& gas,
                   size_t loops=10000,
                   size_t runs=7)
{
    Eigen::VectorXd out;

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
    std::cout << "getFwdRatesOfProgress:  ";
    timeit_base(&Kinetics::getFwdRatesOfProgress, &kin, gas, nReactions);

    std::cout << "fwdRatesOfProgress_ddT: ";
    timeit_vector(&Kinetics::fwdRatesOfProgress_ddT, &kin, gas);

    std::cout << "fwdRatesOfProgress_ddP: ";
    timeit_vector(&Kinetics::fwdRatesOfProgress_ddP, &kin, gas);

    std::cout << "fwdRatesOfProgress_ddC: ";
    timeit_vector(&Kinetics::fwdRatesOfProgress_ddC, &kin, gas);

    std::cout << "fwdRatesOfProgress_ddX: ";
    timeit_matrix(&Kinetics::fwdRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getRevRatesOfProgress:  ";
    timeit_base(&Kinetics::getRevRatesOfProgress, &kin, gas, nReactions);

    std::cout << "revRatesOfProgress_ddT: ";
    timeit_vector(&Kinetics::revRatesOfProgress_ddT, &kin, gas);

    std::cout << "revRatesOfProgress_ddP: ";
    timeit_vector(&Kinetics::revRatesOfProgress_ddP, &kin, gas);

    std::cout << "revRatesOfProgress_ddC: ";
    timeit_vector(&Kinetics::revRatesOfProgress_ddC, &kin, gas);

    std::cout << "revRatesOfProgress_ddX: ";
    timeit_matrix(&Kinetics::revRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getNetRatesOfProgress:  ";
    timeit_base(&Kinetics::getNetRatesOfProgress, &kin, gas, nReactions);

    std::cout << "netRatesOfProgress_ddT: ";
    timeit_vector(&Kinetics::netRatesOfProgress_ddT, &kin, gas);

    std::cout << "netRatesOfProgress_ddP: ";
    timeit_vector(&Kinetics::netRatesOfProgress_ddP, &kin, gas);

    std::cout << "netRatesOfProgress_ddC: ";
    timeit_vector(&Kinetics::netRatesOfProgress_ddC, &kin, gas);

    std::cout << "netRatesOfProgress_ddX: ";
    timeit_matrix(&Kinetics::netRatesOfProgress_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getCreationRates:       ";
    timeit_base(&Kinetics::getCreationRates, &kin, gas, nSpecies);

    std::cout << "creationRates_ddT:      ";
    timeit_vector(&Kinetics::creationRates_ddT, &kin, gas);

    std::cout << "creationRates_ddP:      ";
    timeit_vector(&Kinetics::creationRates_ddP, &kin, gas);

    std::cout << "creationRates_ddC:      ";
    timeit_vector(&Kinetics::creationRates_ddC, &kin, gas);

    std::cout << "creationRates_ddX:      ";
    timeit_matrix(&Kinetics::creationRates_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getDestructionRates:    ";
    timeit_base(&Kinetics::getDestructionRates, &kin, gas, nSpecies);

    std::cout << "destructionRates_ddT:   ";
    timeit_vector(&Kinetics::destructionRates_ddT, &kin, gas);

    std::cout << "destructionRates_ddP:   ";
    timeit_vector(&Kinetics::destructionRates_ddP, &kin, gas);

    std::cout << "destructionRates_ddC:   ";
    timeit_vector(&Kinetics::destructionRates_ddC, &kin, gas);

    std::cout << "destructionRates_ddX:   ";
    timeit_matrix(&Kinetics::destructionRates_ddX, &kin, gas);

    std::cout << std::endl;
    std::cout << "getNetProductionRates:  ";
    timeit_base(&Kinetics::getNetProductionRates, &kin, gas, nSpecies);

    std::cout << "netProductionRates_ddT: ";
    timeit_vector(&Kinetics::netProductionRates_ddT, &kin, gas);

    std::cout << "netProductionRates_ddP: ";
    timeit_vector(&Kinetics::netProductionRates_ddP, &kin, gas);

    std::cout << "netProductionRates_ddC: ";
    timeit_vector(&Kinetics::netProductionRates_ddC, &kin, gas);

    std::cout << "netProductionRates_ddX: ";
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
