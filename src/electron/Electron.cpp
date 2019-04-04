// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/Electron.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/numerics/funcs.h"
#include <iostream>
#include <limits>
#include <set>

namespace Cantera {

Electron::Electron()
    : m_electronCrossSectionTargets(0)
    , m_electronCrossSectionKinds(0)
    , m_ncs(0)
    , m_points(1000)
    , m_kTe(Undef)
    , m_kT(Undef)
    , m_electronCrossSections_ok(false)
    , m_f0_ok(false)
{
    // default energy grid
    m_eps.resize(m_points);
    for (size_t j = 0; j < m_points; j++) {
        m_eps[j] = j / 100.0;
    }

    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

Electron::~Electron()
{
}

void Electron::init(thermo_t* thermo)
{
    m_thermo = thermo;
    m_f0_ok = false;
}

void Electron::update_T()
{
    // signal that temperature-dependent quantities will need to be recomputed
    // before use, and update the local temperature.
    m_kT = Boltzmann * m_thermo->temperature() / ElectronCharge;

    // flag for quantities need to be re-calculated
    m_f0_ok = false;
}

void Electron::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    compositionMap gas_composition = m_thermo->getMoleFractionsByName(0.0);
    m_moleFractions.resize(m_ncs, 0.0);
    for (auto const& x : gas_composition) {
        bool not_found = true;
        for (size_t i = 0; i < m_ncs; i++) {
            if (m_electronCrossSectionTargets[i] == x.first) {
                m_moleFractions[i] = x.second;
                not_found = false;
            }
        }
        if (not_found) {
            if (x.second > 0.01) {
                std::cout << "The mole fraction of species " << x.first
                            << " is more than 0.01 but it has no cross section data."
                            << std::endl;
            }
        }
    }
    // flag for quantities need to be re-calculated
    m_f0_ok = false;
    m_totalCrossSection_ok = false;
}

bool Electron::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    if (ecs->kind == "EFFECTIVE") {
        for (size_t i = 0; i < m_ncs; i++) {
            if (m_electronCrossSectionTargets[i] == ecs->target) {
                if (m_electronCrossSectionKinds[i] == "EFFECTIVE") {
                    throw CanteraError("Electron::addElectronCrossSection",
                                        "Already contains a data of EFFECTIVE cross section for '{}'.",
                                        ecs->target);
                }
            }
        }
    }
    ecs->validate();
    m_electronCrossSectionTargets.push_back(ecs->target);
    m_electronCrossSectionKinds.push_back(ecs->kind);
    m_massRatios.push_back(ecs->mass_ratio);

    // transpose data
    std::vector<std::vector<double>> transdata(2, std::vector<double>(ecs->data.size()));
    for (size_t i = 0; i < ecs->data.size(); i++) {
        for (size_t j = 0; j < 2; j++) {
            transdata[j][i] = ecs->data[i][j];
        }
    }
    m_electronCrossSectionData.push_back(transdata);
    m_ncs++;
    m_electronCrossSections_ok = false;
    return true;
}

void Electron::setupGrid(size_t n, const double* eps)
{
    m_points = n;
    m_eps.resize(n);
    for (size_t j = 0; j < m_points; j++) {
        m_eps[j] = eps[j];
    }
    m_electronCrossSections_ok = false;
    m_totalCrossSection_ok = false;
    m_f0_ok = false;
}

void Electron::setupCrossSections()
{
    m_electronCrossSections.resize(m_ncs, std::vector<double>(m_points));
    for (size_t i = 0; i < m_ncs; i++) {
        vector_fp x = m_electronCrossSectionData[i][0];
        vector_fp y = m_electronCrossSectionData[i][1];
        if (x[0] > 0.0) {
            x.insert(x.begin(), 0.0);
            y.insert(y.begin(), m_electronCrossSectionData[i][1][0]);
        }
        x.push_back(1e8);
        y.push_back(m_electronCrossSectionData[i][1].back());
        for (size_t j = 0; j < m_points; j++) {
            m_electronCrossSections[i][j] = linearInterp(m_eps[j], x, y);
        }
    }
    m_electronCrossSections_ok = true;
}

void Electron::calculateTotalCrossSection()
{
    if (m_electronCrossSections_ok == false) {
        setupCrossSections();
    }
    m_totalCrossSection.resize(m_points, 0.0);
    m_attachCrossSection.resize(m_points, 0.0);
    m_ionizCrossSection.resize(m_points, 0.0);
    for (size_t j = 0; j < m_points; j++) {
        for (size_t i = 0; i < m_ncs; i++) {
            if (m_electronCrossSectionKinds[i] == "EFFECTIVE") {
                m_totalCrossSection[j] += m_moleFractions[i] * m_electronCrossSections[i][j];
            } else if (m_electronCrossSectionKinds[i] == "ATTACHMENT") {
                m_attachCrossSection[j] += m_moleFractions[i] * m_electronCrossSections[i][j];
            } else if (m_electronCrossSectionKinds[i] == "IONIZATION") {
                m_ionizCrossSection[j] += m_moleFractions[i] * m_electronCrossSections[i][j];
            }
        }
    }
    m_totalCrossSection_ok = true;
}

double Electron::netProductionFrequency(const vector_fp& f0)
{
    if (m_totalCrossSection_ok == false) {
        calculateTotalCrossSection();
    }
    double vi = 0.0;
    if (f0.size() != m_points) {
        throw CanteraError("Electron::netProductionFrequency",
                            "The size of input vector must equal to grid points, {}.",
                            m_points);
    }
    for (size_t j = 0; j < m_points - 1; j++) {
        double left = (m_ionizCrossSection[j] - m_attachCrossSection[j]) * m_eps[j] * f0[j];
        double right = (m_ionizCrossSection[j+1] - m_attachCrossSection[j+1]) * m_eps[j+1] * f0[j+1];
        vi += 0.5 * (left + right) * (m_eps[j+1] - m_eps[j]);
    }
    return vi;
}

void Electron::calculateDistributionFunction()
{
    update_T();
    update_C();
    calculateTotalCrossSection();
    m_f0.resize(m_points);
    m_df0.resize(m_points);
    m_f0_ok = true;
}

double Electron::electronDiffusivity(double N)
{
    if (m_f0_ok == false) {
        calculateDistributionFunction();
    }
    double vi = netProductionFrequency(m_f0);
    vector_fp y(m_points, 0.0);
    for (size_t j = 0; j < m_points; j++) {
        if (m_eps[j] != 0.0) {
            y[j] = m_eps[j] * m_f0[j] / (m_totalCrossSection[j] + vi / pow(m_eps[j], 0.5));
        }
    }
    return 1./3. * m_gamma * simpsonQuadrature(m_eps, y) / N;
}

double Electron::electronMobility(double N)
{
    if (m_f0_ok == false) {
        calculateDistributionFunction();
    }
    double vi = netProductionFrequency(m_f0);
    vector_fp y(m_points, 0.0);
    for (size_t j = 0; j < m_points; j++) {
        if (m_eps[j] != 0.0) {
            y[j] = m_eps[j] * m_df0[j] / (m_totalCrossSection[j] + vi / pow(m_eps[j], 0.5));
        }
    }
    return -1./3. * m_gamma * simpsonQuadrature(m_eps, y) / N;
}

}
