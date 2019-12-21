// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

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
    : m_ncs(0)
    , m_points(200)
    , m_kT(Undef)
    , m_E(0.0)
    , m_F(0.0)
    , m_f0_ok(false)
    , m_maxn(100)
    , m_rtol(1e-5)
    , m_delta0(1e14)
    , m_factorM(4.0)
    , m_init_kTe(0.0)
    , m_warn(true)
{
    // default energy grid
    m_gridCenter.resize(m_points);
    m_gridEdge.resize(m_points + 1);
    m_f0.resize(m_points);
    for (size_t j = 0; j < m_points; j++) {
        m_gridCenter[j] = j / 20.0 + 1.0 / 40.0;
        m_gridEdge[j] = j / 20.0;
    }
    m_gridEdge[m_points] = 10.0;
    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

Electron::~Electron()
{
}

void Electron::init(thermo_t* thermo)
{
    m_thermo = thermo;
    m_f0_ok = false;
    calculateElasticCrossSection();
    setGridCache();
    // set up target index
    m_kTargets.resize(m_ncs);
    for (size_t k = 0; k < m_ncs; k++) {
        m_kTargets[k] = m_thermo->speciesIndex(m_targets[k]);
    }
    // set up indices of species which has no cross-section data
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        auto it = std::find(m_kTargets.begin(), m_kTargets.end(), k); 
        if (it == m_kTargets.end()) {
            m_kOthers.push_back(k);
        }
    }
    m_kProducts.clear();
    m_kProducts.resize(m_ncs);
    for (size_t k = 0; k < m_ncs; k++) {
        std::string prdct = m_products[k];
        size_t dp = prdct.find(" + ");
        if (dp != npos) {
            size_t kp0 = m_thermo->speciesIndex(prdct.substr(0, dp));
            prdct.erase(0, dp + 3);
            size_t kp1 = m_thermo->speciesIndex(prdct);
            m_kProducts[k].push_back(kp0);
            m_kProducts[k].push_back(kp1);
        } else {
            size_t kp0 = m_thermo->speciesIndex(prdct);
            m_kProducts[k].push_back(kp0);
            m_kProducts[k].push_back(9999);
        }
    }
    m_moleFractions.resize(m_ncs, 0.0);
}

void Electron::setGridCache()
{
    m_sigma.clear();
    m_sigma.resize(m_ncs);
    m_eps.clear();
    m_eps.resize(m_ncs);
    m_j.clear();
    m_j.resize(m_ncs);
    m_i.clear();
    m_i.resize(m_ncs);
    for (size_t k = 0; k < m_ncs; k++) {
        vector_fp& x = m_crossSections[k][0];
        vector_fp& y = m_crossSections[k][1];
        vector_fp eps1(m_points + 1);
        for (size_t i = 0; i < m_points + 1; i++) {
            eps1[i] = clip(m_shiftFactor[k] * m_gridEdge[i] + m_thresholds[k],
                           m_gridEdge[0] + 1e-9, m_gridEdge[m_points] - 1e-9);
        }

        vector_fp nodes = eps1;
        for (size_t i = 0; i < m_points + 1; i++) {
            if (m_gridEdge[i] >= eps1[0] && m_gridEdge[i] <= eps1[m_points]) {
                nodes.push_back(m_gridEdge[i]);
            }
        }
        for (size_t i = 0; i < x.size(); i++) {
            if (x[i] >= eps1[0] && x[i] <= eps1[m_points]) {
                nodes.push_back(x[i]);
            }
        }

        std::sort(nodes.begin(), nodes.end());
        auto last = std::unique(nodes.begin(), nodes.end());
        nodes.resize(std::distance(nodes.begin(), last));
        vector_fp sigma0(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            sigma0[i] = linearInterp(nodes[i], x, y);
        }

        // search position of cell j
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(m_gridEdge.begin(), m_gridEdge.end(), nodes[i]);
            m_j[k].push_back(low - m_gridEdge.begin() - 1);
        }

        // search position of cell i
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(eps1.begin(), eps1.end(), nodes[i]);
            m_i[k].push_back(low - eps1.begin() - 1);
        }

        // construct sigma
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp sigma{sigma0[i], sigma0[i+1]};
            m_sigma[k].push_back(sigma);
        }

        // construct eps
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp eps{nodes[i], nodes[i+1]};
            m_eps[k].push_back(eps);
        }
    }
}

void Electron::update_T()
{
    // signal that temperature-dependent quantities will need to be recomputed
    // before use, and update the local temperature.
    double kT = Boltzmann * m_thermo->temperature() / ElectronCharge;
    if (m_kT != kT) {
        m_kT = kT;
        m_N = m_thermo->pressure() / Boltzmann / m_thermo->temperature();
        m_f0_ok = false;
    }
}

void Electron::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    vector_fp X(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&X[0]);
    for (size_t k = 0; k < m_ncs; k++) {
        size_t kk = m_kTargets[k];
        if (m_moleFractions[k] != X[kk]) {
            m_moleFractions[k] = X[kk];
            m_f0_ok = false;
        }
    }
    // warn that a specific species needs cross-section data.
    for (size_t k : m_kOthers) {
        if (X[k] > 0.01) {
            writelog("Cantera::Electron::update_C");
            writelog("\n");
            writelog("Warning: The mole fraction of species {} is more than 0.01",
                    m_thermo->speciesName(k));
            writelog(" but it has no data of cross section.");
            writelog("\n");
        }
    }
}

bool Electron::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    ecs->validate();
    m_targets.push_back(ecs->target);
    m_kinds.push_back(ecs->kind);
    m_products.push_back(ecs->product);
    m_massRatios.push_back(ecs->mass_ratio);
    m_thresholds.push_back(ecs->threshold);

    // transpose data
    std::vector<vector_fp> transdata(2, vector_fp(ecs->data.size()));
    for (size_t i = 0; i < ecs->data.size(); i++) {
        for (size_t j = 0; j < 2; j++) {
            transdata[j][i] = ecs->data[i][j];
        }
    }
    m_crossSections.push_back(transdata);

    // shift factor
    if (ecs->kind == "IONIZATION") {
        m_shiftFactor.push_back(2);
    } else {
        m_shiftFactor.push_back(1);
    }

    // scattering-in factor
    if (ecs->kind == "IONIZATION") {
        m_inFactor.push_back(2);
    } else if (ecs->kind == "ATTACHMENT") {
        m_inFactor.push_back(0);
    } else {
        m_inFactor.push_back(1);
    }

    if (ecs->kind == "EFFECTIVE" || ecs->kind == "ELASTIC") {
        for (size_t k = 0; k < m_ncs; k++) {
            if (m_targets[k] == ecs->target)
                if (m_kinds[k] == "ELASTIC" || m_kinds[k] == "EFFECTIVE") {
                    throw CanteraError("Electron::addElectronCrossSection",
                                       "Already contains a data of EFFECTIVE/ELASTIC cross section for '{}'.",
                                       ecs->target);
            }
        }
        m_kElastic.push_back(m_ncs);
    } else {
        m_kInelastic.push_back(m_ncs);
    }

    // add one to number of cross sections
    m_ncs++;

    m_f0_ok = false;

    return true;
}

void Electron::calculateElasticCrossSection()
{
    for (size_t ke : m_kElastic) {
        if (m_kinds[ke] == "EFFECTIVE") {
            // replace effective with elastic
            m_kinds[ke] = "ELASTIC";
            for (size_t k : m_kInelastic) {
                if (m_targets[k] == m_targets[ke]) {
                    vector_fp& x = m_crossSections[k][0];
                    vector_fp& y = m_crossSections[k][1];
                    for (size_t i = 0; i < m_crossSections[ke][0].size(); i++) {
                        m_crossSections[ke][1][i] -= linearInterp(m_crossSections[ke][0][i], x, y);
                    }
                }
            }
            // replace negative values with zero.
            for (size_t i = 0; i < m_crossSections[ke][0].size(); i++) {
                m_crossSections[ke][1][i] = std::max(0.0, m_crossSections[ke][1][i]);
            }
        }
    }
}

void Electron::setupGrid(size_t n, const double* eps)
{
    m_points = n-1;
    m_gridCenter.resize(n-1);
    m_gridEdge.resize(n);
    m_f0.resize(m_points);
    m_gridEdge[n-1] = eps[n-1];
    for (size_t i = 0; i < m_points; i++) {
        m_gridEdge[i] = eps[i];
        m_gridCenter[i] = 0.5 * (eps[i] + eps[i+1]);
    }
    setGridCache();
    m_f0_ok = false;
}

}
