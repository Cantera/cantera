// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/PlasmaPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/numerics/funcs.h"
#include <iostream>
#include <limits>
#include <set>

namespace Cantera {

PlasmaPhase::PlasmaPhase()
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

void PlasmaPhase::initPlasma(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    if (phaseNode.hasKey("cross-sections")) {
        if (phaseNode["cross-sections"].is<std::vector<std::string>>()) {
            // 'cross-sections' is a list of target species names to be added from the current
            // file's 'cross-sections' section
            addElectronCrossSections(phaseNode["cross-sections"], rootNode["cross-sections"]);
        } else if (phaseNode["cross-sections"].is<std::string>()) {
            // 'cross-sections' is a keyword applicable to the current file's 'cross-sections'
            // section
            addElectronCrossSections(phaseNode["cross-sections"], rootNode["cross-sections"]);
        } else if (phaseNode["cross-sections"].is<std::vector<AnyMap>>()) {
            // Each item in 'cross-sections' is a map with one item, where the key is
            // a section in another YAML file, and the value is a
            // list of target species names to read from that section 
            for (const auto& crossSectionsNode : phaseNode["cross-sections"].asVector<AnyMap>()) {
                const std::string& source = crossSectionsNode.begin()->first;
                const auto& names = crossSectionsNode.begin()->second;
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    AnyMap crossSections = AnyMap::fromYamlFile(fileName);
                    addElectronCrossSections(crossSections[node], names);
                } else {
                    throw InputFileError("newPlasma", crossSectionsNode,
                        "Could not find species section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("newPlasma", phaseNode["cross-sections"],
                "Could not parse cross-sections declaration of type '{}'",
                phaseNode["cross-sections"].type_str());
        }
    } else if (rootNode.hasKey("cross-sections")) {
        // By default, add all cross sections from the 'cross-sections' section
        addElectronCrossSections(AnyValue("all"), rootNode["cross-sections"]);
    }

    m_f0_ok = false;
    calculateElasticCrossSection();
    setGridCache();
    // set up target index
    m_kTargets.resize(m_ncs);
    for (size_t k = 0; k < m_ncs; k++) {
        m_kTargets[k] = speciesIndex(target(k));
    }
    // set up indices of species which has no cross-section data
    for (size_t k = 0; k < nSpecies(); k++) {
        auto it = std::find(m_kTargets.begin(), m_kTargets.end(), k); 
        if (it == m_kTargets.end()) {
            m_kOthers.push_back(k);
        }
    }
}

void PlasmaPhase::setGridCache()
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
            eps1[i] = clip(m_shiftFactor[k] * m_gridEdge[i] + threshold(k),
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

void PlasmaPhase::setTemperature(const double temp)
{
    Phase::setTemperature(temp);
    m_f0_ok = false;
}

void PlasmaPhase::compositionChanged()
{
    Phase::compositionChanged();
    m_f0_ok = false;
    // warn that a specific species needs cross-section data.
    for (size_t k : m_kOthers) {
        if (moleFraction(k) > 0.01) {
            writelog("Cantera::PlasmaPhase::update_C");
            writelog("\n");
            writelog("Warning: The mole fraction of species {} is more than 0.01",
                    speciesName(k));
            writelog(" but it has no data of cross section.");
            writelog("\n");
        }
    }
}

bool PlasmaPhase::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    ecs->validate();
    m_ecss.push_back(ecs);

    // transpose data
    std::vector<vector_fp> transdata(2, vector_fp(ecs->data.size()));
    for (size_t i = 0; i < ecs->data.size(); i++) {
        for (size_t j = 0; j < 2; j++) {
            transdata[j][i] = ecs->data[i][j];
        }
    }
    m_crossSections.push_back(transdata);

    // shift factor
    if (ecs->kind == "ionization") {
        m_shiftFactor.push_back(2);
    } else {
        m_shiftFactor.push_back(1);
    }

    // scattering-in factor
    if (ecs->kind == "ionization") {
        m_inFactor.push_back(2);
    } else if (ecs->kind == "attachment") {
        m_inFactor.push_back(0);
    } else {
        m_inFactor.push_back(1);
    }

    if (ecs->kind == "effective" || ecs->kind == "elastic") {
        for (size_t k = 0; k < m_ncs; k++) {
            if (target(k) == ecs->target)
                if (kind(k) == "elastic" || kind(k) == "effective") {
                    throw CanteraError("PlasmaPhase::addElectronCrossSection",
                                       "Already contains a data of effective/ELASTIC cross section for '{}'.",
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

void PlasmaPhase::addElectronCrossSections(const AnyValue& crossSections, const AnyValue& names)
{
    if (names.is<std::vector<std::string>>()) {
        // 'names' is a list of target species names which should be found in 'cross-sections'
        for (const auto& name : names.asVector<std::string>()) {
            for (const auto& item : crossSections.asVector<AnyMap>()) {
                if (item["target"].asString() == name) {
                    addElectronCrossSection(newElectronCrossSection(item));
                }
            }
        }
    } else if (names.is<std::string>() && names.asString() == "all") {
        // The keyword 'all' means to add all cross-sections from this source
        for (const auto& item : crossSections.asVector<AnyMap>()) {
            addElectronCrossSection(newElectronCrossSection(item));
        }
    } else {
        throw InputFileError("addElectronCrossSections", names,
            "Could not parse cross-sections declaration of type '{}'", names.type_str());
    }
}

void PlasmaPhase::calculateElasticCrossSection()
{
    for (size_t ke : m_kElastic) {
        if (kind(ke) == "effective") {
            // substract inelastic from effective
            for (size_t k : m_kInelastic) {
                if (target(k) == target(ke)) {
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

void PlasmaPhase::setupGrid(size_t n, const double* eps)
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
