// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/Electron.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include <iostream>
#include <limits>
#include <set>

namespace Cantera {

Electron::Electron()
    : m_electronCrossSectionTargets(0)
    , m_electronCrossSectionKinds(0)
    , m_ncs(0)
{
}

Electron::~Electron()
{
}

size_t Electron::nElectronCrossSections() const
{
    return m_ncs;
}

bool Electron::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    if (std::find(m_electronCrossSectionTargets.begin(),
                  m_electronCrossSectionTargets.end(),
                  ecs->target) != m_electronCrossSectionTargets.end()) {
        if (std::find(m_electronCrossSectionKinds.begin(),
                      m_electronCrossSectionKinds.end(),
                      ecs->kind) != m_electronCrossSectionKinds.end()) {
            throw CanteraError("Electron::addElectronCrossSection",
                                "Already contains a data of type '{}' for '{}'.",
                                ecs->kind, ecs->target);
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
    boost::math::barycentric_rational<double> cs(transdata[0].data(),
                                                 transdata[1].data(),
                                                 transdata[0].size());
    m_electronCrossSectionFunctions.push_back(cs);
    m_ncs++;

    return true;
}

}
