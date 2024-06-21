//! @file ReactorEdge.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorEdge.h"
#include "cantera/zeroD/ReactorNet.h"
//#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

double ReactorEdge::length() const
{
    return m_length;
}

void ReactorEdge::setLength(double a)
{
    m_length = a;
}

void ReactorEdge::setKinetics(Kinetics* kin) {
    m_kinetics = kin;
    if (kin == nullptr) {
        m_thermo = nullptr;
        return;
    }

    m_thermo = dynamic_cast<EdgePhase*>(&kin->thermo(0));
    if (m_thermo == nullptr) {
        throw CanteraError("ReactorEdge::setKinetics",
            "Specified kinetics manager does not represent an edge "
            "kinetics mechanism.");
    }
    m_cov.resize(m_thermo->nSpecies());
    m_thermo->getCoverages(m_cov.data());
}

void ReactorEdge::setReactor(ReactorBase* reactor)
{
    m_reactor = reactor;
}

void ReactorEdge::setCoverages(const double* cov)
{
    copy(cov, cov + m_cov.size(), m_cov.begin());
}

void ReactorEdge::setCoverages(const Composition& cov)
{
    m_thermo->setCoveragesByName(cov);
    m_thermo->getCoverages(m_cov.data());
}

void ReactorEdge::setCoverages(const string& cov)
{
    m_thermo->setCoveragesByName(cov);
    m_thermo->getCoverages(m_cov.data());
}

void ReactorEdge::getCoverages(double* cov) const
{
    copy(m_cov.begin(), m_cov.end(), cov);
}

void ReactorEdge::syncState()
{
    m_thermo->setTemperature(m_reactor->temperature());
    m_thermo->setCoveragesNoNorm(m_cov.data());
}

}
