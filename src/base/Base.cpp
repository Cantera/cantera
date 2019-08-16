//! @file Base.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Base.h"

namespace Cantera
{

SolutionBase::SolutionBase() :
    m_thermo(nullptr),
    m_kinetics(nullptr),
    m_transport(nullptr),
    m_type("<SolutionBase_type>"),
    m_name("<unique_name>")
{}

SolutionBase::SolutionBase(const std::string& infile,
                           const std::string& phasename) :
    SolutionBase()
{
    // this *may* be a spot to load all pieces of a phase
    throw NotImplementedError("SolutionBase constructor from file");
}

void SolutionBase::setThermoPhase(shared_ptr<ThermoPhase> thermo) {
    m_thermo = thermo;
}

void SolutionBase::setKinetics(shared_ptr<Kinetics> kinetics) {
    m_kinetics = kinetics;
}

void SolutionBase::setTransport(shared_ptr<Transport> transport) {
    m_transport = transport;
}
}
