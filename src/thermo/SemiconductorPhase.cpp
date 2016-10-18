//! @file SemiconductorPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SemiconductorPhase.h"

using namespace std;

namespace Cantera
{
static doublereal JoyceDixon(doublereal r)
{
    return log(r) + 1.0/sqrt(8.0)*r - (3.0/16.0 - sqrt(3.0)/9.0)*r*r;
}

SemiconductorPhase::SemiconductorPhase(std::string infile,
                                       std::string id_) {
    warn_deprecated("class SemiconductorPhase",
                    "To be removed after Cantera 2.3.");
}

void SemiconductorPhase::getChemPotentials(doublereal* mu) const
{
    getActivityConcentrations(m_work.data());
    mu[0] = ec() + RT()*(JoyceDixon(m_work[0]/nc()));
    mu[1] = ev() + RT()*(log(m_work[1]/nv()));
}

// units: kmol/m^3
doublereal SemiconductorPhase::nc() const
{
    doublereal fctr = effectiveMass_e() * Boltzmann * temperature()/
                      (2.0*Pi*Planck_bar*Planck_bar);
    return 2.0*pow(fctr, 1.5)/Avogadro;
}

doublereal SemiconductorPhase::nv() const
{
    doublereal fctr = effectiveMass_h() * Boltzmann * temperature()/
                      (2.0*Pi*Planck_bar*Planck_bar);
    return 2.0*pow(fctr, 1.5)/Avogadro;
}

doublereal SemiconductorPhase::ev() const
{
    return 0.0;
}

doublereal SemiconductorPhase::ec() const
{
    return ev() + bandgap();
}

// private
void SemiconductorPhase::initLengths()
{
    m_work.resize(nSpecies());
}
}
