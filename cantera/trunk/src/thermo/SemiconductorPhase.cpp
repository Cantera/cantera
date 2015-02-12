//! @file SemiconductorPhase.cpp
#include "cantera/thermo/SemiconductorPhase.h"

using namespace std;

namespace Cantera
{
static doublereal JoyceDixon(doublereal r)
{
    return log(r) + 1.0/sqrt(8.0)*r - (3.0/16.0 - sqrt(3.0)/9.0)*r*r;
}


SemiconductorPhase::SemiconductorPhase(std::string infile,
                                       std::string id_) {}

void SemiconductorPhase::getChemPotentials(doublereal* mu) const
{
    getActivityConcentrations(DATA_PTR(m_work));
    mu[0] = ec() + GasConstant*temperature()*(JoyceDixon(m_work[0]/nc()));
    mu[1] = ev() + GasConstant*temperature()*(log(m_work[1]/nv()));
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

/**
 * Energy at the top of the conduction band. By default, energies
 * are referenced to this energy, and so this function simply
 * returns zero.
 */
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
