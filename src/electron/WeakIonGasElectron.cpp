// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/WeakIonGasElectron.h"
#include "cantera/electron/ElectronFactory.h"
#include "cantera/base/utilities.h"

namespace Cantera {

WeakIonGasElectron::WeakIonGasElectron()
{
}

void WeakIonGasElectron::setElectronTemperature(double Te)
{
    m_kTe = Boltzmann * Te / ElectronCharge;
    m_f0_ok = false;
}

void WeakIonGasElectron::calculateDistributionFunction()
{
    Electron::calculateDistributionFunction();
    if (m_kTe == Undef) {
        m_kTe = m_kT;
    }
    for (size_t j = 0; j < m_points; j++) {
        m_f0[j] = 2.0 * pow(1.0/Pi, 0.5) * pow(m_kTe, -3./2.) * std::exp(-m_eps[j]/m_kTe);
        m_df0[j] = -2.0 * pow(1.0/Pi, 0.5) * pow(m_kTe, -5./2.) * std::exp(-m_eps[j]/m_kTe);
    }
}

}
