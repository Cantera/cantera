/**
 * @file WeakIonGasElectron.h
 * Header file for class WeakIonGasElectron.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_WEAKIONGASELECTRON_H
#define CT_WEAKIONGASELECTRON_H

#include "cantera/electron/Electron.h"

namespace Cantera
{
/**
 * This class calculates the properties of electron in a gas.
 * assuming a Maxwell-Boltzmann electron energy distribution function.
 * @ingroup electron
 */
class WeakIonGasElectron: public Electron
{
public:
    WeakIonGasElectron(); //!< Default constructor.

    void setElectronTemperature(double Te);

    virtual void calculateDistributionFunction();
};

}

#endif
