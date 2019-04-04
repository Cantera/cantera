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
 * This class calculates the properties of electron in a weakly ionized gas.
 * Only electron-neutral collisions are considered for calculating the
 * electron energy distribution function.
 *
 * Reference:
 * [1] G. J. M. Hagelaar and L. C. Pitchford
 * "Solving the Boltzmann equation to obtain electron transport
 * coefficients and rate coefficients for fluid models."
 * Plasma Sources Science and Technology 14.4 (2005): 722.
 * doi: https://doi.org/10.1088/0963-0252/14/4/011
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
