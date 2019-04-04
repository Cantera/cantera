/**
 * @file MaxwellBoltzmannElectron.h
 * Header file for class MaxwellBoltzmannElectron.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_MAXWELLBOLTZMANNELECTRON_H
#define CT_MAXWELLBOLTZMANNELECTRON_H

#include "cantera/electron/Electron.h"

namespace Cantera
{
/**
 * This class calculates the properties of electron in a gas.
 * assuming a Maxwell-Boltzmann electron energy distribution function.
 * @ingroup electron
 */
class MaxwellBoltzmannElectron: public Electron
{
public:
    MaxwellBoltzmannElectron(); //!< Default constructor.

    //! set electron temperature
    void setElectronTemperature(double Te);

    //! calculate electron energy distribution function
    virtual void calculateDistributionFunction();
};

}

#endif
