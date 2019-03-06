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

class MaxwellBoltzmannElectron: public Electron
{
public:
    MaxwellBoltzmannElectron(); //!< Default constructor.
};

}

#endif
