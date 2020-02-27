/**
 * @file Inlet1D.h
 *
 * Boundary objects for one-dimensional simulations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BDRY1D_H
#define CT_BDRY1D_H

#pragma message("warning: Inlet1D.h is renamed to Boundary1D.h and will be removed after Cantera 2.5.")

#include "Boundary1D.h"

namespace Cantera
{

/*!
 * Renamed base class for boundaries between one-dimensional spatial domains.
 * @deprecated To be removed after Cantera 2.5.
 */
class Bdry1D : public Boundary1D
{
public:
    Bdry1D() : Boundary1D() {
        warn_deprecated("Bdry1D::Bdry1D()",
                        "To be removed after Cantera 2.5. "
                        "Class renamed to Boundary1D.");
    }
};

}

#endif
