/**
 * @file EdgeKinetics.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EDGEKINETICS_H
#define CT_EDGEKINETICS_H

#include "InterfaceKinetics.h"

namespace Cantera
{
/**
 * Heterogeneous reactions at one-dimensional interfaces between
 * multiple adjacent two-dimensional surfaces.
 * @ingroup kineticsmgr
 */
class EdgeKinetics : public InterfaceKinetics
{
public:
    //! Constructor
    EdgeKinetics() {
        m_nDim = 1;
    }

    string kineticsType() const override {
        return "edge";
    }
};
}

#endif
