/**
 * @file EdgeKinetics.h
 *
 * @ingroup chemkinetics
 * @ingroup electrochem
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_EDGEKINETICS_H
#define CT_EDGEKINETICS_H

#include "InterfaceKinetics.h"

namespace Cantera
{
/**
 * Heterogeneous reactions at one-dimensional interfaces between
 * multiple adjacent two-dimensional surfaces.
 */
class EdgeKinetics : public InterfaceKinetics
{
public:
    //! Constructor
    EdgeKinetics() : InterfaceKinetics() {
        m_nDim = 1;
    }

    EdgeKinetics(const EdgeKinetics& right) :
        InterfaceKinetics(right) {
        *this=right;
    }

    EdgeKinetics& operator=(const EdgeKinetics& right) {
        if (this != &right) {
            InterfaceKinetics::operator=(right);
        }
        return *this;
    }

    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const {
        EdgeKinetics* iK = new EdgeKinetics(*this);
        iK->assignShallowPointers(tpVector);
        return iK;
    }

    virtual int type() const {
        warn_deprecated("EdgeKinetics::type",
                        "To be removed after Cantera 2.3.");
        return cEdgeKinetics;
    }

    virtual std::string kineticsType() const {
        return "Edge";
    }
};
}

#endif
