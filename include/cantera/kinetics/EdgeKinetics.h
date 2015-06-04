/**
 * @file EdgeKinetics.h
 *
 * @ingroup chemkinetics
 * @ingroup electrochem
 */
// Copyright 2001  California Institute of Technology

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
    EdgeKinetics() : InterfaceKinetics() {}

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
        return cEdgeKinetics;
    }

    virtual void finalize();
};
}

#endif
