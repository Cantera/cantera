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

    /**
    * Constructor
     *
     */
    EdgeKinetics() : InterfaceKinetics() {}

    /// Destructor.
    virtual ~EdgeKinetics() {}

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

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This virtual routine can be used to duplicate %Kinetics objects
     *  inherited from %Kinetics even if the application only has
     *  a pointer to %Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy  constructor.
     *
     * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
     *                  m_thermo vector within this object
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const {
        EdgeKinetics* iK = new EdgeKinetics(*this);
        iK->assignShallowPointers(tpVector);
        return iK;
    }

    /**
     * Identifies the subclass of the Kinetics manager type.
     * These are listed in mix_defs.h.
     */
    virtual int type() const {
        return cEdgeKinetics;
    }

    // defined in InterfaceKinetics.cpp
    virtual void finalize();
};
}

#endif
