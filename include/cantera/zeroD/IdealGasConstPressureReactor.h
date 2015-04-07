/**
 *  @file ConstPressureReactor.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_IDEALGASCONSTP_REACTOR_H
#define CT_IDEALGASCONSTP_REACTOR_H

#include "ConstPressureReactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which
 * may be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the other
 * end of the flow device, allowing construction of arbitrary reactor
 * networks.
 */
class IdealGasConstPressureReactor : public ConstPressureReactor
{
public:
    IdealGasConstPressureReactor() {}

    virtual int type() const {
        return IdealGasConstPressureReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "m", "T", the name
    //! of a homogeneous phase species, or the name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;

protected:
    vector_fp m_hk; //!< Species molar enthalpies
};
}

#endif
