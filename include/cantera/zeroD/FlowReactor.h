//! @file FlowReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWREACTOR_H
#define CT_FLOWREACTOR_H

#include "Reactor.h"

namespace Cantera
{

//! Adiabatic flow in a constant-area duct.
class FlowReactor : public Reactor
{
public:
    FlowReactor();

    virtual std::string typeStr() const {
        warn_deprecated("FlowReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "FlowReactor";
    }

    virtual std::string type() const {
        return "FlowReactor";
    }

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void eval(double t, double* LHS, double* RHS);
    virtual void updateState(doublereal* y);

    void setMassFlowRate(double mdot);

    void setTimeConstant(doublereal tau) {
        m_fctr = 1.0/tau;
    }

    double speed() const {
        return m_speed;
    }
    double distance() const {
        return m_dist;
    }

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "X" (position),
    //! "U", the name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;

protected:
    doublereal m_speed, m_dist, m_T;
    doublereal m_fctr;
    doublereal m_rho0, m_speed0, m_P0, m_h0;
};
}

#endif
