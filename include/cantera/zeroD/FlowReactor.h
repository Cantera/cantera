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
    FlowReactor() = default;

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
    double m_speed = 0.0;
    double m_dist = 0.0;
    double m_T = 0.0;
    double m_fctr = 1.0e10;
    double m_rho0 = 0.0;
    double m_speed0 = 0.0;
    double m_P0 = 0.0;
    double m_h0 = 0.0;
};
}

#endif
