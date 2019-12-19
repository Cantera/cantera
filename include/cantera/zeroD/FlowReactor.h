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
        return "FlowReactor";
    }

    /*!
     * @deprecated To be changed after Cantera 2.5.
     */
    virtual int type() const {
        warn_deprecated("FlowReactor::type",
                        "To be changed after Cantera 2.5. "
                        "Return string instead of magic number; use "
                        "FlowReactor::typeStr during transition");
        return FlowReactorType;
    }

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);
    virtual void updateState(doublereal* y);

    void setMassFlowRate(doublereal mdot) {
        m_rho0 = m_thermo->density();
        m_speed = mdot/m_rho0;
        m_speed0 = m_speed;
        m_T = m_thermo->temperature();
        m_P0 = m_thermo->pressure() + m_rho0*m_speed*m_speed;
        m_h0 = m_thermo->enthalpy_mass() + 0.5*m_speed*m_speed;
    }

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
