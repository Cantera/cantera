//! @file FlowReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWREACTOR_H
#define CT_FLOWREACTOR_H

#include "IdealGasReactor.h"

namespace Cantera
{

//! Adiabatic flow in a constant-area duct with homogeneous and heterogeneous reactions
//! @ingroup reactorGroup
class FlowReactor : public IdealGasReactor
{
public:
    FlowReactor(shared_ptr<Solution> sol, const string& name="(none)");
    FlowReactor(shared_ptr<Solution> sol, bool clone, const string& name="(none)");

    string type() const override {
        return "FlowReactor";
    }

    bool isOde() const override {
        return false;
    }

    bool timeIsIndependent() const override {
        return false;
    }

    //! Not implemented; FlowReactor implements getStateDae() instead.
    void getState(double* y) override {
        throw NotImplementedError("FlowReactor::getState");
    }

    void getStateDae(double* y, double* ydot) override;
    void updateState(double* y) override;

    //! Not implemented; FlowReactor implements evalDae() instead.
    void eval(double t, double* LHS, double* RHS) override {
        throw NotImplementedError("FlowReactor::eval");
    }

    void evalDae(double t, double* y, double* ydot, double* residual) override;

    void getConstraints(double* constraints) override;
    vector<size_t> steadyConstraints() const override {
        throw CanteraError("FlowReactor::steadyConstraints",
            "FlowReactor is not compatible with time-dependent steady-state solver.");
    }

    //! Mass flow rate through the reactor [kg/s]
    double massFlowRate() {
        return m_u * m_rho * m_area;
    }

    //! Set the mass flow rate through the reactor [kg/s]
    void setMassFlowRate(double mdot);

    //! The current gas speed in the reactor [m/s]
    double speed() const {
        return m_u;
    }

    //! The cross-sectional area of the reactor [m²]
    double area() const override {
        return m_area;
    }

    //! Sets the area of the reactor [m²]
    void setArea(double area) override;

    //! Return the index in the solution vector for this reactor of the component named
    //! *nm*. Possible values for *nm* are "density", "speed", "pressure",
    //! "temperature" or the name of a homogeneous phase species.
    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

protected:
    //! Density [kg/m^3]. First component of the state vector.
    double m_rho = NAN;
    //! Axial velocity [m/s]. Second component of the state vector.
    double m_u = -1.0;
    //! offset to the species equations
    const size_t m_offset_Y = 4;
    //! reactor area [m^2]
    double m_area = 1.0;
    //! temporary storage for species partial molar enthalpies
    vector<double> m_hk;
};

}

#endif
