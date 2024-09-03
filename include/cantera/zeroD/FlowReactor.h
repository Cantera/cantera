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
    using IdealGasReactor::IdealGasReactor; // inherit constructors

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
    void initialize(double t0=0.0) override;
    void syncState() override;
    void updateState(double* y) override;

    //! Not implemented; FlowReactor implements evalDae() instead.
    void eval(double t, double* LHS, double* RHS) override {
        throw NotImplementedError("FlowReactor::eval");
    }

    void evalDae(double t, double* y, double* ydot, double* residual) override;

    void getConstraints(double* constraints) override;

    //! Set the mass flow rate through the reactor [kg/s]
    void setMassFlowRate(double mdot);

    //! The current gas speed in the reactor [m/s]
    double speed() const {
        return m_u;
    }

    //! The cross-sectional area of the reactor [m^2]
    double area() const {
        return m_area;
    }

    //! Sets the area of the reactor [m^2]
    void setArea(double area);

    //! The ratio of the reactor's surface area to volume ratio [m^-1]
    //! @note If the surface area to volume ratio is unspecified by the user,
    //!       this will be calculated assuming the reactor is a cylinder.
    double surfaceAreaToVolumeRatio() const;

    //! Set the reactor's surface area to volume ratio [m^-1]
    void setSurfaceAreaToVolumeRatio(double sa_to_vol) {
        m_sa_to_vol = sa_to_vol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double inletSurfaceAtol() const {
        return m_ss_atol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceAtol(double atol) {
        m_ss_atol = atol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double inletSurfaceRtol() const {
        return m_ss_rtol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceRtol(double rtol) {
        m_ss_rtol = rtol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double inletSurfaceMaxSteps() const {
        return m_max_ss_steps;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceMaxSteps(int max_steps) {
        m_max_ss_steps = max_steps;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double inletSurfaceMaxErrorFailures() const {
        return m_max_ss_error_fails;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceMaxErrorFailures(int max_fails) {
        m_max_ss_error_fails = max_fails;
    }

    //! Return the index in the solution vector for this reactor of the component named
    //! *nm*. Possible values for *nm* are "density", "speed", "pressure",
    //! "temperature", the name of a homogeneous phase species, or the name of a surface
    //! species.
    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

    void updateSurfaceState(double* y) override;

protected:
    //! Density [kg/m^3]. First component of the state vector.
    double m_rho = NAN;
    //! Axial velocity [m/s]. Second component of the state vector.
    double m_u = -1.0;
    //! Pressure [Pa]. Third component of the state vector.
    double m_P = NAN;
    //! Temperature [K]. Fourth component of the state vector.
    double m_T = NAN;
    //! offset to the species equations
    const size_t m_offset_Y = 4;
    //! reactor area [m^2]
    double m_area = 1.0;
    //! reactor surface area to volume ratio [m^-1]
    double m_sa_to_vol = -1.0;
    //! temporary storage for surface species production rates
    vector<double> m_sdot_temp;
    //! temporary storage for species partial molar enthalpies
    vector<double> m_hk;
    //! steady-state relative tolerance, used to determine initial surface coverages
    double m_ss_rtol = 1e-7;
    //! steady-state absolute tolerance, used to determine initial surface coverages
    double m_ss_atol = 1e-14;
    //! maximum number of steady-state coverage integrator-steps
    int m_max_ss_steps = 20000;
    //! maximum number of steady-state integrator error test failures
    int m_max_ss_error_fails = 10;
};
}

#endif
