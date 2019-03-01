//! @file FlowReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWREACTOR_H
#define CT_FLOWREACTOR_H

#include "IdealGasReactor.h"

namespace Cantera
{

//! Adiabatic flow in a constant-area duct.
class FlowReactor : public IdealGasReactor
{
public:
    //! Note: currently assumes a cylinder
    FlowReactor(double area=1.0,
                double sa_to_vol=-1,
                double ss_atol=1e-14,
                double ss_rtol=1e-7,
                int max_ss_steps = 20000,
                int max_ss_error_fails=10);

    virtual int type() const {
        return FlowReactorType;
    }

    virtual void getState(double* y)
    {
        throw CanteraError("FlowReactor::getState", "Not Implemented!");
    }
    virtual void getState(double* y, double* ydot);

    virtual void initialize(doublereal t0 = 0.0);
    /*!
     * Evaluate the reactor governing equations. Called by ReactorNet::eval.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] params sensitivity parameter vector, length ReactorNet::nparams()
     */
    virtual void evalEqs(double t, double* y,
                         double* ydot, double* params)
    {
        throw CanteraError("FlowReactor::evalEqs", "Not Implemented!");
    }

    /*!
     * Evaluate the reactor governing equations. Called by ReactorNet::eval.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[in] ydot rate of change of solution vector, length neq()
     * @param[in] params sensitivity parameter vector, length ReactorNet::nparams()
     * @param[out] residual resisduals vector, length neq()
     */
    virtual void evalEqs(double t, double* y,
                         double* ydot, double* params,
                         double* residual);

    //! Given a vector of length neq(), mark which variables should be
    //! considered algebraic constraints
    virtual void getConstraints(double* constraints) {
        // mark all variables differential equations unless otherwise specified
        std::fill(constraints, constraints + m_nv, 1.0);
        // the species coverages are algebraic constraints
        std::fill(constraints + m_non_spec_eq + m_nsp, constraints + m_nv, 0.0);
    }


    virtual void syncState();
    virtual void updateState(doublereal* y);

    void setMassFlowRate(doublereal mdot) {
        m_rho = m_thermo->density();
        m_u = mdot/(m_rho * m_area);
    }

    //! The current gas speed in the reactor [m/s]
    double speed() const {
        return m_u;
    }

    //! The area of the reactor [m^2]
    double area() const {
        return m_area;
    }

    //! Sets the area of the reactor [m^2]
    void setArea(double area) {
        m_area = area;
        setMassFlowRate(m_rho * m_u * area);
    }

    //! The ratio of the reactor's surface area to volume ratio [m^-1]
    //! @note If the surface area to volume ratio is unspecified by the user,
    //!       this will be calculated assuming the reactor is a cylinder.
    double surfaceAreaToVolumeRatio() const {
        if (m_sa_to_vol > 0)
            return m_sa_to_vol;

        // assuming a cylinder, volume = Pi * r^2 * L, and perimeter = 2 * Pi * r * L
        // where L is the length, and r is the radius of the reactor
        // hence, perimeter / area = 2 * Pi * r * L / (Pi * L * r^2) = 2 / r
        return 2.0 / sqrt(m_area / Pi);
    }

    //! Set the reactor's surface area to volume ratio [m^-1]
    void setSurfaceAreaToVolumeRatio(double sa_to_vol) {
        m_sa_to_vol = sa_to_vol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setSteadyStateAtol(double atol) {
        m_ss_atol = atol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setSteadyStateRtol(double rtol) {
        m_ss_rtol = rtol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setSteadyStateMaxSteps(int max_steps) {
        m_max_ss_steps = max_steps;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setSteadyStateMaxErrorFailures(int max_fails) {
        m_max_ss_error_fails = max_fails;
    }


    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "X" (position),
    //! "U", the name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;

    virtual void updateSurfaceState(double* y);

protected:
    double m_u, m_T, m_P, m_rho;
    //! offset to the species equations
    const size_t m_non_spec_eq = 4;
    //! reactor area [m^2]
    double m_area;
    //! reactor surface area to volume ratio [m^-1]
    double m_sa_to_vol;
    //! temporary storage for surface species production rates
    vector_fp m_sdot_temp;
    //! temporary storage for species partial molar enthalipes
    vector_fp m_hk;

    //! steady-state relative tolerance, used to determine initial surface coverages
    double m_ss_rtol;
    //! steady-state absolute tolerance, used to determine initial surface coverages
    double m_ss_atol;
    //! maximum number of steady-state coverage integrator-steps
    int m_max_ss_steps;
    //! maximum number of steady-state integrator error test failures
    int m_max_ss_error_fails;
};
}

#endif
