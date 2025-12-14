//! @file ReactorSurface.h Header file for class ReactorSurface

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_SURFACE_H
#define CT_REACTOR_SURFACE_H

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/kinetics/InterfaceKinetics.h"

namespace Cantera
{

class SurfPhase;

//! A surface where reactions can occur that is in contact with the bulk fluid of a
//! Reactor.
//! @ingroup reactorGroup
class ReactorSurface : public ReactorBase
{
public:
    //! Create a new ReactorSurface
    //! @param soln  Thermodynamic and kinetic model representing species and reactions
    //!     on the surface
    //! @param clone  Determines whether to clone `soln` so that the internal state of
    //!     this reactor surface is independent of the original Solution (Interface)
    //!     object and any Solution objects used by other reactors in the network.
    //! @param reactors  One or more reactors whose phases participate in reactions
    //!     occurring on the surface. For the purpose of rate evaluation, the
    //!     temperature of the surface is set equal to the temperature of the first
    //!     reactor specified.
    //! @param name  Name used to identify the surface
    //! @since  Constructor signature including `reactors` and `clone` arguments
    //!     introduced in %Cantera 3.2.
    ReactorSurface(shared_ptr<Solution> soln,
                   const vector<shared_ptr<ReactorBase>>& reactors,
                   bool clone,
                   const string& name="(none)");

    //! String indicating the wall model implemented.
    string type() const override {
        return "ReactorSurface";
    }

    //! Returns the surface area [m²]
    double area() const override;

    //! Set the surface area [m²]
    void setArea(double a) override;

    //! Accessor for the SurfPhase object
    SurfPhase* thermo() {
        return m_surf.get();
    }

    //! Accessor for the InterfaceKinetics object
    Kinetics* kinetics() {
        return m_kinetics.get();
    }

    void addInlet(FlowDevice& inlet) override {
        throw NotImplementedError("ReactorSurface::addInlet",
            "Inlets are undefined for reactors of type '{}'.", type());
    }

    void addOutlet(FlowDevice& outlet) override {
        throw NotImplementedError("ReactorSurface::addOutlet",
            "Outlets are undefined for reactors of type '{}'.", type());
    }

    void addWall(WallBase& w, int lr) override {
        throw NotImplementedError("ReactorSurface::addWall");
    }

    void addSurface(ReactorSurface* surf) override {
        throw NotImplementedError("ReactorSurface::addSurface");
    }

    //! Set the surface coverages. Array `cov` has length equal to the number of
    //! surface species.
    void setCoverages(const double* cov);

    //! Set the surface coverages by name
    void setCoverages(const Composition& cov);

    //! Set the surface coverages by name
    void setCoverages(const string& cov);

    //! Get the surface coverages. Array `cov` should have length equal to the
    //! number of surface species.
    void getCoverages(double* cov) const;

    const vector<double>& netProductionRates() const {
        return m_sdot;
    }

    void getState(double* y) override;
    void initialize(double t0=0.0) override;
    void updateState(double* y) override;
    void eval(double t, double* LHS, double* RHS) override;

    void addSensitivityReaction(size_t rxn) override;
    // void addSensitivitySpeciesEnthalpy(size_t k) override;
    void applySensitivity(double* params) override;
    void resetSensitivity(double* params) override;

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

    // vector<size_t> steadyConstraints() const override;
    // double upperBound(size_t k) const override;
    // double lowerBound(size_t k) const override;
    // void resetBadValues(double* y) override;
    // void setDerivativeSettings(AnyMap& settings) override;

protected:
    double m_area = 1.0;

    shared_ptr<SurfPhase> m_surf;
    shared_ptr<InterfaceKinetics> m_kinetics;
    vector<ReactorBase*> m_reactors;
    vector<double> m_sdot; //!< species production rates for all phases
};

class MoleReactorSurface : public ReactorSurface
{
public:
    using ReactorSurface::ReactorSurface;
    void getState(double* y) override;
    void updateState(double* y) override;
    void eval(double t, double* LHS, double* RHS) override;

protected:
    vector<double> m_cov_tmp;
};

//! A surface in contact with a FlowReactor.
//!
//! May represent the reactor wall or a catalytic surface within the reactor.
//! @ingroup reactorGroup
class FlowReactorSurface : public ReactorSurface
{
public:
    //! @copydoc ReactorSurface::ReactorSurface
    FlowReactorSurface(shared_ptr<Solution> soln,
                   const vector<shared_ptr<ReactorBase>>& reactors,
                   bool clone,
                   const string& name="(none)");

    string type() const override {
        return "FlowReactorSurface";
    }

    void evalDae(double t, double* y, double* ydot, double* residual) override;
    void getStateDae(double* y, double* ydot) override;
    void getConstraints(double* constraints) override;

    //! Surface area per unit length [m]
    //! @note If unspecified by the user, this will be calculated assuming the surface
    //!     is the wall of a cylindrical reactor.
    double area() const override;

    //! Set the reactor surface area per unit length [m].
    //!
    //! If the surface is the wall of the reactor, then this is the perimeter of the
    //! reactor. If the surface represents something like a catalyst monolith, this is
    //! the inverse of the surface area to volume ratio.
    void setArea(double A) override { m_area = A; }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double initialAtol() const {
        return m_ss_atol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInitialAtol(double atol) {
        m_ss_atol = atol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double initialRtol() const {
        return m_ss_rtol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInitialRtol(double rtol) {
        m_ss_rtol = rtol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double initialMaxSteps() const {
        return m_max_ss_steps;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInitialMaxSteps(int max_steps) {
        m_max_ss_steps = max_steps;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    double initialMaxErrorFailures() const {
        return m_max_ss_error_fails;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInitialMaxErrorFailures(int max_fails) {
        m_max_ss_error_fails = max_fails;
    }

protected:
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
