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
                   span<shared_ptr<ReactorBase>> reactors,
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
    InterfaceKinetics* kinetics() {
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

    //! Get the number of Reactor and Reservoir objects adjacent to this surface
    //! @since New in %Cantera 4.0.
    size_t nAdjacent() const {
        return m_reactors.size();
    }

    //! Access an adjacent Reactor or Reservoir
    //! @since New in %Cantera 4.0.
    shared_ptr<ReactorBase> adjacent(size_t n) {
        return m_reactors.at(n);
    }

    //! Set the surface coverages. Array `cov` has length equal to the number of
    //! surface species.
    void setCoverages(span<const double> cov);

    //! Set the surface coverages by name
    void setCoverages(const Composition& cov);

    //! Set the surface coverages by name
    void setCoverages(const string& cov);

    //! Get the surface coverages. Array `cov` should have length equal to the
    //! number of surface species.
    void getCoverages(span<double> cov) const;

    void getState(span<double> y) override;
    void initialize(double t0=0.0) override;
    vector<size_t> initializeSteady() override;
    void updateState(span<const double> y) override;
    void eval(double t, span<double> LHS, span<double> RHS) override;
    void evalSteady(double t, span<double> LHS, span<double> RHS) override;

    void addSensitivityReaction(size_t rxn) override;
    // void addSensitivitySpeciesEnthalpy(size_t k) override;
    void applySensitivity(span<const double> params) override;
    void resetSensitivity(span<const double> params) override;

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    void resetBadValues(span<double> y) override;
    // void setDerivativeSettings(AnyMap& settings) override;

protected:
    double m_area = 1.0;

    shared_ptr<SurfPhase> m_surf;
    shared_ptr<InterfaceKinetics> m_kinetics;
    vector<shared_ptr<ReactorBase>> m_reactors;
};

//! A surface where the state variables are the total number of moles of each species.
//!
//! This class provides the approximate Jacobian elements for interactions between
//! itself and the IdealGasMoleReactor and ConstPressureIdealGasMoleReactor classes
//! needed to work with the AdaptivePreconditioner class.
//!
//! @ingroup reactorGroup
//! @since  New in %Cantera 4.0.
class MoleReactorSurface : public ReactorSurface
{
public:
    using ReactorSurface::ReactorSurface;
    string type() const override { return "MoleReactorSurface"; }
    void initialize(double t0=0.0) override;
    void getState(span<double> y) override;
    void updateState(span<const double> y) override;
    void eval(double t, span<double> LHS, span<double> RHS) override;
    void evalSteady(double t, span<double> LHS, span<double> RHS) override;
    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    void resetBadValues(span<double> y) override;
    void getJacobianElements(vector<Eigen::Triplet<double>>& trips) override;

protected:
    //! Temporary storage for coverages
    vector<double> m_cov_tmp;

    //! Temporary storage for d(moles)/d(moles) scaling factor
    vector<double> m_f_species;

    //! Temporary storage for d(energy)/d(moles) scaling factors
    vector<double> m_f_energy;

    //! Mapping from InterfaceKinetics species index to ReactorNet state vector index.
    //! Set to -1 for species not included in the reactor network state vector.
    vector<Eigen::SparseMatrix<double>::StorageIndex> m_kin2net;

    //! Mapping from InterfaceKinetics species index to the corresponding reactor.
    //! Set to `nullptr` for surface species.
    vector<ReactorBase*> m_kin2reactor;

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
                   span<shared_ptr<ReactorBase>> reactors,
                   bool clone,
                   const string& name="(none)");

    string type() const override {
        return "FlowReactorSurface";
    }

    bool timeIsIndependent() const override { return false; }

    void evalDae(double t, span<const double> y, span<const double> ydot,
                 span<double> residual) override;
    void getStateDae(span<double> y, span<double> ydot) override;
    void getConstraints(span<double> constraints) override;

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
