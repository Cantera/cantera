//! @file ReactorSurface.h Header file for class ReactorSurface

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_SURFACE_H
#define CT_REACTOR_SURFACE_H

#include "cantera/zeroD/ReactorBase.h"

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
    shared_ptr<Kinetics> m_kinetics;
    vector<ReactorBase*> m_reactors;
    vector<double> m_sdot; //!< species production rates for all phases
};

}

#endif
