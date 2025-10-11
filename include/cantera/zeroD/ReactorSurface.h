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

    //! @deprecated To be removed after Cantera 3.2. Replaced by constructor where
    //!    contents and adjacent reactors are specified
    ReactorSurface(shared_ptr<Solution> sol, const string& name="(none)");

    //! @deprecated To be removed after Cantera 3.2. Replaced by constructor where
    //!    contents and adjacent reactors are specified
    ReactorSurface(shared_ptr<Solution> sol, bool clone, const string& name="(none)");

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
        return m_surf;
    }

    //! Accessor for the InterfaceKinetics object
    Kinetics* kinetics() {
        return m_kinetics;
    }

    //! Set the InterfaceKinetics object for this surface
    //! @deprecated To be removed after %Cantera 3.2. Use constructor with
    //!     Solution object instead.
    void setKinetics(Kinetics* kin);

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

    //! Set the reactor that this Surface interacts with
    //! @deprecated To be removed after %Cantera 3.2. Superseded by constructor taking
    //!     a list of adjacent reactors.
    void setReactor(ReactorBase* reactor);

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

    //! Set the coverages and temperature in the surface phase object to the
    //! values for this surface. The temperature is set to match the bulk phase
    //! of the attached Reactor.
    //! @since  Prior to %Cantera 3.2, this operation was performed by syncState()
    void restoreState() override;

    //! Set the coverages for this ReactorSurface based on the attached SurfPhase.
    //! @since  Behavior changed in %Cantera 3.2 for consistency with
    //!     ReactorBase::syncState(). Previously, this method performed the inverse
    //!     operation of setting the ReactorSurface state based on the SurfPhase and
    //!     attached Reactor.
    void syncState() override;

    void addSensitivityReaction(size_t rxn) override;

    //! Set reaction rate multipliers. `params` is the global vector of
    //! sensitivity parameters. This function is called within
    //! ReactorNet::eval() before the reaction rates are evaluated.
    void setSensitivityParameters(const double* params);

    //! Set reaction rate multipliers back to their initial values. This
    //! function is called within ReactorNet::eval() after all rates have been
    //! evaluated.
    void resetSensitivityParameters();

protected:
    void setThermo(ThermoPhase& thermo) override {}

    //! Set the InterfaceKinetics object for this surface.
    //! Method is needed to prevent compiler warnings by disambiguating from the
    //! non-protected variant.
    //! @since New in %Cantera 3.2.
    //! @deprecated To be removed after %Cantera 3.2. Use constructor with
    //!     Solution object instead.
    void setKinetics(Kinetics& kin) override;

    double m_area = 1.0;

    SurfPhase* m_surf = nullptr;
    Kinetics* m_kinetics = nullptr;
    vector<ReactorBase*> m_reactors;
    vector<double> m_cov;
};

}

#endif
