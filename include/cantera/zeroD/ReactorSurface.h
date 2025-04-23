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
    ReactorSurface(shared_ptr<Solution> sol, const string& name="(none)");
    using ReactorBase::ReactorBase; // inherit constructors

    //! String indicating the wall model implemented.
    string type() const override {
        return "ReactorSurface";
    }

    //! Returns the surface area [m^2]
    double area() const;

    //! Set the surface area [m^2]
    void setArea(double a);

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
    ReactorBase* m_reactor = nullptr;
    vector<double> m_cov;
};

}

#endif
