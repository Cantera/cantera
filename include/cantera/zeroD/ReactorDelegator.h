//! @file ReactorDelegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORDELEGATOR_H
#define CT_REACTORDELEGATOR_H

#include "Reactor.h"
#include "cantera/base/Delegator.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"

namespace Cantera
{

//! An abstract base class for providing access to protected capabilities
//! Reactor objects from delegate methods, which would normally only be able to
//! access public Reactor members.
//!
//! Actual implementations of these methods are found in the templated
//! ReactorDelegator class. The purpose of this base class is so these methods
//! can be accessed by casting a Reactor* to a ReactorAccessor* without needing
//! to know the specific kind of Reactor at compilation time.
class ReactorAccessor
{
public:
    //! Set the number of equations represented by this reactor
    virtual void setNEq(size_t n) = 0;

    //! Get the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual double vdot() const = 0;

    //! Set the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual void setVdot(double v) = 0;

    //! Get the net heat transfer rate (for example, through walls) into the
    //! reactor [W]. This value is initialized and calculated as part of
    //! Reactor::evalWalls().
    virtual double qdot() const = 0;

    //! Set the net heat transfer rate (for example, through walls) into the
    //! reactor [W]. For a value set using this method to affect the calculations done
    //! by Reactor::eval, this method should be called in either a "replace" or "after"
    //! delegate for Reactor::evalWalls().
    virtual void setQdot(double q) = 0;

    //! Set the state of the thermo object to correspond to the state of the reactor
    virtual void restoreThermoState() = 0;

    //! Set the state of the thermo object for surface *n* to correspond to the
    //! state of that surface
    virtual void restoreSurfaceState(size_t n) = 0;
};

//! Delegate methods of the Reactor class to external functions
template <class R>
class ReactorDelegator : public Delegator, public R, public ReactorAccessor
{
public:
    ReactorDelegator() {
        install("initialize", m_initialize, [this](double t0) { R::initialize(t0); });
        install("syncState", m_syncState, [this]() { R::syncState(); });
        install("getState", m_getState,
            [this](std::array<size_t, 1> sizes, double* y) { R::getState(y); });
        install("updateState", m_updateState,
            [this](std::array<size_t, 1> sizes, double* y) { R::updateState(y); });
        install("updateSurfaceState", m_updateSurfaceState,
            [this](std::array<size_t, 1> sizes, double* y) { R::updateSurfaceState(y); });
        install("getSurfaceInitialConditions", m_getSurfaceInitialConditions,
            [this](std::array<size_t, 1> sizes, double* y) {
                R::getSurfaceInitialConditions(y);
            }
        );
        install("updateConnected", m_updateConnected,
            [this](bool updatePressure) { R::updateConnected(updatePressure); });
        install("eval", m_eval,
            [this](std::array<size_t, 2> sizes, double t, double* LHS, double* RHS) {
                R::eval(t, LHS, RHS);
            }
        );
        install("evalWalls", m_evalWalls, [this](double t) { R::evalWalls(t); });
        install("evalSurfaces", m_evalSurfaces,
            [this](std::array<size_t, 3> sizes, double* LHS, double* RHS, double* sdot) {
                R::evalSurfaces(LHS, RHS, sdot);
            }
        );
        install("componentName", m_componentName,
            [this](size_t k) { return R::componentName(k); });
        install("componentIndex", m_componentIndex,
            [this](const std::string& nm) { return R::componentIndex(nm); });
        install("speciesIndex", m_speciesIndex,
            [this](const std::string& nm) { return R::speciesIndex(nm); });
    }

    // Overrides of Reactor methods

    virtual void initialize(double t0) override {
        m_initialize(t0);
    }

    virtual void syncState() override {
        m_syncState();
    }

    virtual void getState(double* y) override {
        std::array<size_t, 1> sizes{R::neq()};
        m_getState(sizes, y);
    }

    virtual void updateState(double* y) override {
        std::array<size_t, 1> sizes{R::neq()};
        m_updateState(sizes, y);
    }

    virtual void updateSurfaceState(double* y) override {
        std::array<size_t, 1> sizes{R::m_nv_surf};
        m_updateSurfaceState(sizes, y);
    }

    virtual void getSurfaceInitialConditions(double* y) override {
        std::array<size_t, 1> sizes{R::m_nv_surf};
        m_getSurfaceInitialConditions(sizes, y);
    }

    virtual void updateConnected(bool updatePressure) override {
        m_updateConnected(updatePressure);
    }

    virtual void eval(double t, double* LHS, double* RHS) override {
        std::array<size_t, 2> sizes{R::neq(), R::neq()};
        m_eval(sizes, t, LHS, RHS);
    }

    virtual void evalWalls(double t) override {
        m_evalWalls(t);
    }

    virtual void evalSurfaces(double* LHS, double* RHS, double* sdot) override
    {
        std::array<size_t, 3> sizes{R::m_nv_surf, R::m_nv_surf, R::m_nsp};
        m_evalSurfaces(sizes, LHS, RHS, sdot);
    }

    virtual std::string componentName(size_t k) override {
        return m_componentName(k);
    }

    virtual size_t componentIndex(const std::string& nm) const override {
        return m_componentIndex(nm);
    }

    virtual size_t speciesIndex(const std::string& nm) const override {
        return m_speciesIndex(nm);
    }

    // Public access to protected Reactor variables needed by derived classes

    virtual void setNEq(size_t n) override {
        R::m_nv = n;
    }

    virtual double vdot() const override {
        return R::m_vdot;
    }

    virtual void setVdot(double v) override {
        R::m_vdot = v;
    }

    virtual double qdot() const override {
        return R::m_Qdot;
    }

    virtual void setQdot(double q) override {
        R::m_Qdot = q;
    }

    virtual void restoreThermoState() override {
        R::m_thermo->restoreState(R::m_state);
    }

    virtual void restoreSurfaceState(size_t n) override {
        R::m_surfaces.at(n)->syncState();
    }

private:
    std::function<void(double)> m_initialize;
    std::function<void()> m_syncState;
    std::function<void(std::array<size_t, 1>, double*)> m_getState;
    std::function<void(std::array<size_t, 1>, double*)> m_updateState;
    std::function<void(std::array<size_t, 1>, double*)> m_updateSurfaceState;
    std::function<void(std::array<size_t, 1>, double*)> m_getSurfaceInitialConditions;
    std::function<void(bool)> m_updateConnected;
    std::function<void(std::array<size_t, 2>, double, double*, double*)> m_eval;
    std::function<void(double)> m_evalWalls;
    std::function<void(std::array<size_t, 3>, double*, double*, double*)> m_evalSurfaces;
    std::function<std::string(size_t)> m_componentName;
    std::function<size_t(const std::string&)> m_componentIndex;
    std::function<size_t(const std::string&)> m_speciesIndex;
};

}
#endif
