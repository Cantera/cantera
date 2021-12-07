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

    //! Get the net heat transfer rate (for example, through walls) into the reactor [W]
    virtual double qdot() const = 0;

    //! Set the net heat transfer rate (for example, through walls) into the reactor [W]
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
        m_initialize = [this](double t0) { R::initialize(t0); };
        m_syncState = [this]() { R::syncState(); };
        m_getState = [this](std::array<size_t, 1> sizes, double* y) { R::getState(y); };
        m_updateState = [this](std::array<size_t, 1> sizes, double* y) { R::updateState(y); };
        m_updateSurfaceState = [this](std::array<size_t, 1> sizes, double* y) { R::updateSurfaceState(y); };
        m_getSurfaceInitialConditions = [this](std::array<size_t, 1> sizes, double* y) {
            R::getSurfaceInitialConditions(y);
        };
        m_updateConnected = [this](bool updatePressure) {
            R::updateConnected(updatePressure);
        };
        m_eval = [this](std::array<size_t, 2> sizes, double t, double* LHS, double* RHS) { R::eval(t, LHS, RHS); };
        m_evalWalls = [this](double t) { R::evalWalls(t); };
        m_evalSurfaces = [this](std::array<size_t, 1> sizes, double t, double* ydot) {
            return R::evalSurfaces(t, ydot);
        };
        m_componentName = [this](size_t k) { return R::componentName(k); };
        m_componentIndex = [this](const std::string& nm) {
            return R::componentIndex(nm);
        };
        m_speciesIndex = [this](const std::string& nm) {
            return R::speciesIndex(nm);
        };
    }

    virtual void setDelegate(const std::string& name,
                             const std::function<void()>& func,
                             const std::string& when) override
    {
        if (name == "syncState") {
            m_syncState = makeDelegate(func, when,
                [this]() { R::syncState(); }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature void(double)", name);
        }
    }

    virtual void setDelegate(const std::string& name,
                             const std::function<void(bool)>& func,
                             const std::string& when) override
    {
        if (name == "updateConnected") {
            m_updateConnected = makeDelegate(
                func, when,
                [this](bool updatePressure) {
                    R::updateConnected(updatePressure);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature void(bool)", name);
        }
    }

    virtual void setDelegate(const std::string& name,
                             const std::function<void(double)>& func,
                             const std::string& when) override
    {
        if (name == "initialize") {
            m_initialize = makeDelegate(func, when,
                [this](double t0) { R::initialize(t0); }
            );
        } else if (name == "evalWalls") {
            m_evalWalls = makeDelegate(func, when,
                [this](double t) {
                    R::evalWalls(t);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature void(double)", name);
        }
    }

    // For functions with the signature void(double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 1>, double*)>& func,
        const std::string& when) override
    {
        if (name == "getState") {
            m_getState = makeDelegate<1>(func,
                when,
                [this](double* y) {
                    R::getState(y);
                }
            );
        } else if (name == "updateState") {
            m_updateState = makeDelegate<1>(func,
                when,
                [this](double* y) {
                    R::updateState(y);
                }
            );
        } else if (name == "updateSurfaceState") {
            m_updateSurfaceState = makeDelegate<1>(func,
                when,
                [this](double* y) {
                    R::updateSurfaceState(y);
                }
            );
        } else if (name == "getSurfaceInitialConditions") {
            m_getSurfaceInitialConditions = makeDelegate<1>(func,
                when,
                [this](double* y) {
                    R::getSurfaceInitialConditions(y);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature"
                " void(array<size_t, 1>, double*)", name);
        }
    }

    // For functions with the signature double(double, double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(double&, std::array<size_t, 1>, double, double*)>& func,
        const std::string& when) override
    {
        if (name == "evalSurfaces") {
            m_evalSurfaces = makeDelegate<1>(
                func,
                when,
                [this](double t, double* ydot) {
                    return R::evalSurfaces(t, ydot);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature"
                " void(double&, array<size_t, 1>, double, double*)", name);
        }
    }

    // For functions with the signature string(size_t)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(std::string&, size_t)>& func,
        const std::string& when) override
    {
        if (name == "componentName") {
            m_componentName = makeDelegate(func, when,
                [this](size_t k) {
                    return R::componentName(k);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature"
                " void(string&, size_t)", name);
        }
    }

    // For functions with the signature size_t(string)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(size_t&, const std::string&)>& func,
        const std::string& when) override
    {
        if (name == "componentIndex") {
            m_componentIndex = makeDelegate(func, when,
                [this](const std::string& nm) {
                    return R::componentIndex(nm);
                }
            );
        } else if (name == "speciesIndex") {
            m_speciesIndex = makeDelegate(
                func, when,
                [this](const std::string& nm) {
                    return R::speciesIndex(nm);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature"
                " void(size_t&, string)", name);
        }
    }

    // For functions with the signature void(double, double*, double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 2>, double, double*, double*)>& func,
        const std::string& when) override
    {
        if (name == "eval") {
            m_eval = makeDelegate<2>(func,
                when,
                [this](double t, double* LHS, double* RHS) {
                    R::eval(t, LHS, RHS);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature"
                "void(array<size_t, 2>, double, double*, double*)", name);
        }
    }

    // Overrides of Reactor methods

    virtual void initialize(double t0)  override{
        m_initialize(t0);
    }

    virtual void syncState()  override{
        m_syncState();
    }

    virtual void getState(double* y)  override{
        std::array<size_t, 1> sizes{R::neq()};
        m_getState(sizes, y);
    }

    virtual void updateState(double* y)  override{
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

    virtual void updateConnected(bool updatePressure) override{
        m_updateConnected(updatePressure);
    }

    virtual void eval(double t, double* LHS, double* RHS) override {
        std::array<size_t, 2> sizes{R::neq(), R::neq()};
        m_eval(sizes, t, LHS, RHS);
    }

    virtual void evalWalls(double t) override {
        m_evalWalls(t);
    }

    virtual double evalSurfaces(double t, double* ydot) override {
        std::array<size_t, 1> sizes{R::m_nv_surf};
        return m_evalSurfaces(sizes, t, ydot);
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
        ReactorSurface* surf = R::m_surfaces.at(n);
        surf->thermo()->setTemperature(R::m_state[0]);
        surf->syncCoverages();
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
    std::function<double(std::array<size_t, 1>, double, double*)> m_evalSurfaces;
    std::function<std::string(size_t)> m_componentName;
    std::function<size_t(const std::string&)> m_componentIndex;
    std::function<size_t(const std::string&)> m_speciesIndex;
};

}
#endif
