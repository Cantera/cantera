//! @file ReactorDelegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORDELEGATOR_H
#define CT_REACTORDELEGATOR_H

#include "Reactor.h"
#include "cantera/base/Delegator.h"

namespace Cantera
{

class ReactorAccessor
{
public:
    //! Set the number of equations represented by this reactor
    virtual void setNEq(size_t n) = 0;

    //! Get the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual double vdot() const = 0;

    //! Set the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual void setVdot(double v) = 0;

    //! Get the net heat transfer rate (for example, through walls) [W]
    virtual double qdot() const = 0;

    //! Set the net heat transfer rate (for example, through walls) [W]
    virtual void setQdot(double q) = 0;
};

template <class R>
class ReactorDelegator : public Delegator, public R, public ReactorAccessor
{
public:
    ReactorDelegator() {
        m_initialize = [this](double t0) { R::initialize(t0); };
        m_syncState = [this]() { R::syncState(); };
        m_getState = [this](double* y) { R::getState(y); };
        m_updateState = [this](double* y) { R::updateState(y); };
        m_updateSurfaceState = [this](double* y) { R::updateSurfaceState(y); };
        m_getSurfaceInitialConditions = [this](double* y) {
            R::getSurfaceInitialConditions(y);
        };
        m_updateConnected = [this](bool updatePressure) {
            R::updateConnected(updatePressure);
        };
        m_evalEqs = [this](double t, double* y, double* ydot, double* params) {
            R::evalEqs(t, y, ydot, params);
        };
        m_evalWalls = [this](double t) { R::evalWalls(t); };
        m_evalSurfaces = [this](double t, double* ydot) {
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
                [this]() {
                    return std::array<size_t, 1>{R::neq()};
                },
                when,
                [this](double* y) {
                    R::getState(y);
                }
            );
        } else if (name == "updateState") {
            m_updateState = makeDelegate<1>(func,
                [this]() {
                    return std::array<size_t, 1>{R::neq()};
                },
                when,
                [this](double* y) {
                    R::updateState(y);
                }
            );
        } else if (name == "updateSurfaceState") {
            m_updateSurfaceState = makeDelegate<1>(func,
                [this]() {
                    return std::array<size_t, 1>{R::m_nv_surf};
                },
                when,
                [this](double* y) {
                    R::updateSurfaceState(y);
                }
            );
        } else if (name == "getSurfaceInitialConditions") {
            m_getSurfaceInitialConditions = makeDelegate<1>(func,
                [this]() {
                    return std::array<size_t, 1>{R::m_nv_surf};
                },
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

    // For functions with the signature void(double, double*, double*, double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
        const std::string& when) override
    {
        if (name == "evalEqs") {
            m_evalEqs = makeDelegate<3>(func,
                [this]() {
                    return std::array<size_t, 3>{R::neq(), R::neq(), R::nSensParams()};
                },
                when,
                [this](double t, double* y, double* ydot, double* params) {
                    R::evalEqs(t, y, ydot, params);
                }
            );
        } else {
            throw NotImplementedError("ReactorDelegator::setDelegate",
                "For function named '{}' with signature",
                "void(array<size_t, 3>, double, double*, double*, double*)", name);
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
                [this]() {
                    return std::array<size_t, 1>{R::m_nv_surf};
                },
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

    // Overrides of Reactor methods

    virtual void initialize(double t0)  override{
        m_initialize(t0);
    }

    virtual void syncState()  override{
        m_syncState();
    }

    virtual void getState(double* y)  override{
        m_getState(y);
    }

    virtual void updateState(double* y)  override{
        m_updateState(y);
    }

    virtual void updateSurfaceState(double* y) override {
        m_updateSurfaceState(y);
    }

    virtual void getSurfaceInitialConditions(double* y) override {
        m_getSurfaceInitialConditions(y);
    }

    virtual void updateConnected(bool updatePressure) override{
        m_updateConnected(updatePressure);
    }

    virtual void evalEqs(double t, double* y, double* ydot, double* params) override {
        m_evalEqs(t, y, ydot, params);
    }

    virtual void evalWalls(double t) override {
        m_evalWalls(t);
    }

    virtual double evalSurfaces(double t, double* ydot) override {
        return m_evalSurfaces(t, ydot);
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
        return R::m_Q;
    }

    virtual void setQdot(double q) override {
        R::m_Q = q;
    }

private:
    std::function<void(double)> m_initialize;
    std::function<void()> m_syncState;
    std::function<void(double*)> m_getState;
    std::function<void(double*)> m_updateState;
    std::function<void(double*)> m_updateSurfaceState;
    std::function<void(double*)> m_getSurfaceInitialConditions;
    std::function<void(bool)> m_updateConnected;
    std::function<void(double, double*, double*, double*)> m_evalEqs;
    std::function<void(double)> m_evalWalls;
    std::function<double(double, double*)> m_evalSurfaces;
    std::function<std::string(size_t)> m_componentName;
    std::function<size_t(const std::string&)> m_componentIndex;
    std::function<size_t(const std::string&)> m_speciesIndex;
};

}
#endif
