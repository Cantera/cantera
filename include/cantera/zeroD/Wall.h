//! @file Wall.h Header file for base class WallBase.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WALL_H
#define CT_WALL_H

#include "cantera/base/ctexceptions.h"
#include "cantera/zeroD/ReactorBase.h"
#include "ConnectorNode.h"

namespace Cantera
{

class Func1;

/**
 * Base class for 'walls' (walls, pistons, etc.) connecting reactors.
 * @ingroup connectorGroup
 */
class WallBase : public ConnectorNode
{
public:
    WallBase(shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1,
             const string& name="(none)");
    using ConnectorNode::ConnectorNode;  // inherit constructors

    string type() const override {
        return "WallBase";
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors at current reactor
    //! network time.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, constant volume), but may be overloaded.
     * @since New in %Cantera 3.0.
     */
    virtual double expansionRate() {
        return 0.0;
    }

    //! Heat flow rate through the wall (W) at current reactor network time.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, an adiabatic wall), but may be overloaded.
     * @since New in %Cantera 3.0.
     */
    virtual double heatRate() {
        return 0.0;
    }

    //! Area in (m^2).
    double area() {
        return m_area;
    }

    //! Set the area [m^2].
    virtual void setArea(double a);

    //! Called just before the start of integration
    virtual void initialize() {}

    //! True if the wall is correctly configured and ready to use.
    virtual bool ready() {
        return (m_left != 0 && m_right != 0);
    }

    //! Return a reference to the Reactor or Reservoir to the left of the wall.
    ReactorBase& left() const {
        return *m_left;
    }

    //! Return a reference to the Reactor or Reservoir to the right of the wall.
    ReactorBase& right() {
        return *m_right;
    }

    //! Set current reactor network time
    /*!
     * @since New in %Cantera 3.0.
     */
    void setSimTime(double time) {
        m_time = time;
    }

protected:
    ReactorBase* m_left = nullptr;
    ReactorBase* m_right = nullptr;

    //! current reactor network time
    double m_time = 0.0;

    double m_area = 1.0;
};

//! Represents a wall between between two ReactorBase objects.
/*!
 * Walls can move (changing the volume of the adjacent reactors) and allow heat
 * transfer between reactors.
 * @ingroup connectorGroup
 */
class Wall : public WallBase
{
public:
    using WallBase::WallBase;  // inherit constructors

    //! String indicating the wall model implemented. Usually
    //! corresponds to the name of the derived class.
    string type() const override {
        return "Wall";
    }

    //! Wall velocity @f$ v(t) @f$ at current reactor network time.
    //! @since New in %Cantera 3.0.
    double velocity() const;

    //! Set the wall velocity to a specified function of time, @f$ v(t) @f$.
    //! @since  Changed in %Cantera 3.2. Previous version used a raw pointer.
    void setVelocity(shared_ptr<Func1> f) {
        if (f) {
            m_vf = f.get();
        }
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * The volume rate of change is given by
     * @f[
     *     \dot V = K A (P_{left} - P_{right}) + F(t)
     * @f]
     * where *K* is the specified expansion rate coefficient, *A* is the wall area,
     * and and *F(t)* is a specified function evaluated at the current network time.
     * Positive values for `expansionRate` correspond to increases in the volume of
     * reactor on left, and decreases in the volume of the reactor on the right.
     * @since New in %Cantera 3.0.
     */
    double expansionRate() override;

    //! Heat flux function @f$ q_0(t) @f$ evaluated at current reactor network time.
    //! @since New in %Cantera 3.0.
    double heatFlux() const;

    //! Specify the heat flux function @f$ q_0(t) @f$.
    //! @since  Changed in %Cantera 3.2. Previous version used a raw pointer.
    void setHeatFlux(shared_ptr<Func1> q) {
        m_qf = q.get();
    }

    //! Heat flow rate through the wall (W).
    /*!
     * The heat flux is given by
     * @f[
     *     Q = h A (T_{left} - T_{right}) + A G(t)
     * @f]
     * where *h* is the heat transfer coefficient, *A* is the wall area, and
     * *G(t)* is a specified function of time evaluated at the current network
     * time. Positive values denote a flux from left to right.
     * @since New in %Cantera 3.0.
     */
    double heatRate() override;

    //! Set the thermal resistance of the wall [K*m^2/W].
    void setThermalResistance(double Rth) {
        m_rrth = 1.0/Rth;
    }

    //! Set the overall heat transfer coefficient [W/m^2/K].
    void setHeatTransferCoeff(double U) {
        m_rrth = U;
    }

    //! Get the overall heat transfer coefficient [W/m^2/K].
    double getHeatTransferCoeff() const {
        return m_rrth;
    }

    //! Set the emissivity.
    void setEmissivity(double epsilon) {
        if (epsilon > 1.0 || epsilon < 0.0) {
            throw CanteraError("WallBase::setEmissivity",
                               "emissivity must be between 0.0 and 1.0");
        }
        m_emiss = epsilon;
    }

    //! Get the emissivity.
    double getEmissivity() const {
        return m_emiss;
    }

    //! Set the expansion rate coefficient.
    void setExpansionRateCoeff(double k) {
        m_k = k;
    }

    //! Get the expansion rate coefficient
    double getExpansionRateCoeff() const {
        return m_k;
    }

protected:

    //! expansion rate coefficient
    double m_k = 0.0;

    //! heat transfer coefficient
    double m_rrth = 0.0;

    //! emissivity
    double m_emiss = 0.0;

    //! Velocity function
    Func1* m_vf = nullptr;

    //! Heat flux function
    Func1* m_qf = nullptr;
};

}

#endif
