//! @file Wall.h Header file for class Wall.

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_WALL_H
#define CT_WALL_H

#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/Func1.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/zeroD/ReactorSurface.h"

namespace Cantera
{

class Kinetics;
class SurfPhase;

//! Represents a wall between between two ReactorBase objects.
/*!
 * Walls can move (changing the volume of the adjacent reactors) and allow heat
 * transfer between reactors.
 */
class Wall
{
public:
    Wall();

    virtual ~Wall() {}
    Wall(const Wall&) = delete;
    Wall& operator=(const Wall&) = delete;

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * The volume rate of change is given by
     * \f[
     *     \dot V = K A (P_{left} - P_{right}) + F(t)
     * \f]
     * where *K* is the specified expansion rate coefficient, *A* is the wall
     * area, and *F(t)* is a specified function of time. Positive values for
     * `vdot` correspond to increases in the volume of reactor on left, and
     * decreases in the volume of the reactor on the right.
     */
    virtual doublereal vdot(doublereal t);

    //! Heat flow rate through the wall (W).
    /*!
     * The heat flux is given by
     * \f[
     *     Q = h A (T_{left} - T_{right}) + A G(t)
     * \f]
     * where *h* is the heat transfer coefficient, *A* is the wall area, and
     * *G(t)* is a specified function of time. Positive values denote a flux
     * from left to right.
     */
    virtual doublereal Q(doublereal t);

    //! Area in m^2.
    doublereal area() {
        return m_area;
    }

    //! Set the area [m^2].
    void setArea(doublereal a) {
        m_area = a;
        m_surf[0].setArea(a);
        m_surf[1].setArea(a);
    }

    //! Get the area [m^2]
    double getArea() const {
        return m_area;
    }

    void setThermalResistance(doublereal Rth) {
        m_rrth = 1.0/Rth;
    }

    //! Set the overall heat transfer coefficient [W/m^2/K].
    void setHeatTransferCoeff(doublereal U) {
        m_rrth = U;
    }

    //! Get the overall heat transfer coefficient [W/m^2/K].
    double getHeatTransferCoeff() const {
        return m_rrth;
    }

    //! Set the emissivity.
    void setEmissivity(doublereal epsilon) {
        if (epsilon > 1.0 || epsilon < 0.0) {
            throw CanteraError("Wall::setEmissivity",
                               "emissivity must be between 0.0 and 1.0");
        }
        m_emiss = epsilon;
    }

    double getEmissivity() const {
        return m_emiss;
    }

    //! Set the wall velocity to a specified function of time
    void setVelocity(Func1* f=0) {
        if (f) {
            m_vf = f;
        }
    }

    //! Set the expansion rate coefficient.
    void setExpansionRateCoeff(doublereal k) {
        m_k = k;
    }

    //! Get the expansion rate coefficient
    double getExpansionRateCoeff() const {
        return m_k;
    }

    //! Specify the heat flux function \f$ q_0(t) \f$.
    void setHeatFlux(Func1* q) {
        m_qf = q;
    }

    //! Install the wall between two reactors or reservoirs
    bool install(ReactorBase& leftReactor, ReactorBase& rightReactor);

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
    const ReactorBase& right() {
        return *m_right;
    }

protected:
    ReactorBase* m_left;
    ReactorBase* m_right;

    std::vector<ReactorSurface> m_surf;

    doublereal m_area, m_k, m_rrth;
    doublereal m_emiss;
    Func1* m_vf;
    Func1* m_qf;
};

}

#endif
