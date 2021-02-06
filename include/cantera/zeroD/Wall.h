//! @file Wall.h Header file for base class WallBase.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WALL_H
#define CT_WALL_H

#include "cantera/base/ctexceptions.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/ReactorBase.h"

namespace Cantera
{

class Kinetics;
class SurfPhase;
class Func1;

//! Magic numbers
//! @deprecated To be removed after Cantera 2.5.
const int WallType = 1;

/**
 * Base class for 'walls' (walls, pistons, etc.) connecting reactors.
 * @ingroup reactor0
 */
class WallBase
{
public:
    WallBase();

    virtual ~WallBase() {}
    WallBase(const WallBase&) = delete;
    WallBase& operator=(const WallBase&) = delete;

    //! String indicating the wall model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual std::string type() const {
        return "WallBase";
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (i.e. constant volume), but may be overloaded.
     */
    virtual double vdot(double t) {
        return 0.0;
    }

    //! Heat flow rate through the wall (W).
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (i.e. adiabatic wall), but may be overloaded.
     */
    virtual double Q(double t) {
        return 0.0;
    }

    //! Area in (m^2).
    double area() {
        return m_area;
    }

    //! Set the area [m^2].
    virtual void setArea(double a);

    //! Get the area [m^2]
    /*!
     * Redundant function (same as WallBase::area()).
     * @deprecated To be removed after Cantera 2.5.
     */
    double getArea() const {
        warn_deprecated("WallBase::getArea",
                        "To be removed after Cantera 2.5. "
                        "Replace with WallBase::area.");
        return m_area;
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

    double m_area;
};

//! Represents a wall between between two ReactorBase objects.
/*!
 * Walls can move (changing the volume of the adjacent reactors) and allow heat
 * transfer between reactors.
 */
class Wall : public WallBase
{
public:
    Wall();

    //! String indicating the wall model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual std::string type() const {
        return "Wall";
    }

    //! Set the wall velocity to a specified function of time, i.e. \f$ v(t) \f$.
    void setVelocity(Func1* f=0) {
        if (f) {
            m_vf = f;
        }
    }

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
    virtual double vdot(double t);

    //! Specify the heat flux function \f$ q_0(t) \f$.
    void setHeatFlux(Func1* q) {
        m_qf = q;
    }

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
    virtual double Q(double t);

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
    double m_k;

    //! heat transfer coefficient
    double m_rrth;

    //! emissivity
    double m_emiss;

    //! Velocity function
    Func1* m_vf;

    //! Heat flux function
    Func1* m_qf;
};

}

#endif
