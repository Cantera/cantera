/**
 *  @file Wall.h
 *  Header file for class Wall.
 */
// Copyright 2001-2004  California Institute of Technology

#ifndef CT_WALL_H
#define CT_WALL_H

#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

// forward references
class ReactorBase;
class Kinetics;
class SurfPhase;

//! Represents a wall between between two ReactorBase objects.
/*!
 * Walls can move (changing the volume of the adjacent reactors), allow heat
 * transfer between reactors, and provide a location for surface reactions to
 * take place.
 */
class Wall
{
public:
    Wall();

    virtual ~Wall() {}

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*! The volume rate of change is given by
     * \f[ \dot V = K A (P_{left} - P_{right}) + F(t) \f]
     * where *K* is the specified expansion rate coefficient, *A* is the wall
     * area, and *F(t)* is a specified function of time. Positive values for
     * `vdot` correspond to increases in the volume of reactor on left, and
     * decreases in the volume of the reactor on the right.
     */
    virtual doublereal vdot(doublereal t);

    //! Heat flow rate through the wall (W).
    /*!
     * The heat flux is given by
     * \f[ Q = h A (T_{left} - T_{right}) + A G(t) \f]
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
        if (epsilon > 1.0 || epsilon < 0.0)
            throw Cantera::CanteraError("Wall::setEmissivity",
                                        "emissivity must be between 0.0 and 1.0");
        m_emiss = epsilon;
    }

    double getEmissivity() const {
        return m_emiss;
    }

    //! Set the wall velocity to a specified function of time
    void setVelocity(Cantera::Func1* f=0) {
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
    void setHeatFlux(Cantera::Func1* q) {
        m_qf = q;
    }

    //! Install the wall between two reactors or reservoirs
    bool install(ReactorBase& leftReactor, ReactorBase& rightReactor);

    //! Called just before the start of integration
    virtual void initialize();

    //! True if the wall is correctly configured and ready to use.
    virtual bool ready() {
        return (m_left != 0 && m_right != 0);
    }

    //! Return a reference to the Reactor or Reservoir to the left
    //! of the wall.
    ReactorBase& left() const {
        return *m_left;
    }

    //! Return a reference to the Reactor or Reservoir to the
    //! right of the wall.
    const ReactorBase& right() {
        return *m_right;
    }

    //! Specify the heterogeneous reaction mechanisms for each side of the
    //! wall. Passing a null pointer indicates that there is no reaction
    //! mechanism for the corresponding wall surface.
    void setKinetics(Cantera::Kinetics* leftMechanism,
                     Cantera::Kinetics* rightMechanism);

    //! Return a pointer to the surface phase object for the left
    //! (`leftright=0`) or right (`leftright=1`) wall surface.
    Cantera::SurfPhase* surface(int leftright) {
        return m_surf[leftright];
    }

    //! Return a pointer to the surface kinetics object for the left
    //! (`leftright=0`) or right (`leftright=1`) wall surface.
    Cantera::Kinetics* kinetics(int leftright) {
        return m_chem[leftright];
    }

    //! Set the surface coverages on the left (`leftright = 0`) or right
    //! (`leftright = 1`) surface to the values in array `cov`.
    void setCoverages(int leftright, const doublereal* cov);

    //! Set the surface coverages on the left (`leftright = 0`) or right
    //! (`leftright = 1`) surface to the values in array `cov`.
    void setCoverages(int leftright, const compositionMap& cov);

    //! Set the surface coverages on the left (`leftright = 0`) or right
    //! (`leftright = 1`) surface to the values in array `cov`.
    void setCoverages(int leftright, const std::string& cov);

    //! Write the coverages of the left or right surface into array `cov`.
    void getCoverages(int leftright, doublereal* cov);

    //! Set the coverages in the surface phase object to the
    //! values for this wall surface.
    void syncCoverages(int leftright);

    //! Number of sensitivity parameters associated with reactions on the left
    //! (`lr = 0`) or right (`lr = 1`) side of the wall.
    size_t nSensParams(int lr) const {
        if (lr == 0) {
            return m_pleft.size();
        } else {
            return m_pright.size();
        }
    }
    void addSensitivityReaction(int leftright, size_t rxn);
    void setSensitivityParameters(int lr, double* params);
    void resetSensitivityParameters(int lr);

protected:
    ReactorBase* m_left;
    ReactorBase* m_right;
    Cantera::Kinetics* m_chem[2];
    Cantera::SurfPhase* m_surf[2];
    size_t m_nsp[2];
    doublereal m_area, m_k, m_rrth;
    doublereal m_emiss;
    Cantera::Func1* m_vf;
    Cantera::Func1* m_qf;
    Cantera::vector_fp m_leftcov, m_rightcov;

    std::vector<size_t> m_pleft, m_pright;
    Cantera::vector_fp m_leftmult_save, m_rightmult_save;
};

}

#endif
