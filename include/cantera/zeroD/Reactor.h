//! @file Reactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#include "ReactorBase.h"
#include "cantera/numerics/eigen_sparse.h"


namespace Cantera
{

class Solution;
class AnyMap;

/**
 * Class Reactor is a general-purpose class for stirred reactors. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may be
 * connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of
 * the flow device, allowing construction of arbitrary reactor networks.
 *
 * The reactor class integrates the same governing equations no matter what
 * type of reactor is simulated. The differences among reactor types are
 * completely specified by the attached flow devices and the time-dependent
 * user-specified boundary conditions.
 *
 * If an instance of class Reactor is used directly, it will simulate an
 * adiabatic, constant volume reactor with gas-phase chemistry but no surface
 * chemistry. Other reactor types may be simulated by deriving a class from
 * Reactor.  This method allows specifying the following in terms of the
 * instantaneous reactor state:
 *
 *  - rate of change of the total volume (m^3/s)
 *  - surface heat loss rate (W)
 *  - species surface production rates (kmol/s)
 *
 * @ingroup ZeroD
 */
class Reactor : public ReactorBase
{
public:
    Reactor();

    virtual std::string type() const {
        return "Reactor";
    }

    /**
     * Insert something into the reactor. The 'something' must belong to a class
     * that is a subclass of both ThermoPhase and Kinetics.
     */
    template<class G>
    void insert(G& contents) {
        setThermoMgr(contents);
        setKineticsMgr(contents);
    }

    void insert(shared_ptr<Solution> sol);

    virtual void setKineticsMgr(Kinetics& kin);

    void setChemistry(bool cflag=true) {
        m_chem = cflag;
    }

    //! Returns `true` if changes in the reactor composition due to chemical reactions are enabled.
    bool chemistryEnabled() const {
        return m_chem;
    }

    void setEnergy(int eflag=1) {
        if (eflag > 0) {
            m_energy = true;
        } else {
            m_energy = false;
        }
    }

    //! Returns `true` if solution of the energy equation is enabled.
    bool energyEnabled() const {
        return m_energy;
    }

    //! Number of equations (state variables) for this reactor
    size_t neq() {
        if (!m_nv) {
            initialize();
        }
        return m_nv;
    }

    //! Get the the current state of the reactor.
    /*!
     *  @param[out] y state vector representing the initial state of the reactor
     */
    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    //! Evaluate the reactor governing equations. Called by ReactorNet::eval.
    //! @param[in] t time.
    //! @param[out] LHS pointer to start of vector of left-hand side
    //! coefficients for governing equations, length m_nv, default values 1
    //! @param[out] RHS pointer to start of vector of right-hand side
    //! coefficients for governing equations, length m_nv, default values 0
    virtual void eval(double t, double* LHS, double* RHS);

    virtual void syncState();

    //! Set the state of the reactor to correspond to the state vector *y*.
    virtual void updateState(doublereal* y);

    //! Number of sensitivity parameters associated with this reactor
    //! (including walls)
    virtual size_t nSensParams();

    //! Add a sensitivity parameter associated with the reaction number *rxn*
    //! (in the homogeneous phase).
    virtual void addSensitivityReaction(size_t rxn);

    //! Add a sensitivity parameter associated with the enthalpy formation of
    //! species *k* (in the homogeneous phase)
    virtual void addSensitivitySpeciesEnthalpy(size_t k);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "volume",
    //! "int_energy", the name of a homogeneous phase species, or the name of a
    //! surface species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

    //! Set absolute step size limits during advance
    //! @param limits array of step size limits with length neq
    void setAdvanceLimits(const double* limits);

    //! Check whether Reactor object uses advance limits
    //! @returns           True if at least one limit is set, False otherwise
    bool hasAdvanceLimits() {
        return !m_advancelimits.empty();
    }

    //! Retrieve absolute step size limits during advance
    //! @param[out] limits array of step size limits with length neq
    //! @returns           True if at least one limit is set, False otherwise
    bool getAdvanceLimits(double* limits);

    //! Set individual step size limit for component name *nm*
    //! @param nm component name
    //! @param limit value for step size limit
    void setAdvanceLimit(const std::string& nm, const double limit);

    //! Method to calculate the reactor specific jacobian
    //! @param t current time of the simulation
    //! @param y pointer to state vector
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    virtual Eigen::SparseMatrix<double> jacobian(double t, double* y) {
        throw NotImplementedError("Reactor::jacobian");
    }

    //! Use this to set the kinetics objects derivative settings
    virtual void setDerivativeSettings(AnyMap& settings);

    //! Set reaction rate multipliers based on the sensitivity variables in
    //! *params*.
    virtual void applySensitivity(double* params);

    //! Reset the reaction rate multipliers
    virtual void resetSensitivity(double* params);

protected:
    //! Return the index in the solution vector for this reactor of the species
    //! named *nm*, in either the homogeneous phase or a surface phase, relative
    //! to the start of the species terms. Used to implement componentIndex for
    //! specific reactor implementations.
    virtual size_t speciesIndex(const std::string& nm) const;

    //! Evaluate terms related to Walls. Calculates #m_vdot and #m_Qdot based on
    //! wall movement and heat transfer.
    //! @param t     the current time
    virtual void evalWalls(double t);

    //! Evaluate terms related to surface reactions.
    //! @param[out] LHS   Multiplicative factor on the left hand side of ODE for surface
    //!                   species coverages
    //! @param[out] RHS   Right hand side of ODE for surface species coverages
    //! @param[out] sdot  array of production rates of bulk phase species on surfaces
    //!                   [kmol/s]
    virtual void evalSurfaces(double* LHS, double* RHS, double* sdot);

    //! Update the state of SurfPhase objects attached to this reactor
    virtual void updateSurfaceState(double* y);

    //! Update the state information needed by connected reactors and flow
    //! devices. Called from updateState().
    //! @param updatePressure  Indicates whether to update #m_pressure. Should
    //!     `true` for reactors where the pressure is a dependent property,
    //!     calculated from the state, and `false` when the pressure is constant
    //!     or an independent variable.
    virtual void updateConnected(bool updatePressure);

    //! Get initial conditions for SurfPhase objects attached to this reactor
    virtual void getSurfaceInitialConditions(double* y);

    //! Pointer to the homogeneous Kinetics object that handles the reactions
    Kinetics* m_kin;

    doublereal m_vdot; //!< net rate of volume change from moving walls [m^3/s]

    double m_Qdot; //!< net heat transfer into the reactor, through walls [W]

    doublereal m_mass; //!< total mass
    vector_fp m_work;

    //! Production rates of gas phase species on surfaces [kmol/s]
    vector_fp m_sdot;

    vector_fp m_wdot; //!< Species net molar production rates
    vector_fp m_uk; //!< Species molar internal energies
    bool m_chem;
    bool m_energy;
    size_t m_nv;
    size_t m_nv_surf; //!!< Number of variables associated with reactor surfaces

    vector_fp m_advancelimits; //!< Advance step limit

    // Data associated each sensitivity parameter
    std::vector<SensitivityParameter> m_sensParams;

    //! Vector of triplets representing the jacobian
    std::vector<Eigen::Triplet<double>> m_jac_trips;
};
}

#endif
