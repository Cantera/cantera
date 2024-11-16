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
 * See the [Science Reference](../reference/reactors/controlreactor.html) for
 * the governing equations of class Reactor.
 *
 * @ingroup reactorGroup
 */
class Reactor : public ReactorBase
{
public:
    Reactor(shared_ptr<Solution> sol, const string& name="(none)");
    using ReactorBase::ReactorBase; // inherit constructors

    string type() const override {
        return "Reactor";
    }

    //! Indicate whether the governing equations for this reactor type are a system of
    //! ODEs or DAEs. In the first case, this class implements the eval() method. In the
    //! second case, this class implements the evalDae() method.
    virtual bool isOde() const {
        return true;
    }

    //! Indicates whether the governing equations for this reactor are functions of time
    //! or a spatial variable. All reactors in a network must have the same value.
    virtual bool timeIsIndependent() const {
        return true;
    }

    /**
     * Insert something into the reactor. The 'something' must belong to a class
     * that is a subclass of both ThermoPhase and Kinetics.
     * @deprecated Unused; to be removed after %Cantera 3.1.
     */
    template<class G>
    void insert(G& contents) {
        warn_deprecated("Reactor::insert", "Unused; to be removed after Cantera 3.1.");
        setThermo(contents);
        setKinetics(contents);
    }

    using ReactorBase::insert;

    void setChemistry(bool cflag=true) override {
        m_chem = cflag;
    }

    //! Returns `true` if changes in the reactor composition due to chemical reactions
    //! are enabled.
    bool chemistryEnabled() const {
        return m_chem;
    }

    void setEnergy(int eflag=1) override {
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
    virtual void getState(double* y);

    //! Get the current state and derivative vector of the reactor for a DAE solver
    /*!
     *  @param[out] y     state vector representing the initial state of the reactor
     *  @param[out] ydot  state vector representing the initial derivatives of the
     *                    reactor
     */
    virtual void getStateDae(double* y, double* ydot) {
        throw NotImplementedError("Reactor::getStateDae(y, ydot)");
    }

    void initialize(double t0=0.0) override;

    //! Evaluate the reactor governing equations. Called by ReactorNet::eval.
    //! @param[in] t time.
    //! @param[out] LHS pointer to start of vector of left-hand side
    //! coefficients for governing equations, length m_nv, default values 1
    //! @param[out] RHS pointer to start of vector of right-hand side
    //! coefficients for governing equations, length m_nv, default values 0
    virtual void eval(double t, double* LHS, double* RHS);

    /**
     * Evaluate the reactor governing equations. Called by ReactorNet::eval.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[in] ydot rate of change of solution vector, length neq()
     * @param[out] residual residuals vector, length neq()
     */
    virtual void evalDae(double t, double* y, double* ydot, double* residual) {
        throw NotImplementedError("Reactor::evalDae");
    }

    //! Given a vector of length neq(), mark which variables should be
    //! considered algebraic constraints
    virtual void getConstraints(double* constraints) {
        throw NotImplementedError("Reactor::getConstraints");
    }

    void syncState() override;

    //! Set the state of the reactor to correspond to the state vector *y*.
    virtual void updateState(double* y);

    //! Number of sensitivity parameters associated with this reactor
    //! (including walls)
    virtual size_t nSensParams() const;

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
    virtual size_t componentIndex(const string& nm) const;

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual string componentName(size_t k);

    //! Set absolute step size limits during advance
    //! @param limits array of step size limits with length neq
    void setAdvanceLimits(const double* limits);

    //! Check whether Reactor object uses advance limits
    //! @returns           True if at least one limit is set, False otherwise
    bool hasAdvanceLimits() const {
        return !m_advancelimits.empty();
    }

    //! Retrieve absolute step size limits during advance
    //! @param[out] limits array of step size limits with length neq
    //! @returns           True if at least one limit is set, False otherwise
    bool getAdvanceLimits(double* limits) const;

    //! Set individual step size limit for component name *nm*
    //! @param nm component name
    //! @param limit value for step size limit
    void setAdvanceLimit(const string& nm, const double limit);

    //! Calculate the Jacobian of a specific Reactor specialization.
    //! @warning Depending on the particular implementation, this may return an
    //! approximate Jacobian intended only for use in forming a preconditioner for
    //! iterative solvers.
    //! @ingroup derivGroup
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    virtual Eigen::SparseMatrix<double> jacobian() {
        throw NotImplementedError("Reactor::jacobian");
    }

    //! Calculate the reactor-specific Jacobian using a finite difference method.
    //!
    //! This method is used only for informational purposes. Jacobian calculations
    //! for the full reactor system are handled internally by CVODES.
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    Eigen::SparseMatrix<double> finiteDifferenceJacobian();

    //! Use this to set the kinetics objects derivative settings
    virtual void setDerivativeSettings(AnyMap& settings);

    //! Set reaction rate multipliers based on the sensitivity variables in
    //! *params*.
    virtual void applySensitivity(double* params);

    //! Reset the reaction rate multipliers
    virtual void resetSensitivity(double* params);

    //! Return a false if preconditioning is not supported or true otherwise.
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    //!
    //! @since New in %Cantera 3.0
    //!
    virtual bool preconditionerSupported() const {return false;};

protected:
    void setKinetics(Kinetics& kin) override;

    //! Return the index in the solution vector for this reactor of the species
    //! named *nm*, in either the homogeneous phase or a surface phase, relative
    //! to the start of the species terms. Used to implement componentIndex for
    //! specific reactor implementations.
    virtual size_t speciesIndex(const string& nm) const;

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

    virtual void evalSurfaces(double* RHS, double* sdot);

    //! Update the state of SurfPhase objects attached to this reactor
    virtual void updateSurfaceState(double* y);

    //! Update the state information needed by connected reactors, flow devices,
    //! and reactor walls. Called from updateState().
    //! @param updatePressure  Indicates whether to update #m_pressure. Should
    //!     `true` for reactors where the pressure is a dependent property,
    //!     calculated from the state, and `false` when the pressure is constant
    //!     or an independent variable.
    virtual void updateConnected(bool updatePressure);

    //! Get initial conditions for SurfPhase objects attached to this reactor
    virtual void getSurfaceInitialConditions(double* y);

    //! Pointer to the homogeneous Kinetics object that handles the reactions
    Kinetics* m_kin = nullptr;

    double m_vdot = 0.0; //!< net rate of volume change from moving walls [m^3/s]

    double m_Qdot = 0.0; //!< net heat transfer into the reactor, through walls [W]

    double m_mass = 0.0; //!< total mass
    vector<double> m_work;

    //! Production rates of gas phase species on surfaces [kmol/s]
    vector<double> m_sdot;

    vector<double> m_wdot; //!< Species net molar production rates
    vector<double> m_uk; //!< Species molar internal energies
    bool m_chem = false;
    bool m_energy = true;
    size_t m_nv = 0;
    size_t m_nv_surf; //!!< Number of variables associated with reactor surfaces

    vector<double> m_advancelimits; //!< Advance step limit

    // Data associated each sensitivity parameter
    vector<SensitivityParameter> m_sensParams;

    //! Vector of triplets representing the jacobian
    vector<Eigen::Triplet<double>> m_jac_trips;
};
}

#endif
