//! @file ReactorBase.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORBASE_H
#define CT_REACTORBASE_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

//! @defgroup zerodGroup Zero-Dimensional Reactor Networks
//!
//! @details See the [Reactor Science](../reference/reactors/index.html) section of the
//! %Cantera website for a description of the governing equations for specific reactor
//! types and the methods used for solving networks of interconnected reactors.

class FlowDevice;
class WallBase;
class ReactorNet;
class ReactorSurface;
class Kinetics;
class ThermoPhase;
class Solution;
class AnyMap;

enum class SensParameterType {
    reaction,
    enthalpy
};

struct SensitivityParameter
{
    size_t local; //!< local parameter index
    size_t global; //!< global parameter index
    double value; //!< nominal value of the parameter
    SensParameterType type; //!< type of sensitivity parameter
};

/**
 * Base class for reactor objects. Allows using any substance model, with arbitrary
 * inflow, outflow, heat loss/gain, surface chemistry, and volume change, whenever
 * defined.
 * @ingroup reactorGroup
 */
class ReactorBase : public std::enable_shared_from_this<ReactorBase>
{
public:
    //! Instantiate a ReactorBase object with Solution contents.
    //! @param sol  Solution object to be set.
    //! @param name  Name of the reactor.
    //! @since New in %Cantera 3.1.
    ReactorBase(shared_ptr<Solution> sol, const string& name="(none)");

    //! Instantiate a ReactorBase object with Solution contents.
    //! @param sol  Solution object representing the contents of this reactor
    //! @param clone  Determines whether to clone `sol` so that the internal state of
    //!     this reactor is independent of the original Solution object and any Solution
    //!     objects used by other reactors in the network.
    //! @param name  Name of the reactor.
    //! @since Added the `clone` argument in %Cantera 3.2. If not specified, the default
    //!     behavior in %Cantera 3.2 is not to clone the Solution object. This will
    //!     change after %Cantera 3.2 to default to `true`.
    ReactorBase(shared_ptr<Solution> sol, bool clone, const string& name="(none)");

    virtual ~ReactorBase();
    ReactorBase(const ReactorBase&) = delete;
    ReactorBase& operator=(const ReactorBase&) = delete;

    //! String indicating the reactor model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual string type() const {
        return "ReactorBase";
    }

    //! Return the name of this reactor
    string name() const {
        return m_name;
    }

    //! Set the name of this reactor
    void setName(const string& name) {
        m_name = name;
    }

    //! Set the default name of a reactor. Returns `false` if it was previously set.
    bool setDefaultName(map<string, int>& counts);

    //! Access the Solution object used to represent the contents of this reactor.
    //! @since New in %Cantera 3.2
    shared_ptr<Solution> phase() { return m_solution; }

    //! Access the Solution object used to represent the contents of this reactor.
    //! @since New in %Cantera 3.2
    shared_ptr<const Solution> phase() const { return m_solution; }

    //! Indicates whether the governing equations for this reactor are functions of time
    //! or a spatial variable. All reactors in a network must have the same value.
    virtual bool timeIsIndependent() const {
        return true;
    }

    //! Number of equations (state variables) for this reactor
    size_t neq() {
        return m_nv;
    }

    //! @name Methods to set up a simulation
    //! @{

    //! Set the initial reactor volume.
    virtual void setInitialVolume(double vol) {
        throw NotImplementedError("ReactorBase::setInitialVolume",
            "Volume is undefined for reactors of type '{}'.", type());
    }

    //! Returns an area associated with a reactor [m²].
    //! Examples: surface area of ReactorSurface or cross section area of FlowReactor.
    virtual double area() const {
        throw NotImplementedError("ReactorBase::area",
            "Area is undefined for reactors of type '{}'.", type());
    }

    //! Set an area associated with a reactor [m²].
    //! Examples: surface area of ReactorSurface or cross section area of FlowReactor.
    virtual void setArea(double a) {
        throw NotImplementedError("ReactorBase::setArea",
            "Area is undefined for reactors of type '{}'.", type());
    }

    //! Returns `true` if changes in the reactor composition due to chemical reactions
    //! are enabled.
    //! @since New in %Cantera 3.2.
    virtual bool chemistryEnabled() const {
        throw NotImplementedError("ReactorBase::chemistryEnabled",
            "Not implemented for reactor type '{}'.", type());
    }

    //! Enable or disable changes in reactor composition due to chemical reactions.
    //! @since New in %Cantera 3.2.
    virtual void setChemistryEnabled(bool cflag = true) {
        throw NotImplementedError("ReactorBase::setChemistryEnabled",
            "Not implemented for reactor type '{}'.", type());
    }

    //! Returns `true` if solution of the energy equation is enabled.
    //! @since New in %Cantera 3.2.
    virtual bool energyEnabled() const {
        throw NotImplementedError("ReactorBase::energyEnabled",
            "Not implemented for reactor type '{}'.", type());
    }

    //! Set the energy equation on or off.
    //! @since New in %Cantera 3.2.
    virtual void setEnergyEnabled(bool eflag = true) {
        throw NotImplementedError("ReactorBase::setEnergyEnabled",
            "Not implemented for reactor type '{}'.", type());
    }

    //! Connect an inlet FlowDevice to this reactor
    virtual void addInlet(FlowDevice& inlet);

    //! Connect an outlet FlowDevice to this reactor
    virtual void addOutlet(FlowDevice& outlet);

    //! Return a reference to the *n*-th inlet FlowDevice connected to this reactor.
    FlowDevice& inlet(size_t n = 0);

    //! Return a reference to the *n*-th outlet FlowDevice connected to this reactor.
    FlowDevice& outlet(size_t n = 0);

    //! Return the number of inlet FlowDevice objects connected to this reactor.
    size_t nInlets() {
        return m_inlet.size();
    }

    //! Return the number of outlet FlowDevice objects connected to this reactor.
    size_t nOutlets() {
        return m_outlet.size();
    }

    //! Return the number of Wall objects connected to this reactor.
    size_t nWalls() {
        return m_wall.size();
    }

    //! Insert a Wall between this reactor and another reactor.
    /*!
     *  `lr` = 0 if this reactor is to the left of the wall and `lr` = 1 if
     *  this reactor is to the right of the wall. This method is called
     *  automatically for both the left and right reactors by WallBase::install.
     */
    virtual void addWall(WallBase& w, int lr);

    //! Return a reference to the *n*-th Wall connected to this reactor.
    WallBase& wall(size_t n);

    //! Add a ReactorSurface object to a Reactor object.
    //! @attention This method should generally not be called directly by users.
    //!     Reactor and ReactorSurface objects should be connected by providing adjacent
    //!     reactors to the newReactorSurface factory function.
    virtual void addSurface(ReactorSurface* surf);

    //! Return a reference to the *n*-th ReactorSurface connected to this reactor.
    ReactorSurface* surface(size_t n);

    //! Return the number of surfaces in a reactor
    virtual size_t nSurfs() const {
        return m_surfaces.size();
    }

    /**
     * Initialize the reactor. Called automatically by ReactorNet::initialize.
     */
    virtual void initialize(double t0 = 0.0) {
        throw NotImplementedError("ReactorBase::initialize");
    }

    //! Get the current state of the reactor.
    /*!
     *  @param[out] y state vector representing the initial state of the reactor
     */
    virtual void getState(double* y) {
        throw NotImplementedError("ReactorBase::getState");
    }

    //! Get the current state and derivative vector of the reactor for a DAE solver
    /*!
     *  @param[out] y     state vector representing the initial state of the reactor
     *  @param[out] ydot  state vector representing the initial derivatives of the
     *                    reactor
     */
    virtual void getStateDae(double* y, double* ydot) {
        throw NotImplementedError("ReactorBase::getStateDae(y, ydot)");
    }

    //! Evaluate the reactor governing equations. Called by ReactorNet::eval.
    //! @param[in] t time.
    //! @param[out] LHS pointer to start of vector of left-hand side
    //! coefficients for governing equations, length m_nv, default values 1
    //! @param[out] RHS pointer to start of vector of right-hand side
    //! coefficients for governing equations, length m_nv, default values 0
    virtual void eval(double t, double* LHS, double* RHS) {
        throw NotImplementedError("ReactorBase::eval");
    }

    /**
     * Evaluate the reactor governing equations. Called by ReactorNet::eval.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[in] ydot rate of change of solution vector, length neq()
     * @param[out] residual residuals vector, length neq()
     */
    virtual void evalDae(double t, double* y, double* ydot, double* residual) {
        throw NotImplementedError("ReactorBase::evalDae");
    }

    //! Given a vector of length neq(), mark which variables should be
    //! considered algebraic constraints
    virtual void getConstraints(double* constraints) {
        throw NotImplementedError("ReactorBase::getConstraints");
    }

    //! Get the indices of equations that are algebraic constraints when solving the
    //! steady-state problem.
    //!
    //! @warning  This method is an experimental part of the %Cantera API and may be
    //!     changed or removed without notice.
    //! @since New in %Cantera 3.2.
    virtual vector<size_t> steadyConstraints() const {
        throw NotImplementedError("ReactorBase::steadyConstraints");
    }

    //! Set the state of the reactor to correspond to the state vector *y*.
    virtual void updateState(double* y) {
        throw NotImplementedError("ReactorBase::updateState");
    }

    //! Add a sensitivity parameter associated with the enthalpy formation of
    //! species *k*.
    virtual void addSensitivitySpeciesEnthalpy(size_t k) {
        throw NotImplementedError("ReactorBase::addSensitivitySpeciesEnthalpy");
    }

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*.
    virtual size_t componentIndex(const string& nm) const {
        throw NotImplementedError("ReactorBase::componentIndex");
    }

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual string componentName(size_t k) {
        throw NotImplementedError("ReactorBase::componentName");
    }

    //! Get the upper bound on the k-th component of the local state vector.
    virtual double upperBound(size_t k) const {
        throw NotImplementedError("ReactorBase::upperBound");
    }

    //! Get the lower bound on the k-th component of the local state vector.
    virtual double lowerBound(size_t k) const {
        throw NotImplementedError("ReactorBase::lowerBound");
    }

    //! Reset physically or mathematically problematic values, such as negative species
    //! concentrations.
    //! @param[inout] y  current state vector, to be updated; length neq()
    virtual void resetBadValues(double* y) {
        throw NotImplementedError("ReactorBase::resetBadValues");
    }

    //! Get Jacobian elements for this reactor within the full reactor network.
    //!
    //! Indices within `trips` are global indices within the full reactor network. The
    //! reactor is responsible for providing all elements of the Jacobian in the rows
    //! corresponding to its state variables, that is, all derivatives of its state
    //! variables with respect to all state variables in the network.
    //!
    //! @warning  This method is an experimental part of the %Cantera API and may be
    //! changed or removed without notice.
    virtual void getJacobianElements(vector<Eigen::Triplet<double>>& trips) {};

    //! Calculate the Jacobian of a specific reactor specialization.
    //! @warning Depending on the particular implementation, this may return an
    //! approximate Jacobian intended only for use in forming a preconditioner for
    //! iterative solvers.
    //! @ingroup derivGroup
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    virtual Eigen::SparseMatrix<double> jacobian();

    //! Use this to set the kinetics objects derivative settings
    virtual void setDerivativeSettings(AnyMap& settings) {
        throw NotImplementedError("ReactorBase::setDerivativeSettings");
    }

    //! Set reaction rate multipliers based on the sensitivity variables in
    //! *params*.
    virtual void applySensitivity(double* params) {
        throw NotImplementedError("ReactorBase::applySensitivity");
    }

    //! Reset the reaction rate multipliers
    virtual void resetSensitivity(double* params) {
        throw NotImplementedError("ReactorBase::resetSensitivity");
    }

    //! Return a false if preconditioning is not supported or true otherwise.
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    //!
    //! @since New in %Cantera 3.0
    //!
    virtual bool preconditionerSupported() const { return false; };

    //! @}

    //! Set the state of the reactor to the associated ThermoPhase object.
    //! This method will trigger integrator reinitialization.
    //! @deprecated To be removed after %Cantera 4.0. Use ReactorNet::reinitialize to
    //!     indicate a change in state that requires integrator reinitialization.
    virtual void syncState();

    //! Update state information needed by connected reactors, flow devices, and walls.
    //!
    //! Called from updateState() for normal reactor types, and from
    //! ReactorNet::updateState for Reservoir.
    //!
    //! @param updatePressure  Indicates whether to update #m_pressure. Should
    //!     `true` for reactors where the pressure is a dependent property,
    //!     calculated from the state, and `false` when the pressure is constant
    //!     or an independent variable.
    virtual void updateConnected(bool updatePressure);

    //! Return the residence time (s) of the contents of this reactor, based
    //! on the outlet mass flow rates and the mass of the reactor contents.
    double residenceTime();

    //! @name Solution components
    //!
    //! The values returned are those after the last call to ReactorNet::advance
    //! or ReactorNet::step.
    //! @{

    //! Returns the current volume (m^3) of the reactor.
    double volume() const {
        return m_vol;
    }

    //! Returns the current density (kg/m^3) of the reactor's contents.
    double density() const;

    //! Returns the current temperature (K) of the reactor's contents.
    double temperature() const;

    //! Returns the current enthalpy (J/kg) of the reactor's contents.
    double enthalpy_mass() const {
        return m_enthalpy;
    }

    //! Returns the current pressure (Pa) of the reactor.
    double pressure() const {
        return m_pressure;
    }

    //! Returns the mass (kg) of the reactor's contents.
    double mass() const {
        return m_mass;
    }

    //! Return the vector of species mass fractions.
    const double* massFractions() const;

    //! Return the mass fraction of the *k*-th species.
    double massFraction(size_t k) const;

    //! @}

    //! The ReactorNet that this reactor belongs to.
    ReactorNet& network();

    //! Set the ReactorNet that this reactor belongs to.
    void setNetwork(ReactorNet* net);

    //! Get the starting offset for this reactor's state variables within the global
    //! state vector of the ReactorNet.
    size_t offset() const { return m_offset; }

    //! Set the starting offset for this reactor's state variables within the global
    //! state vector of the ReactorNet.
    void setOffset(size_t offset) { m_offset = offset; }

    //! Offset of the first species in the local state vector
    size_t speciesOffset() const {
        return m_nv - m_nsp;
    }

    //! Get scaling factors for the Jacobian matrix terms proportional to
    //! @f$ d\dot{n}_k/dC_j @f$.
    //!
    //! Used to determine contribution of surface phases to the Jacobian.
    //!
    //! @param f_species  Scaling factor for derivatives appearing in the species
    //!     equations. Equal to $1/V$.
    //! @param f_energy  Scaling factor for each species term appearing in the energy
    //!     equation.
    virtual void getJacobianScalingFactors(double& f_species, double* f_energy) {
        throw NotImplementedError("ReactorBase::getJacobianScalingFactors");
    }

    //! Add a sensitivity parameter associated with the reaction number *rxn*
    virtual void addSensitivityReaction(size_t rxn) {
        throw NotImplementedError("ReactorBase::addSensitivityReaction");
    }

    //! Number of sensitivity parameters associated with this reactor.
    virtual size_t nSensParams() const {
        return m_sensParams.size();
    }

protected:
    explicit ReactorBase(const string& name="(none)");

    //! Number of homogeneous species in the mixture
    size_t m_nsp = 0;

    ThermoPhase* m_thermo = nullptr;
    double m_vol = 0.0; //!< Current volume of the reactor [m^3]
    double m_mass = 0.0; //!< Current mass of the reactor [kg]
    double m_enthalpy = 0.0; //!< Current specific enthalpy of the reactor [J/kg]
    double m_pressure = 0.0; //!< Current pressure in the reactor [Pa]
    vector<double> m_state;
    vector<FlowDevice*> m_inlet, m_outlet;

    vector<WallBase*> m_wall;
    vector<ReactorSurface*> m_surfaces;

    //! Vector of length nWalls(), indicating whether this reactor is on the left (0)
    //! or right (1) of each wall.
    vector<int> m_lr;
    string m_name;  //!< Reactor name.
    bool m_defaultNameSet = false;  //!< `true` if default name has been previously set.
    size_t m_nv = 0; //!< Number of state variables for this reactor

    ReactorNet* m_net = nullptr; //!< The ReactorNet that this reactor is part of
    size_t m_offset = 0; //!< Offset into global ReactorNet state vector

    //! Composite thermo/kinetics/transport handler
    shared_ptr<Solution> m_solution;

    // Data associated each sensitivity parameter
    vector<SensitivityParameter> m_sensParams;
};
}

#endif
