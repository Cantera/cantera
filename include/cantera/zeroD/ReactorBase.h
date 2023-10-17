//! @file ReactorBase.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORBASE_H
#define CT_REACTORBASE_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

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
class ReactorBase
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

    //! @}

    //! Set the state of the Phase object associated with this reactor to the
    //! reactor's current state.
    virtual void restoreState();

    //! Set the state of the reactor to the associated ThermoPhase object.
    //! This method is the inverse of restoreState() and will trigger integrator
    //! reinitialization.
    virtual void syncState();

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
    double density() const {
        if (m_state.empty()) {
            throw CanteraError("ReactorBase::density",
                               "Reactor state empty and/or contents not defined.");
        }
        return m_state[1];
    }

    //! Returns the current temperature (K) of the reactor's contents.
    double temperature() const {
        if (m_state.empty()) {
            throw CanteraError("ReactorBase::temperature",
                               "Reactor state empty and/or contents not defined.");
        }
        return m_state[0];
    }

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
    const double* massFractions() const {
        if (m_state.empty()) {
            throw CanteraError("ReactorBase::massFractions",
                               "Reactor state empty and/or contents not defined.");
        }
        return m_state.data() + 2;
    }

    //! Return the mass fraction of the *k*-th species.
    double massFraction(size_t k) const {
        if (m_state.empty()) {
            throw CanteraError("ReactorBase::massFraction",
                               "Reactor state empty and/or contents not defined.");
        }
        return m_state[k+2];
    }

    //! @}

    //! The ReactorNet that this reactor belongs to.
    ReactorNet& network();

    //! Set the ReactorNet that this reactor belongs to.
    void setNetwork(ReactorNet* net);

    //! Add a sensitivity parameter associated with the reaction number *rxn*
    virtual void addSensitivityReaction(size_t rxn) {
        throw NotImplementedError("ReactorBase::addSensitivityReaction");
    }

    //! Number of sensitivity parameters associated with this reactor.
    virtual size_t nSensParams() const {
        return m_sensParams.size();
    }

    /*! Calculate the derivative of temperature with respect to the temperature in the
     * heat transfer equation based on the reactor specific equation of state.
     * This function should also transform the state of the derivative to that
     * appropriate for the jacobian's state/
     * @warning This function is an experimental part of the %Cantera API and may be changed
     * or removed without notice.
     * @since New in %Cantera 3.0.
     */
    virtual double temperatureDerivative() {
        throw NotImplementedError("Reactor::temperatureDerivative");
    }

    /*! Calculate the derivative of temperature with respect to the temperature in the
     * heat transfer radiation equation based on the reactor specific equation of state.
     * This function should also transform the state of the derivative to that
     * appropriate for the jacobian's state/
     * @warning This function is an experimental part of the %Cantera API and may be changed
     * or removed without notice.
     * @since New in %Cantera 3.0.
     */
    virtual double temperatureRadiationDerivative() {
        throw NotImplementedError("Reactor::temperatureRadiationDerivative");
    }

    /*! Calculate the derivative of T with respect to the ith species in the heat
     * transfer equation based on the reactor specific equation of state.
     * @param index index of the species the derivative is with respect too
     * @warning This function is an experimental part of the %Cantera API and may be changed
     * or removed without notice.
     * @since New in %Cantera 3.0.
     */
    virtual double moleDerivative(size_t index) {
        throw NotImplementedError("Reactor::moleDerivative");
    }

    /*! Calculate the derivative of T with respect to the ith species in the heat
     * transfer radiation equation based on the reactor specific equation of state.
     * @param index index of the species the derivative is with respect too
     * @warning This function is an experimental part of the %Cantera API and may be changed
     * or removed without notice.
     * @since New in %Cantera 3.0.
     */
    virtual double moleRadiationDerivative(size_t index) {
        throw NotImplementedError("Reactor::moleRadiationDerivative");
    }

    //! Return the index associated with energy of the system
    virtual size_t energyIndex() const { return m_eidx; };

    //! Return the offset between species and state variables
    virtual size_t speciesOffset() const { return m_sidx; };

protected:
    explicit ReactorBase(const string& name="(none)");

    //! Number of homogeneous species in the mixture
    size_t m_nsp = 0;

    //! species offset in the state vector
    const size_t m_sidx = 3;

    //! index of state variable associated with energy
    const size_t m_eidx = 1;

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

    //! The ReactorNet that this reactor is part of
    ReactorNet* m_net = nullptr;

    //! Composite thermo/kinetics/transport handler
    shared_ptr<Solution> m_solution;

    // Data associated each sensitivity parameter
    vector<SensitivityParameter> m_sensParams;
};
}

#endif
