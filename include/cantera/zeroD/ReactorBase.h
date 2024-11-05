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
 * Base class for stirred reactors. Allows using any substance model, with
 * arbitrary inflow, outflow, heat loss/gain, surface chemistry, and volume
 * change.
 * @ingroup reactorGroup
 */
class ReactorBase
{
public:
    explicit ReactorBase(const string& name="(none)");
    //! Instantiate a ReactorBase object with Solution contents.
    //! @param sol  Solution object to be set.
    //! @param name  Name of the reactor.
    //! @since New in %Cantera 3.1.
    ReactorBase(shared_ptr<Solution> sol, const string& name="(none)");
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

    //! Set the Solution specifying the ReactorBase content.
    //! @param sol  Solution object to be set.
    //! @since New in %Cantera 3.1.
    void setSolution(shared_ptr<Solution> sol);

    //! @name Methods to set up a simulation
    //! @{

    //! Set the initial reactor volume. By default, the volume is 1.0 m^3.
    void setInitialVolume(double vol) {
        m_vol = vol;
    }

    //! @deprecated To be removed after %Cantera 3.1. Superseded by setSolution.
    void insert(shared_ptr<Solution> sol);

    //! Specify the mixture contained in the reactor. Note that a pointer to
    //! this substance is stored, and as the integration proceeds, the state of
    //! the substance is modified.
    //! @deprecated To be removed after %Cantera 3.1. Superseded by setSolution.
    void setThermoMgr(ThermoPhase& thermo);

    //! @deprecated To be removed after %Cantera 3.1. Superseded by setSolution.
    void setKineticsMgr(Kinetics& kin);

    //! Enable or disable changes in reactor composition due to chemical reactions.
    virtual void setChemistry(bool cflag = true) {
        throw NotImplementedError("ReactorBase::setChemistry");
    }

    //! Set the energy equation on or off.
    virtual void setEnergy(int eflag = 1) {
        throw NotImplementedError("ReactorBase::setEnergy");
    }

    //! Connect an inlet FlowDevice to this reactor
    void addInlet(FlowDevice& inlet);

    //! Connect an outlet FlowDevice to this reactor
    void addOutlet(FlowDevice& outlet);

    //! Return a reference to the *n*-th inlet FlowDevice connected to this
    //! reactor.
    FlowDevice& inlet(size_t n = 0);

    //! Return a reference to the *n*-th outlet FlowDevice connected to this
    //! reactor.
    FlowDevice& outlet(size_t n = 0);

    //! Return the number of inlet FlowDevice objects connected to this reactor.
    size_t nInlets() {
        return m_inlet.size();
    }

    //! Return the number of outlet FlowDevice objects connected to this
    //! reactor.
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
    void addWall(WallBase& w, int lr);

    //! Return a reference to the *n*-th Wall connected to this reactor.
    WallBase& wall(size_t n);

    virtual void addSurface(ReactorSurface* surf);

    //! Return a reference to the *n*-th ReactorSurface connected to this
    //! reactor
    ReactorSurface* surface(size_t n);

    //! Return the number of surfaces in a reactor
    virtual size_t nSurfs() {
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
    void restoreState();

    //! Set the state of the reactor to correspond to the state of the
    //! associated ThermoPhase object. This is the inverse of restoreState().
    //! Calling this will trigger integrator reinitialization.
    virtual void syncState();

    //! return a reference to the contents.
    ThermoPhase& contents() {
        if (!m_thermo) {
            throw CanteraError("ReactorBase::contents",
                               "Reactor contents not defined.");
        }
        return *m_thermo;
    }

    const ThermoPhase& contents() const {
        if (!m_thermo) {
            throw CanteraError("ReactorBase::contents",
                               "Reactor contents not defined.");
        }
        return *m_thermo;
    }

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

    //! Returns the current internal energy (J/kg) of the reactor's contents.
    double intEnergy_mass() const {
        return m_intEnergy;
    }

    //! Returns the current pressure (Pa) of the reactor.
    double pressure() const {
        return m_pressure;
    }

    //! Returns the mass (kg) of the reactor's contents.
    double mass() const {
        return m_vol * density();
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

protected:
    //! Specify the mixture contained in the reactor. Note that a pointer to
    //! this substance is stored, and as the integration proceeds, the state of
    //! the substance is modified.
    //! @since New in %Cantera 3.1.
    virtual void setThermo(ThermoPhase& thermo);

    //! Specify the kinetics manager for the reactor. Called by setSolution().
    //! @since New in %Cantera 3.1.
    virtual void setKinetics(Kinetics& kin) {
        throw NotImplementedError("ReactorBase::setKinetics");
    }

    //! Number of homogeneous species in the mixture
    size_t m_nsp = 0;

    ThermoPhase* m_thermo = nullptr;
    double m_vol = 1.0; //!< Current volume of the reactor [m^3]
    double m_enthalpy = 0.0; //!< Current specific enthalpy of the reactor [J/kg]
    double m_intEnergy = 0.0; //!< Current internal energy of the reactor [J/kg]
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
};
}

#endif
