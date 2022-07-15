//! @file ReactorBase.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORBASE_H
#define CT_REACTORBASE_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

//! Namespace for classes implementing zero-dimensional reactor networks.
namespace Cantera
{

//! @defgroup ZeroD Zero-dimensional reactor networks
//!
//! See https://cantera.org/science/reactors.html for a description of the governing
//! equations for specific reactor types and the methods used for solving networks of
// interconnected reactors.

class FlowDevice;
class WallBase;
class ReactorNet;
class ReactorSurface;
class Kinetics;
class ThermoPhase;

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
 * @ingroup ZeroD
 */
class ReactorBase
{
public:
    explicit ReactorBase(const std::string& name = "(none)");
    virtual ~ReactorBase() {}
    ReactorBase(const ReactorBase&) = delete;
    ReactorBase& operator=(const ReactorBase&) = delete;

    //! String indicating the reactor model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual std::string type() const {
        return "ReactorBase";
    }

    //! Return the name of this reactor
    std::string name() const {
        return m_name;
    }

    //! Set the name of this reactor
    void setName(const std::string& name) {
        m_name = name;
    }

    //! @name Methods to set up a simulation
    //! @{

    //! Set the initial reactor volume. By default, the volume is 1.0 m^3.
    void setInitialVolume(doublereal vol) {
        m_vol = vol;
    }

    //! Specify the mixture contained in the reactor. Note that a pointer to
    //! this substance is stored, and as the integration proceeds, the state of
    //! the substance is modified.
    virtual void setThermoMgr(ThermoPhase& thermo);

    //! Specify chemical kinetics governing the reactor.
    virtual void setKineticsMgr(Kinetics& kin) {
        throw NotImplementedError("ReactorBase::setKineticsMgr");
    }

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

    void addSurface(ReactorSurface* surf);

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
    virtual void initialize(doublereal t0 = 0.0) {
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
    doublereal residenceTime();

    //! @name Solution components
    //!
    //! The values returned are those after the last call to ReactorNet::advance
    //! or ReactorNet::step.
    //! @{

    //! Returns the current volume (m^3) of the reactor.
    doublereal volume() const {
        return m_vol;
    }

    //! Returns the current density (kg/m^3) of the reactor's contents.
    doublereal density() const {
        return m_state[1];
    }

    //! Returns the current temperature (K) of the reactor's contents.
    doublereal temperature() const {
        return m_state[0];
    }

    //! Returns the current enthalpy (J/kg) of the reactor's contents.
    doublereal enthalpy_mass() const {
        return m_enthalpy;
    }

    //! Returns the current internal energy (J/kg) of the reactor's contents.
    doublereal intEnergy_mass() const {
        return m_intEnergy;
    }

    //! Returns the current pressure (Pa) of the reactor.
    doublereal pressure() const {
        return m_pressure;
    }

    //! Returns the mass (kg) of the reactor's contents.
    doublereal mass() const {
        return m_vol * density();
    }

    //! Return the vector of species mass fractions.
    const doublereal* massFractions() const {
        return m_state.data() + 2;
    }

    //! Return the mass fraction of the *k*-th species.
    doublereal massFraction(size_t k) const {
        return m_state[k+2];
    }

    //! @}

    //! The ReactorNet that this reactor belongs to.
    ReactorNet& network();

    //! Set the ReactorNet that this reactor belongs to.
    void setNetwork(ReactorNet* net);

protected:
    //! Number of homogeneous species in the mixture
    size_t m_nsp;

    ThermoPhase* m_thermo;
    double m_vol; //!< Current volume of the reactor [m^3]
    double m_enthalpy; //!< Current specific enthalpy of the reactor [J/kg]
    double m_intEnergy; //!< Current internal energy of the reactor [J/kg]
    double m_pressure; //!< Current pressure in the reactor [Pa]
    vector_fp m_state;
    std::vector<FlowDevice*> m_inlet, m_outlet;

    std::vector<WallBase*> m_wall;
    std::vector<ReactorSurface*> m_surfaces;

    //! Vector of length nWalls(), indicating whether this reactor is on the left (0)
    //! or right (1) of each wall.
    vector_int m_lr;
    std::string m_name;

    //! The ReactorNet that this reactor is part of
    ReactorNet* m_net;
};
}

#endif
