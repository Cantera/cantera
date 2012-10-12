/**
 *  @file ReactorBase.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_REACTORBASE_H
#define CT_REACTORBASE_H

#include "cantera/thermo/ThermoPhase.h"

/// Namespace for classes implementing zero-dimensional reactor networks.
namespace Cantera
{
class FlowDevice;
class Wall;

const int ReservoirType = 1;
const int ReactorType = 2;
const int FlowReactorType = 3;
const int ConstPressureReactorType = 4;

/**
 * Base class for stirred reactors.
 * Allows using any substance model, with arbitrary
 * inflow, outflow, heat loss/gain, surface chemistry, and
 * volume change.
 */
class ReactorBase
{

public:

    explicit ReactorBase(const std::string& name = "(none)");
    virtual ~ReactorBase() {}

    //-----------------------------------------------------

    virtual int type() const {
        return 0;
    }
    std::string name() const {
        return m_name;
    }
    void setName(const std::string& name) {
        m_name = name;
    }

    /** @name Methods to set up a simulation. */
    //@{


    /**
     * Set the initial reactor volume. By default, the volume is
     * 1.0 m^3.
     */
    void setInitialVolume(doublereal vol) {
        m_vol = vol;
        m_vol0 = vol;
    }

    /**
     * Specify the mixture contained in the reactor. Note that
     * a pointer to this substance is stored, and as the integration
     * proceeds, the state of the substance is modified.
     */
    void setThermoMgr(thermo_t& thermo);

    void addInlet(FlowDevice& inlet);
    void addOutlet(FlowDevice& outlet);
    FlowDevice& inlet(size_t n = 0);
    FlowDevice& outlet(size_t n = 0);

    size_t nInlets() {
        return m_inlet.size();
    }
    size_t nOutlets() {
        return m_outlet.size();
    }
    size_t nWalls() {
        return m_wall.size();
    }

    void addWall(Wall& w, int lr);
    Wall& wall(size_t n);

    /**
     * Initialize the reactor. Must be called after specifying the
     *  (and if necessary the inlet mixture) and before
     * calling advance.
     */
    virtual void initialize(doublereal t0 = 0.0) {
        tilt();
    }

    virtual void start() {}

    //@}

    //! Set the state of the Phase object associated with this reactor to the
    //! reactor's current state.
    void restoreState() {
        if (!m_thermo) {
            throw CanteraError("ReactorBase::restoreState", "No phase defined.");
        }
        m_thermo->restoreState(m_state);
    }

    /// return a reference to the contents.
    thermo_t& contents() {
        return *m_thermo;
    }

    const thermo_t& contents() const {
        return *m_thermo;
    }

    doublereal residenceTime();

    /**
     * @name Solution components.
     * The values returned are those after the last call to advance
     * or step.
     */
    //@{

    //! Returns the current volume of the reactor
    /*!
     * @return  Return the volume in m**3
     */
    doublereal volume() const {
        return m_vol;
    }
    doublereal density() const {
        return m_state[1];
    }
    doublereal temperature() const {
        return m_state[0];
    }
    doublereal enthalpy_mass() const {
        return m_enthalpy;
    }
    doublereal intEnergy_mass() const {
        return m_intEnergy;
    }
    doublereal pressure() const {
        return m_pressure;
    }
    doublereal mass() const {
        return m_vol * density();
    }
    const doublereal* massFractions() const {
        return DATA_PTR(m_state) + 2;
    }
    doublereal massFraction(size_t k) const {
        return m_state[k+2];
    }

    //@}

    int error(const std::string& msg) const {
        writelog("Error: "+msg);
        return 1;
    }

protected:

    //! Number of homogeneous species in the mixture
    size_t m_nsp;

    thermo_t*  m_thermo;
    doublereal m_vol, m_vol0;
    bool m_init;
    size_t m_nInlets, m_nOutlets;
    bool m_open;
    doublereal m_enthalpy;
    doublereal m_intEnergy;
    doublereal m_pressure;
    vector_fp m_state;
    std::vector<FlowDevice*> m_inlet, m_outlet;
    std::vector<Wall*> m_wall;
    vector_int m_lr;
    size_t m_nwalls;
    std::string m_name;
    double m_rho0;

private:

    void tilt(const std::string& method="") const {
        throw CanteraError("ReactorBase::"+method,
                           "ReactorBase method called!");
    }
};
}

#endif

