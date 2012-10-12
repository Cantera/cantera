/**
 *  @file Reactor.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#include "ReactorBase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

/**
 * Class Reactor is a general-purpose class for stirred
 * reactors. The reactor may have an arbitrary number of inlets
 * and outlets, each of which may be connected to a "flow device"
 * such as a mass flow controller, a pressure regulator,
 * etc. Additional reactors may be connected to the other end of
 * the flow device, allowing construction of arbitrary reactor
 * networks.
 *
 * The reactor class integrates the same governing equations no
 * matter what type of reactor is simulated. The differences
 * among reactor types are completely specified by the attached
 * flow devices and the time-dependent user-specified boundary
 * conditions.
 *
 * If an instance of class Reactor is used directly, it will
 * simulate an adiabatic, constant volume reactor with gas-phase
 * chemistry but no surface chemistry. Other reactor types may be
 * simulated by deriving a class from Reactor and overloading
 * method getParams.  This method allows specifying the following
 * in terms of the instantaneous reactor state:
 *
 *  - rate of change of the total volume (m^3/s)
 *  - surface heat loss rate (W)
 *  - species surface production rates (kmol/s)
 *
 * class Reactor inherits from both ReactorBase and
 * FuncEval. ReactorBase provides the basic reactor-like methods
 * that FlowDevice instances can access to determine their mass
 * flow rate. Class FuncEval is the class used to define a system
 * of ODE's to be integrated.
 */

class Reactor : public ReactorBase
{

public:

    //!  Default constructor.
    Reactor();

    /**
     * Destructor. Deletes the integrator.
     */
    virtual ~Reactor() {}

    virtual int type() const {
        return ReactorType;
    }

    /**
     * Insert something into the reactor. The 'something' must
     * belong to a class that is a subclass of both ThermoPhase
     * and Kinetics.
     */
    template<class G>
    void insert(G& contents) {
        setThermoMgr(contents);
        setKineticsMgr(contents);
    }

    void setKineticsMgr(Kinetics& kin) {
        m_kin = &kin;
        if (m_kin->nReactions() == 0) {
            disableChemistry();
        }
    }

    void disableChemistry() {
        m_chem = false;
    }
    void enableChemistry() {
        m_chem = true;
    }

    /// Set the energy equation on or off.
    void setEnergy(int eflag = 1) {
        if (eflag > 0) {
            m_energy = true;
        } else {
            m_energy = false;
        }
    }

    /// Returns 'true' if solution of the energy equation is enabled.
    bool energyEnabled() const {
        return m_energy;
    }

    // overloaded methods of class FuncEval
    virtual size_t neq() {
        return m_nv;
    }

    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    /**
     * Set the mixture to a state consistent with solution
     * vector y.
     */
    virtual void updateState(doublereal* y);

    virtual size_t nSensParams();
    virtual void addSensitivityReaction(size_t rxn);

    virtual std::string sensParamID(int p) {
        return m_pname[p];
    }

    virtual size_t componentIndex(const std::string& nm) const;

protected:
    //! Pointer to the homogeneous Kinetics object that handles the reactions
    Kinetics*   m_kin;

    //! Tolerance on the temperature
    doublereal m_vdot, m_Q;
    vector_fp m_work;
    vector_fp m_sdot;            // surface production rates
    bool m_chem;
    bool m_energy;
    size_t m_nv;

    size_t m_nsens;
    std::vector<size_t> m_pnum;
    std::vector<std::string> m_pname;
    std::vector<size_t> m_nsens_wall;
    vector_fp m_mult_save;

private:
};
}

#endif

