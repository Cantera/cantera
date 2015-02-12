/**
 *  @file ReactorNet.h
 */
// Copyright 2004  California Institute of Technology

#ifndef CT_REACTORNET_H
#define CT_REACTORNET_H

#include "Reactor.h"
#include "cantera/numerics/FuncEval.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/base/Array.h"

namespace Cantera
{

//! A class representing a network of connected reactors.
/*!
 *  This class is used to integrate the time-dependent governing equations for
 *  a network of reactors (Reactor, ConstPressureReactor) connected by various
 *  means, e.g. Wall, MassFlowController, Valve, PressureController.
 */
class ReactorNet : public Cantera::FuncEval
{
public:
    ReactorNet();
    virtual ~ReactorNet();

    /** @name Methods to set up a simulation. */
    //@{

    /**
     * Set initial time. Default = 0.0 s. Restarts integration
     * from this time using the current mixture state as the
     * initial condition.
     */
    void setInitialTime(doublereal time) {
        m_time = time;
        m_integrator_init = false;
    }

    //! Set the maximum time step.
    void setMaxTimeStep(double maxstep) {
        m_maxstep = maxstep;
        m_init = false;
    }

    //! Set the maximum number of error test failures permitted by the CVODES
    //! integrator in a single time step.
    void setMaxErrTestFails(int nmax) {
        m_maxErrTestFails = nmax;
        m_init = false;
    }

    //! Set the relative and absolute tolerances for the integrator.
    void setTolerances(doublereal rtol, doublereal atol) {
        if (rtol >= 0.0) {
            m_rtol = rtol;
        }
        if (atol >= 0.0) {
            m_atols = atol;
        }
        m_init = false;
    }

    //! Set the relative and absolute tolerances for integrating the
    //! sensitivity equations.
    void setSensitivityTolerances(doublereal rtol, doublereal atol) {
        if (rtol >= 0.0) {
            m_rtolsens = rtol;
        }
        if (atol >= 0.0) {
            m_atolsens = atol;
        }
        m_init = false;
    }

    //! Current value of the simulation time.
    doublereal time() {
        return m_time;
    }

    //! Relative tolerance.
    doublereal rtol() {
        return m_rtol;
    }

    //! Absolute integration tolerance
    doublereal atol() {
        return m_atols;
    }

    //! Relative sensitivity tolerance
    doublereal rtolSensitivity() const {
        return m_rtolsens;
    }

    //! Absolute sensitivity tolerance
    doublereal atolSensitivity() const {
        return m_atolsens;
    }

    /**
     * Advance the state of all reactors in time. Take as many internal
     * timesteps as necessary to reach *time*.
     * @param time Time to advance to (s).
     */
    void advance(doublereal time);

    //! Advance the state of all reactors in time. Take a single timestep
    //! toward *time*.
    double step(doublereal time);

    //@}

    //! Add the reactor *r* to this reactor network.
    void addReactor(Reactor& r);

    //! Add the reactor *r* to this reactor network.
    /**
     *  @deprecated To be removed after Cantera 2.2. Use addReactor(Reactor&)
     *  instead.
     */
    void addReactor(Reactor* r, bool iown = false);

    //! Return a reference to the *n*-th reactor in this network. The reactor
    //! indices are determined by the order in which the reactors were added
    //! to the reactor network.
    Reactor& reactor(int n) {
        return *m_reactors[n];
    }

    //! Returns `true` if verbose logging output is enabled.
    bool verbose() const {
        return m_verbose;
    }

    //! Enable or disable verbose logging while setting up and integrating the
    //! reactor network.
    void setVerbose(bool v = true) {
        m_verbose = v;
    }

    //! Return a reference to the integrator.
    Integrator& integrator() {
        return *m_integ;
    }

    //! Update the state of all the reactors in the network to correspond to
    //! the values in the solution vector *y*.
    void updateState(doublereal* y);

    //! Return the sensitivity of the *k*-th solution component with respect to
    //! the *p*-th sensitivity parameter.
    /*!
     *  The normalized sensitivity coefficient \f$ S_{ki} \f$ of solution
     *  variable \f$ y_k \f$ with respect to sensitivity parameter \f$ p_i \f$
     *  is defined as:
     *
     *  \f[ S_{ki} = \frac{p_i}{y_k} \frac{\partial y_k}{\partial p_i} \f]
     *
     *  For reaction sensitivities, the parameter is a multiplier on the forward
     *  rate constant (and implicitly on the reverse rate constant for
     *  reversible reactions).
     */
    double sensitivity(size_t k, size_t p) {
        if (!m_init) {
            initialize();
        }
        return m_integ->sensitivity(k, m_sensIndex[p])/m_integ->solution(k);
    }

    //! Return the sensitivity of the component named *component* with respect to
    //! the *p*-th sensitivity parameter.
    //! @copydetails ReactorNet::sensitivity(size_t, size_t)
    double sensitivity(const std::string& component, size_t p, int reactor=0) {
        size_t k = globalComponentIndex(component, reactor);
        return sensitivity(k, p);
    }

    //! Evaluate the Jacobian matrix for the reactor network.
    /*!
     *  @param[in] t Time at which to evaluate the Jacobian
     *  @param[in] y Global state vector at time *t*
     *  @param[out] ydot Time derivative of the state vector evaluated at *t*.
     *  @param[in] p sensitivity parameter vector (unused?)
     *  @param[out] j Jacobian matrix, size neq() by neq().
     */
    void evalJacobian(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p, Array2D* j);

    // overloaded methods of class FuncEval
    virtual size_t neq() {
        return m_nv;
    }
    virtual void eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p);
    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);
    virtual size_t nparams() {
        return m_ntotpar;
    }

    //! Return the index corresponding to the component named *component* in the
    //! reactor with index *reactor* in the global state vector for the
    //! reactor network.
    size_t globalComponentIndex(const std::string& component, size_t reactor=0);

    //! Used by Reactor and Wall objects to register the addition of
    //! sensitivity reactions so that the ReactorNet can keep track of the
    //! order in which sensitivity parameters are added.
    void registerSensitivityReaction(void* reactor, size_t reactionIndex,
                                     const std::string& name, int leftright=0);

    //! The name of the p-th sensitivity parameter added to this ReactorNet.
    const std::string& sensitivityParameterName(size_t p) {
        return m_paramNames.at(p);
    }

    //! Reinitialize the integrator. Used to solve a new problem (different
    //! initial conditions) but with the same configuration of the reactor
    //! network. Can be called manually, or automatically after calling
    //! setInitialTime or modifying a reactor's contents.
    void reinitialize();

    //! Called to trigger integrator reinitialization before further
    //! integration.
    void setNeedsReinit() {
        m_integrator_init = false;
    }

protected:
    /**
     * Initialize the reactor network. Called automatically the first time
     * advance or step is called.
     */
    void initialize();

    std::vector<Reactor*> m_reactors;
    Integrator* m_integ;
    doublereal m_time;
    bool m_init;
    bool m_integrator_init; //! True if integrator initialization is current
    size_t m_nv;

    //! m_start[n] is the starting point in the state vector for reactor n
    std::vector<size_t> m_start;

    vector_fp m_atol;
    doublereal m_rtol, m_rtolsens;
    doublereal m_atols, m_atolsens;
    doublereal m_maxstep;
    int m_maxErrTestFails;
    bool m_verbose;
    size_t m_ntotpar;
    std::vector<size_t> m_nparams;

    //! Names corresponding to each sensitivity parameter
    std::vector<std::string> m_paramNames;

    //! Structure used to determine the order of sensitivity parameters
    //! m_sensOrder[Reactor or Wall, leftright][reaction number] = parameter index
    std::map<std::pair<void*, int>, std::map<size_t, size_t> > m_sensOrder;

    //! Mapping from the order in which sensitivity parameters were added to
    //! the ReactorNet to the order in which they occur in the integrator
    //! output.
    std::vector<size_t> m_sensIndex;

    vector_fp m_ydot;

    std::vector<bool> m_iown;
};
}

#endif
