/**
 * @file OneDim.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ONEDIM_H
#define CT_ONEDIM_H

#include "Domain1D.h"
#include "MultiJac.h"

namespace Cantera
{

class Func1;
class MultiNewton;
class AnyMap;

/**
 * Container class for multiple-domain 1D problems. Each domain is
 * represented by an instance of Domain1D.
 * @ingroup onedGroup
 */
class OneDim
{
public:
    OneDim();

    //! Construct a OneDim container for the domains in the list *domains*.
    OneDim(vector<shared_ptr<Domain1D>>& domains);

    virtual ~OneDim();
    OneDim(const OneDim&) = delete;
    OneDim& operator=(const OneDim&) = delete;

    //! Add a domain. Domains are added left-to-right.
    void addDomain(shared_ptr<Domain1D> d);

    //! Return a reference to the Jacobian evaluator of an OneDim object.
    //! @ingroup derivGroup
    MultiJac& jacobian();

    //! Return a reference to the Newton iterator.
    MultiNewton& newton();

    /**
     * Solve F(x) = 0, where F(x) is the multi-domain residual function.
     *
     * @param x0         Starting estimate of solution.
     * @param x1         Final solution satisfying F(x1) = 0.
     * @param loglevel   Controls amount of diagnostic output.
     *
     * @returns
     * - 1 for success
     * - -2 failure (maximum number of damping steps was reached)
     * - -3 failure (solution was up against the bounds
     */
    int solve(double* x0, double* x1, int loglevel);

    //! Number of domains.
    size_t nDomains() const {
        return m_dom.size();
    }

    //! Return a reference to domain i.
    Domain1D& domain(size_t i) const {
        return *m_dom[i];
    }

    size_t domainIndex(const string& name);

    //! Check that the specified domain index is in range.
    //! Throws an exception if n is greater than nDomains()-1
    void checkDomainIndex(size_t n) const {
        if (n >= m_dom.size()) {
            throw IndexError("OneDim::checkDomainIndex", "domains", n,
                             m_dom.size()-1);
        }
    }

    //! Check that an array size is at least nDomains().
    //! Throws an exception if nn is less than nDomains(). Used before calls
    //! which take an array pointer.
    void checkDomainArraySize(size_t nn) const {
        if (m_dom.size() > nn) {
            throw ArraySizeError("OneDim::checkDomainArraySize", nn,
                                 m_dom.size());
        }
    }

    //! The index of the start of domain i in the solution vector.
    size_t start(size_t i) const {
        if (m_dom[i]->nComponents()) {
            return m_dom[i]->loc();
        } else {
            // Special case for domains with no solution components to avoid
            // spurious out-of-bounds memory access
            return 0;
        }
    }

    //! Total solution vector length;
    size_t size() const {
        return m_size;
    }

    //! Pointer to left-most domain (first added).
    Domain1D* left() {
        return m_dom[0].get();
    }

    //! Pointer to right-most domain (last added).
    Domain1D* right() {
        return m_dom.back().get();
    }

    //! Number of solution components at global point jg.
    size_t nVars(size_t jg) {
        return m_nvars[jg];
    }

    //! Location in the solution vector of the first component of global point
    //! jg.
    size_t loc(size_t jg) {
        return m_loc[jg];
    }

    //! Return the domain, local point index, and component name for the i-th
    //! component of the global solution vector
    std::tuple<string, size_t, string> component(size_t i);

    //! Jacobian bandwidth.
    size_t bandwidth() const {
        return m_bw;
    }

    /**
     * Initialize all domains. On the first call, this methods calls the init
     * method of each domain, proceeding from left to right. Subsequent calls
     * do nothing.
     */
    void init();

    //! Total number of points.
    size_t points() {
        return m_pts;
    }

    /**
     * Steady-state max norm (infinity norm) of the residual evaluated using
     * solution x. On return, array r contains the steady-state residual
     * values. Used only for diagnostic output.
     */
    double ssnorm(double* x, double* r);

    //! Reciprocal of the time step.
    double rdt() const {
        return m_rdt;
    }

    //! Prepare for time stepping beginning with solution *x* and timestep *dt*.
    void initTimeInteg(double dt, double* x);

    //! True if transient mode.
    bool transient() const {
        return (m_rdt != 0.0);
    }

    //! True if steady mode.
    bool steady() const {
        return (m_rdt == 0.0);
    }

    /**
     * Prepare to solve the steady-state problem. After invoking this method,
     * subsequent calls to solve() will solve the steady-state problem. Sets
     * the reciprocal of the time step to zero, and, if it was previously non-
     * zero, signals that a new Jacobian will be needed.
     */
    void setSteadyMode();

    /**
     * Evaluate the multi-domain residual function
     *
     * @param j       if j != npos, only evaluate residual for points j-1, j,
     *                and j + 1; otherwise, evaluate at all grid points.
     * @param x       solution vector
     * @param r       on return, contains the residual vector
     * @param rdt     Reciprocal of the time step. if omitted, then
     *                  the default value is used.
     * @param count   Set to zero to omit this call from the statistics
     */
    void eval(size_t j, double* x, double* r, double rdt=-1.0, int count = 1);

    //! Return a pointer to the domain global point *i* belongs to.
    /*!
     * The domains are scanned right-to-left, and the first one with starting
     * location less or equal to i is returned.
     */
    Domain1D* pointDomain(size_t i);

    //! Call after one or more grids has changed size, for example after being refined.
    virtual void resize();

    vector<int>& transientMask() {
        return m_mask;
    }

    /**
     * Take time steps using Backward Euler.
     *
     * @param nsteps number of steps
     * @param dt initial step size
     * @param x current solution vector
     * @param r solution vector after time stepping
     * @param loglevel controls amount of printed diagnostics
     * @returns size of last timestep taken
     */
    double timeStep(int nsteps, double dt, double* x, double* r, int loglevel);

    void resetBadValues(double* x);

    //! Write statistics about the number of iterations and Jacobians at each
    //! grid level
    /*!
     *  @param printTime  Boolean that indicates whether time should be printed
     *                    out The default is true. It's turned off for test
     *                    problems where we don't want to print any times
     */
    void writeStats(int printTime = 1);

    // options
    void setMinTimeStep(double tmin) {
        m_tmin = tmin;
    }
    void setMaxTimeStep(double tmax) {
        m_tmax = tmax;
    }

    /**
     * Sets a factor by which the time step is reduced if the time stepping
     * fails. The default value is 0.5.
     *
     * @param tfactor factor time step is multiplied by if time stepping fails
     */
    void setTimeStepFactor(double tfactor) {
        m_tfactor = tfactor;
    }

    //! Set the maximum number of timeteps allowed before successful
    //! steady-state solve
    void setMaxTimeStepCount(int nmax) {
        m_nsteps_max = nmax;
    }

    //! Return the maximum number of timeteps allowed before successful
    //! steady-state solve
    int maxTimeStepCount() const {
        return m_nsteps_max;
    }

    void setJacAge(int ss_age, int ts_age=-1);

    /**
     * Save statistics on function and Jacobian evaluation, and reset the
     * counters. Statistics are saved only if the number of Jacobian
     * evaluations is greater than zero. The statistics saved are:
     *
     * - number of grid points
     * - number of Jacobian evaluations
     * - CPU time spent evaluating Jacobians
     * - number of non-Jacobian function evaluations
     * - CPU time spent evaluating functions
     * - number of time steps
     */
    void saveStats();

    //! Clear saved statistics
    void clearStats();

    //! Return total grid size in each call to solve()
    const vector<size_t>& gridSizeStats() {
        saveStats();
        return m_gridpts;
    }

    //! Return CPU time spent evaluating Jacobians in each call to solve()
    const vector<double>& jacobianTimeStats() {
        saveStats();
        return m_jacElapsed;
    }

    //! Return CPU time spent on non-Jacobian function evaluations in each call
    //! to solve()
    const vector<double>& evalTimeStats() {
        saveStats();
        return m_funcElapsed;
    }

    //! Return number of Jacobian evaluations made in each call to solve()
    const vector<int>& jacobianCountStats() {
        saveStats();
        return m_jacEvals;
    }

    //! Return number of non-Jacobian function evaluations made in each call to
    //! solve()
    const vector<int>& evalCountStats() {
        saveStats();
        return m_funcEvals;
    }

    //! Return number of time steps taken in each call to solve()
    const vector<int>& timeStepStats() {
        saveStats();
        return m_timeSteps;
    }

    //! Set a function that will be called every time #eval is called.
    //! Can be used to provide keyboard interrupt support in the high-level
    //! language interfaces.
    void setInterrupt(Func1* interrupt) {
        m_interrupt = interrupt;
    }

    //! Set a function that will be called after each successful timestep. The
    //! function will be called with the size of the timestep as the argument.
    //! Intended to be used for observing solver progress for debugging
    //! purposes.
    void setTimeStepCallback(Func1* callback) {
        m_time_step_callback = callback;
    }

protected:
    void evalSSJacobian(double* x, double* xnew);

    double m_tmin = 1e-16; //!< minimum timestep size
    double m_tmax = 1e+08; //!< maximum timestep size

    //! factor time step is multiplied by  if time stepping fails ( < 1 )
    double m_tfactor = 0.5;

    shared_ptr<vector<double>> m_state; //!< Solution vector

    unique_ptr<MultiJac> m_jac; //!< Jacobian evaluator
    unique_ptr<MultiNewton> m_newt; //!< Newton iterator
    double m_rdt = 0.0; //!< reciprocal of time step
    bool m_jac_ok = false; //!< if true, Jacobian is current

    size_t m_bw = 0; //!< Jacobian bandwidth
    size_t m_size = 0; //!< solution vector size

    vector<shared_ptr<Domain1D>> m_dom;
    vector<shared_ptr<Domain1D>> m_connect;
    vector<shared_ptr<Domain1D>> m_bulk;

    bool m_init = false;
    vector<size_t> m_nvars;
    vector<size_t> m_loc;
    vector<int> m_mask;
    size_t m_pts = 0;

    // options
    int m_ss_jac_age = 20;
    int m_ts_jac_age = 20;

    //! Function called at the start of every call to #eval.
    Func1* m_interrupt = nullptr;

    //! User-supplied function called after each successful timestep.
    Func1* m_time_step_callback = nullptr;

    //! Number of time steps taken in the current call to solve()
    int m_nsteps = 0;

    //! Maximum number of timesteps allowed per call to solve()
    int m_nsteps_max = 5000;

private:
    // statistics
    int m_nevals = 0;
    double m_evaltime = 0;
    vector<size_t> m_gridpts;
    vector<int> m_jacEvals;
    vector<double> m_jacElapsed;
    vector<int> m_funcEvals;
    vector<double> m_funcElapsed;

    //! Number of time steps taken in each call to solve() (for example, for each
    //! successive grid refinement)
    vector<int> m_timeSteps;
};

}

#endif
