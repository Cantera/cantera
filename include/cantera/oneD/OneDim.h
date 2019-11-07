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

/**
 * Container class for multiple-domain 1D problems. Each domain is
 * represented by an instance of Domain1D.
 * @ingroup onedim
 */
class OneDim
{
public:
    OneDim();

    //! Construct a OneDim container for the domains in the list *domains*.
    OneDim(std::vector<Domain1D*> domains);
    virtual ~OneDim();
    OneDim(const OneDim&) = delete;
    OneDim& operator=(const OneDim&) = delete;

    /// Add a domain. Domains are added left-to-right.
    void addDomain(Domain1D* d);

    //! Return a reference to the Jacobian evaluator.
    MultiJac& jacobian();

    /// Return a reference to the Newton iterator.
    MultiNewton& newton();

    /**
     * Solve F(x) = 0, where F(x) is the multi-domain residual function.
     * @param x0         Starting estimate of solution.
     * @param x1         Final solution satisfying F(x1) = 0.
     * @param loglevel   Controls amount of diagnostic output.
     */
    int solve(doublereal* x0, doublereal* x1, int loglevel);

    /// Number of domains.
    size_t nDomains() const {
        return m_dom.size();
    }

    /// Return a reference to domain i.
    Domain1D& domain(size_t i) const {
        return *m_dom[i];
    }

    size_t domainIndex(const std::string& name);

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

    /// The index of the start of domain i in the solution vector.
    size_t start(size_t i) const {
        return m_dom[i]->loc();
    }

    /// Total solution vector length;
    size_t size() const {
        return m_size;
    }

    /// Pointer to left-most domain (first added).
    Domain1D* left() {
        return m_dom[0];
    }

    /// Pointer to right-most domain (last added).
    Domain1D* right() {
        return m_dom.back();
    }

    /// Number of solution components at global point jg.
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
    std::tuple<std::string, size_t, std::string> component(size_t i);

    /// Jacobian bandwidth.
    size_t bandwidth() const {
        return m_bw;
    }

    /*!
     * Initialize all domains. On the first call, this methods calls the init
     * method of each domain, proceeding from left to right. Subsequent calls
     * do nothing.
     */
    void init();

    /// Total number of points.
    size_t points() {
        return m_pts;
    }

    /**
     * Steady-state max norm (infinity norm) of the residual evaluated using
     * solution x. On return, array r contains the steady-state residual
     * values. Used only for diagnostic output.
     */
    doublereal ssnorm(doublereal* x, doublereal* r);

    /// Reciprocal of the time step.
    doublereal rdt() const {
        return m_rdt;
    }

    //! Prepare for time stepping beginning with solution *x* and timestep *dt*.
    void initTimeInteg(doublereal dt, doublereal* x);

    /// True if transient mode.
    bool transient() const {
        return (m_rdt != 0.0);
    }

    /// True if steady mode.
    bool steady() const {
        return (m_rdt == 0.0);
    }

    /*!
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
    void eval(size_t j, double* x, double* r, doublereal rdt=-1.0,
              int count = 1);

    //! Return a pointer to the domain global point *i* belongs to.
    /*!
     * The domains are scanned right-to-left, and the first one with starting
     * location less or equal to i is returned.
     */
    Domain1D* pointDomain(size_t i);

    //! Call after one or more grids has changed size, e.g. after being refined.
    virtual void resize();

    vector_int& transientMask() {
        return m_mask;
    }

    /*!
     * Take time steps using Backward Euler.
     *
     * @param nsteps number of steps
     * @param dt initial step size
     * @param x current solution vector
     * @param r solution vector after time stepping
     * @param loglevel controls amount of printed diagnostics
     * @returns size of last timestep taken
     */
    double timeStep(int nsteps, double dt, double* x,
                    double* r, int loglevel);

    void resetBadValues(double* x);

    //! Write statistics about the number of iterations and Jacobians at each
    //! grid level
    /*!
     *  @param printTime  Boolean that indicates whether time should be printed
     *                    out The default is true. It's turned off for test
     *                    problems where we don't want to print any times
     */
    void writeStats(int printTime = 1);

    void save(const std::string& fname, std::string id,
              const std::string& desc, doublereal* sol, int loglevel);

    // options
    void setMinTimeStep(doublereal tmin) {
        m_tmin = tmin;
    }
    void setMaxTimeStep(doublereal tmax) {
        m_tmax = tmax;
    }
    void setTimeStepFactor(doublereal tfactor) {
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
    const std::vector<size_t>& gridSizeStats() {
        saveStats();
        return m_gridpts;
    }

    //! Return CPU time spent evaluating Jacobians in each call to solve()
    const vector_fp& jacobianTimeStats() {
        saveStats();
        return m_jacElapsed;
    }

    //! Return CPU time spent on non-Jacobian function evaluations in each call
    //! to solve()
    const vector_fp& evalTimeStats() {
        saveStats();
        return m_funcElapsed;
    }

    //! Return number of Jacobian evaluations made in each call to solve()
    const vector_int& jacobianCountStats() {
        saveStats();
        return m_jacEvals;
    }

    //! Return number of non-Jacobian function evaluations made in each call to
    //! solve()
    const vector_int& evalCountStats() {
        saveStats();
        return m_funcEvals;
    }

    //! Return number of time steps taken in each call to solve()
    const vector_int& timeStepStats() {
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
    void evalSSJacobian(doublereal* x, doublereal* xnew);

    doublereal m_tmin; //!< minimum timestep size
    doublereal m_tmax; //!< maximum timestep size

    //! factor time step is multiplied by  if time stepping fails ( < 1 )
    doublereal m_tfactor;

    std::unique_ptr<MultiJac> m_jac; //!< Jacobian evaluator
    std::unique_ptr<MultiNewton> m_newt; //!< Newton iterator
    doublereal m_rdt; //!< reciprocal of time step
    bool m_jac_ok; //!< if true, Jacobian is current

    size_t m_bw; //!< Jacobian bandwidth
    size_t m_size; //!< solution vector size

    std::vector<Domain1D*> m_dom, m_connect, m_bulk;

    bool m_init;
    std::vector<size_t> m_nvars;
    std::vector<size_t> m_loc;
    vector_int m_mask;
    size_t m_pts;
    doublereal m_solve_time;

    // options
    int m_ss_jac_age, m_ts_jac_age;

    //! Function called at the start of every call to #eval.
    Func1* m_interrupt;

    //! User-supplied function called after each successful timestep.
    Func1* m_time_step_callback;

    //! Number of time steps taken in the current call to solve()
    int m_nsteps;

    //! Maximum number of timesteps allowed per call to solve()
    int m_nsteps_max;

private:
    // statistics
    int m_nevals;
    doublereal m_evaltime;
    std::vector<size_t> m_gridpts;
    vector_int m_jacEvals;
    vector_fp m_jacElapsed;
    vector_int m_funcEvals;
    vector_fp m_funcElapsed;

    //! Number of time steps taken in each call to solve() (e.g. for each
    //! successive grid refinement)
    vector_int m_timeSteps;
};

}

#endif
