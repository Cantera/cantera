/**
 * @file OneDim.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ONEDIM_H
#define CT_ONEDIM_H

#include "Domain1D.h"
#include "MultiJac.h"
#include "cantera/numerics/SystemJacobian.h"
#include "cantera/numerics/SteadyStateSystem.h"

namespace Cantera
{

class AnyMap;

/**
 * Container class for multiple-domain 1D problems. Each domain is
 * represented by an instance of Domain1D.
 * @ingroup onedGroup
 */
class OneDim : public SteadyStateSystem
{
public:
    //! Default constructor
    OneDim() = default;

    //! Construct a OneDim container for the domains in the list *domains*.
    OneDim(vector<shared_ptr<Domain1D>>& domains);

    //! Add a domain. Domains are added left-to-right.
    void addDomain(shared_ptr<Domain1D> d);

    //! @deprecated To be removed after Cantera 3.2. Use linearSolver() instead.
    shared_ptr<SystemJacobian> getJacobian() {
        warn_deprecated("OneDim::getJacobian",
                        "To be removed after Cantera 3.2. Use linearSolver() instead.");
        return m_jac;
    }

    //! Compute the weighted norm of a step vector
    //!
    //! The weighted norm of a step vector @f$ \Delta x @f$ is defined as
    //! @f[
    //!   ||\Delta x|| = \sqrt{ \frac{1}{N}
    //!       \sum_{d,n,j} \left( \frac{\Delta x_{d,n,j}}{w_{d,n}} \right)^2 }
    //! @f]
    //! where the error weight for solution component @f$ n @f$ in domain @f$ d @f$ is
    //! given by
    //! @f[
    //!   w_{d,n} = \frac{\epsilon_{{\rm r};d,n}}{J_d} \sum_j |x_{d,n,j}|
    //!           + \epsilon_{{\rm a};d,n}
    //! @f]
    //! Here, @f$ \epsilon_{{\rm r};d,n} @f$ is the relative error tolerance for
    //! component @f$ n @f$ in domain @f$ d @f$, and multiplies the average magnitude of
    //! solution component @f$ n @f$ in the domain. The second term,
    //! @f$ \epsilon_{{\rm a};d,n} @f$, is the absolute error tolerance for component
    //! @f$ n @f$ in domain @f$ d @f$. @f$ N @f$ is the total number of state variables
    //! across all domains and @f$ J_d @f$ is the number of grid points in domain
    //! @f$ d @f$.
    double weightedNorm(const double* step) const override;

    //! Return a reference to the Jacobian evaluator of an OneDim object.
    //! @deprecated To be removed after Cantera 3.2. Superseded by linearSolver()
    //! @ingroup derivGroup
    MultiJac& jacobian();

    //! Number of domains.
    size_t nDomains() const {
        return m_dom.size();
    }

    //! Return a reference to domain i.
    Domain1D& domain(size_t i) const {
        return *m_dom[i];
    }

    //! Get the index of the domain named `name`.
    size_t domainIndex(const string& name);

    //! Check that the specified domain index is in range.
    //! Throws an exception if n is greater than nDomains()-1
    void checkDomainIndex(size_t n) const {
        if (n >= m_dom.size()) {
            throw IndexError("OneDim::checkDomainIndex", "domains", n, m_dom.size());
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

    //! Pointer to left-most domain (first added).
    Domain1D* left() {
        return m_dom[0].get();
    }

    //! Pointer to right-most domain (last added).
    Domain1D* right() {
        return m_dom.back().get();
    }

    //! Number of solution components at global point `jg`.
    size_t nVars(size_t jg) {
        return m_nvars[jg];
    }

    //! Location in the solution vector of the first component of global point `jg`.
    size_t loc(size_t jg) {
        return m_loc[jg];
    }

    //! Return the domain, local point index, and component name for the i-th
    //! component of the global solution vector
    std::tuple<string, size_t, string> component(size_t i) const;

    string componentName(size_t i) const override;
    pair<string, string> componentTableHeader() const override;
    string componentTableLabel(size_t i) const override;
    double upperBound(size_t i) const override;
    double lowerBound(size_t i) const override;

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

    void initTimeInteg(double dt, double* x) override;
    void setSteadyMode() override;

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
    void eval(size_t j, double* x, double* r, double rdt=-1.0, int count=1);

    void eval(double* x, double* r, double rdt=-1.0, int count=1) override {
        return eval(npos, x, r, rdt, count);
    }

    void evalJacobian(double* x0) override;

    //! Return a pointer to the domain global point *i* belongs to.
    /*!
     * The domains are scanned right-to-left, and the first one with starting
     * location less or equal to i is returned.
     */
    Domain1D* pointDomain(size_t i);

    //! Call after one or more grids has changed size, for example after being refined.
    virtual void resize();

    void resetBadValues(double* x) override;

    //! Write statistics about the number of iterations and Jacobians at each
    //! grid level
    /*!
     *  @param printTime  Boolean that indicates whether time should be printed
     *                    out The default is true. It's turned off for test
     *                    problems where we don't want to print any times
     */
    void writeStats(int printTime = 1);

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

protected:
    //! All domains comprising the system
    vector<shared_ptr<Domain1D>> m_dom;

    //! All connector and boundary domains
    vector<shared_ptr<Domain1D>> m_connect;

    //! All bulk/flow domains
    vector<shared_ptr<Domain1D>> m_bulk;

    //! Indicates whether one-time initialization for each domain has been completed.
    bool m_init = false;

    //! Number of variables at each point, across all domains. Length points().
    //! Accessed with nVars().
    vector<size_t> m_nvars;

    //! Location in the state vector of the first component of each point, across all
    //! domains. Accessed with loc().
    vector<size_t> m_loc;

    //! Total number of points.
    size_t m_pts = 0;

private:
    //! @name Statistics
    //! Solver stats are collected after successfully solving on a particular grid.
    //! @{
    int m_nevals = 0; //!< Number of calls to eval()
    double m_evaltime = 0; //!< Total time [s] spent in eval()

    vector<size_t> m_gridpts; //!< Number of grid points in this grid
    vector<int> m_jacEvals; //!< Number of Jacobian evaluations on this grid
    vector<double> m_jacElapsed; //!< Time [s] spent evaluating Jacobians on this grid

    //! Number of residual function evaluations on this grid (not counting evaluations
    //! used to construct Jacobians).
    vector<int> m_funcEvals;

    //! Time [s] spent on residual function evaluations on this grid (not counting
    //! evaluations used to construct Jacobians).
    vector<double> m_funcElapsed;

    //! Number of time steps taken in each call to solve() (for example, for each
    //! successive grid refinement)
    vector<int> m_timeSteps;

    //! @}
};

}

#endif
