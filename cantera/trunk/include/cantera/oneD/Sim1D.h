/**
 * @file Sim1D.h
 */

#ifndef CT_SIM1D_H
#define CT_SIM1D_H

#include "OneDim.h"
#include "cantera/numerics/funcs.h"

namespace Cantera
{

/**
 * One-dimensional simulations. Class Sim1D extends class OneDim by storing
 * the solution vector, and by adding a hybrid Newton/time-stepping solver.
 * @ingroup onedim
 */
class Sim1D : public OneDim
{
public:

    //! Default constructor.
    /*!
     *  This constructor is provided to make the class default-constructible,
     *  but is not meant to be used in most applications.  Use the next
     *  constructor
     */
    Sim1D();

    /**
     * Standard constructor.
     * @param domains A vector of pointers to the domains to be linked together.
     * The domain pointers must be entered in left-to-right order --- i.e.,
     * the pointer to the leftmost domain is domain[0], the pointer to the
     * domain to its right is domain[1], etc.
     */
    Sim1D(std::vector<Domain1D*>& domains);

    /**
     * @name Setting initial values
     *
     * These methods are used to set the initial values of
     * solution components.
     */
    //@{

    /// Set initial guess based on equilibrium
    void setInitialGuess(const std::string& component, vector_fp& locs,
                         vector_fp& vals);

    /**
     * Set a single value in the solution vector.
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param localPoint grid point within the domain, beginning with 0 for
     *     the leftmost grid point in the domain.
     * @param value the value.
     */
    void setValue(size_t dom, size_t comp, size_t localPoint,  doublereal value);

    /**
     * Get one entry in the solution vector.
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param localPoint grid point within the domain, beginning with 0 for
     *     the leftmost grid point in the domain.
     */
    doublereal value(size_t dom, size_t comp, size_t localPoint) const;

    doublereal workValue(size_t dom, size_t comp, size_t localPoint) const;

    /**
     * Specify a profile for one component of one domain.
     * @param dom domain number, beginning with 0 for the leftmost domain.
     * @param comp component number
     * @param pos A vector of relative positions, beginning with 0.0 at the
     *     left of the domain, and ending with 1.0 at the right of the domain.
     * @param values A vector of values corresponding to the relative position
     *     locations.
     *
     * Note that the vector pos and values can have lengths different than the
     * number of grid points, but their lengths must be equal. The values at
     * the grid points will be linearly interpolated based on the (pos,
     * values) specification.
     */
    void setProfile(size_t dom, size_t comp, const vector_fp& pos,
                    const vector_fp& values);

    /// Set component 'comp' of domain 'dom' to value 'v' at all points.
    void setFlatProfile(size_t dom, size_t comp, doublereal v);

    //@}

    void save(const std::string& fname, const std::string& id,
              const std::string& desc, int loglevel=1);

    void saveResidual(const std::string& fname, const std::string& id,
                      const std::string& desc, int loglevel=1);

    /// Print to stream s the current solution for all domains.
    void showSolution(std::ostream& s);
    void showSolution();

    const doublereal* solution() {
        return DATA_PTR(m_x);
    }

    void setTimeStep(doublereal stepsize, size_t n, integer* tsteps);

    //void setMaxTimeStep(doublereal tmax) { m_maxtimestep = tmax; }

    void solve(int loglevel = 0, bool refine_grid = true);

    void eval(doublereal rdt=-1.0, int count = 1) {
        OneDim::eval(npos, DATA_PTR(m_x), DATA_PTR(m_xnew), rdt, count);
    }

    /// Refine the grid in all domains.
    int refine(int loglevel=0);

    //! Add node for fixed temperature point of freely propagating flame
    int setFixedTemperature(doublereal t);

    void setAdiabaticFlame(void);

    /**
     * Set grid refinement criteria. If dom >= 0, then the settings
     * apply only to the specified domain.  If dom < 0, the settings
     * are applied to each domain.  @see Refiner::setCriteria.
     */
    void setRefineCriteria(int dom = -1, doublereal ratio = 10.0,
                           doublereal slope = 0.8, doublereal curve = 0.8, doublereal prune = -0.1);
    void setMaxGridPoints(int dom = -1, int npoints = 300);

    //! Set the minimum grid spacing in the specified domain(s).
    /*!
     *  @param dom Domain index. If dom == -1, the specified spacing
                   is applied to all domains.
        @param gridmin The minimum allowable grid spacing [m]
    */
    void setGridMin(int dom, double gridmin);

    //! Initialize the solution with a previously-saved solution.
    void restore(const std::string& fname, const std::string& id, int loglevel=2);

    void getInitialSoln();

    void setSolution(const doublereal* soln) {
        std::copy(soln, soln + m_x.size(), DATA_PTR(m_x));
    }

    const doublereal* solution() const {
        return DATA_PTR(m_x);
    }

    doublereal jacobian(int i, int j);

    void evalSSJacobian();

protected:
    //! the solution vector
    vector_fp m_x;

    //! a work array used to hold the residual or the new solution
    vector_fp m_xnew;

    //! timestep
    doublereal m_tstep;

    //! array of number of steps to take before re-attempting the steady-state
    //! solution
    vector_int m_steps;

private:
    /// Calls method _finalize in each domain.
    void finalize();

    /*! Wrapper around the Newton solver.
     * @return 0 if successful, -1 on failure
     */
    int newtonSolve(int loglevel);
};

}
#endif
