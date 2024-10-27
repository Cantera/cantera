/**
 * @file Sim1D.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SIM1D_H
#define CT_SIM1D_H

#include "OneDim.h"

namespace Cantera
{

/**
 * One-dimensional simulations. Class Sim1D extends class OneDim by storing
 * the solution vector, and by adding a hybrid Newton/time-stepping solver.
 * @ingroup onedGroup
 */
class Sim1D : public OneDim
{
public:
    //! Default constructor.
    /*!
     * This constructor is provided to make the class default-constructible, but
     * is not meant to be used in most applications.  Use the next constructor
     */
    Sim1D() {}

    /**
     * Standard constructor.
     * @param domains  A vector of shared pointers to the domains to be linked together.
     *     The domain pointers must be entered in left-to-right order --- that is,
     *     the pointer to the leftmost domain is domain[0], the pointer to the
     *     domain to its right is domain[1], etc.
     */
    Sim1D(vector<shared_ptr<Domain1D>>& domains);

    //! @name Setting initial values
    //!
    //! These methods are used to set the initial values of solution components.
    //! @{

    //! Set initial guess for one component for all domains
    /**
     * @param component  component name
     * @param locs  A vector of relative positions, beginning with 0.0 at the
     *     left of the domain, and ending with 1.0 at the right of the domain.
     * @param vals  A vector of values corresponding to the relative position
     *     locations.
     */
    void setInitialGuess(const string& component, vector<double>& locs,
                         vector<double>& vals);

    /**
     * Set a single value in the solution vector.
     * @param dom  domain number, beginning with 0 for the leftmost domain.
     * @param comp  component number
     * @param localPoint  grid point within the domain, beginning with 0 for
     *     the leftmost grid point in the domain.
     * @param value  the value.
     */
    void setValue(size_t dom, size_t comp, size_t localPoint, double value);

    /**
     * Get one entry in the solution vector.
     * @param dom  domain number, beginning with 0 for the leftmost domain.
     * @param comp  component number
     * @param localPoint  grid point within the domain, beginning with 0 for
     *     the leftmost grid point in the domain.
     */
    double value(size_t dom, size_t comp, size_t localPoint) const;

    //! Get an entry in the work vector, which may contain either a new system state
    //! or the current residual of the system.
    //! @param dom  domain index
    //! @param comp  component index
    //! @param localPoint  grid point within the domain
    double workValue(size_t dom, size_t comp, size_t localPoint) const;

    /**
     * Specify a profile for one component of one domain.
     * @param dom  domain number, beginning with 0 for the leftmost domain.
     * @param comp  component number
     * @param pos  A vector of relative positions, beginning with 0.0 at the
     *     left of the domain, and ending with 1.0 at the right of the domain.
     * @param values  A vector of values corresponding to the relative position
     *     locations.
     *
     * Note that the vector pos and values can have lengths different than the
     * number of grid points, but their lengths must be equal. The values at
     * the grid points will be linearly interpolated based on the (pos,
     * values) specification.
     */
    void setProfile(size_t dom, size_t comp, const vector<double>& pos,
                    const vector<double>& values);

    //! Set component 'comp' of domain 'dom' to value 'v' at all points.
    void setFlatProfile(size_t dom, size_t comp, double v);

    //! @}

    //! @name Logging, saving and restoring of solutions
    //!
    //! @{

    /**
     * Output information on current solution for all domains to stream.
     * @param s  Output stream
     * @since New in %Cantera 3.0.
     * @deprecated To be removed after Cantera 3.1.
     */
    void show(std::ostream& s);

    /**
     * Show logging information on current solution for all domains.
     * @since New in %Cantera 3.0.
     */
    void show();

    /**
     * Save current simulation data to a container file or CSV format.
     *
     * In order to save the content of a Sim1D object, individual domains are
     * converted to SolutionArray objects and saved using the SolutionArray::save()
     * method. For HDF and YAML output, all domains are written to a single container
     * file with shared header information. Simulation settings of individual domains
     * are preserved as meta data of the corresponding SolutionArray objects.
     * For CSV files, only state and auxiliary data of the main 1D domain are saved.
     *
     * The complete state of the current object can be restored from HDF and YAML
     * container files using the restore() method, while individual domains can be
     * loaded using SolutionArray::restore() for further analysis. While CSV do not
     * contain complete information, they can still be used for setting initial states
     * of individual simulation objects for some %Cantera API's.
     *
     * @param fname  Name of output file (CSV, YAML or HDF)
     * @param name  Identifier of storage location within the container file; this
     *      node/group contains header information and multiple subgroups holding
     *      domain-specific SolutionArray data (YAML/HDF only)
     * @param desc  Custom comment describing the dataset to be stored (YAML/HDF only)
     * @param overwrite  Force overwrite if file/name exists; optional (default=false)
     * @param compression  Compression level (0-9); optional (default=0; HDF only)
     * @param basis  Output mass ("Y"/"mass") or mole ("X"/"mole") fractions;
     *      if not specified (default=""), the native basis of the underlying
     *      ThermoPhase manager is used - @see nativeState (CSV only)
     */
    void save(const string& fname, const string& name, const string& desc,
              bool overwrite=false, int compression=0, const string& basis="");

    /**
     * Save the residual of the current solution to a container file.
     * @param fname  Name of output container file
     * @param name  Identifier of solution within the container file
     * @param desc  Description of the solution
     * @param overwrite  Force overwrite if name exists; optional (default=false)
     * @param compression  Compression level (optional; HDF only)
     */
    void saveResidual(const string& fname, const string& name,
                      const string& desc, bool overwrite=false, int compression=0);

    /**
     * Retrieve data and settings from a previously saved simulation.
     *
     * This method restores a simulation object from YAML or HDF data previously saved
     * using the save() method.
     *
     * @param fname  Name of container file (YAML or HDF)
     * @param name  Identifier of location within the container file; this node/group
     *      contains header information and subgroups with domain-specific SolutionArray
     *      data
     * @return  AnyMap containing header information
     */
    AnyMap restore(const string& fname, const string& name);

    /**
     * Deletes a `debug_sim1d.yaml` file if it exists. Used to clear the file for
     * successive calls to the solve() method.
     */
    void clearDebugFile();

    /**
     * Write solver debugging information to a YAML file based on the specified log
     * level.
     *
     * This method writes solver debug information to a specified YAML file
     * (`debug_sim1d.yaml`). The section headers are formatted according to the provided
     * `header_suffix` and `attempt_counter` arguments. Depending on the log level, the
     * method will save either the solution information or the residual information
     * for each attempted solution.
     *
     * @param header_suffix  Header used to construct a unique section in the YAML file
     *                       where the information will be written to.
     * @param message  A string that is written to the `description` tag in the YAML
     *                 file.
     * @param loglevel  Controls the type of output that will be written. A `loglevel`
     *                  greater than 6 saves the solution, and a `loglevel` greater
     *                  than 7 saves the residual additionally.
     * @param attempt_counter  An integer counter used to uniquely identify the attempt
     *                         which is included in the file header to differentiate
     *                         between multiple solution attempts.
     */
     void writeDebugInfo(const string& header_suffix, const string& message, int loglevel,
                         int attempt_counter);


    //! @}

    //! Set the number of time steps to try when the steady Newton solver is
    //! unsuccessful.
    //! @param stepsize  Initial time step size [s]
    //! @param n  Length of `tsteps` array
    //! @param tsteps  A sequence of time step counts to take after subsequent failures
    //!     of the steady-state solver. The last value in `tsteps` will be used again
    //!     after further unsuccessful solution attempts.
    void setTimeStep(double stepsize, size_t n, const int* tsteps);

    /**
     * Performs the hybrid Newton steady/time-stepping solution.
     *
     * The solver attempts to solve the steady-state problem first. If the steady-state
     * solver fails, the time-stepping solver is used to take multiple time steps to
     * move the solution closer to the steady-state solution. The steady-state solver is
     * called again after the timesteps to make further progress towards the steady-state
     * solution. This process is repeated until the steady-state solver converges or the
     * maximum number of timesteps is reached.
     *
     * At the end of a successful solve, if the `refine_grid` flag is set, the grid will be
     * analyzed and refined if necessary. If the grid is refined, the solution process
     * described above is repeated with the new grid. This process is repeated until the
     * grid no longer needs refinement based on the refine criteria.
     *
     * @param loglevel  Controls the amount of diagnostic output.
     * @param refine_grid  If `true`, the grid will be refined
     */
    void solve(int loglevel = 0, bool refine_grid = true);

    void eval(double rdt=-1.0, int count = 1) {
        OneDim::eval(npos, m_state->data(), m_xnew.data(), rdt, count);
    }

    //! Evaluate the governing equations and return the vector of residuals
    void getResidual(double rdt, double* resid) {
        OneDim::eval(npos, m_state->data(), resid, rdt, 0);
    }

    //! Refine the grid in all domains.
    //!
    //! @returns If positive, the number of new grid points added. If negative, the
    //!     number of grid points removed. If zero, the grid is unchanged.
    //!
    //! @since Changed in %Cantera 3.1. Previously, the return value was zero if points
    //!     were removed but not added.
    int refine(int loglevel=0);

    //! Add node for fixed temperature point of freely propagating flame
    int setFixedTemperature(double t);

    //! Return temperature at the point used to fix the flame location
    double fixedTemperature();

    //! Return location of the point where temperature is fixed
    double fixedTemperatureLocation();

    /**
     * Set the left control point location using the specified temperature.
     * This is used when two-point flame control is active.
     *
     * The provided temperature will be used to locate the closest grid point to
     * that temperature, which will serve to locate the left control point's
     * coordinate. Starting from the left boundary, the first grid point that is
     * equal to or exceeds the specified temperature will be used to locate the
     * left control point's coordinate.
     */
    void setLeftControlPoint(double temperature);

    /**
     * Set the right control point location using the specified temperature.
     * This is used when two-point flame control is active.
     *
     * The provided temperature will be used to locate the closest grid point to
     * that temperature, which will serve to locate the right control point's
     * coordinate. Starting from the right boundary, the first grid point that is
     * equal to or exceeds the specified temperature will be used to locate the
     * right control point's coordinate.
     */
    void setRightControlPoint(double temperature);

    /**
     * Set grid refinement criteria. If dom >= 0, then the settings
     * apply only to the specified domain.  If dom < 0, the settings
     * are applied to each domain.  @see Refiner::setCriteria.
     */
    void setRefineCriteria(int dom = -1, double ratio = 10.0,
                           double slope = 0.8, double curve = 0.8,
                           double prune = -0.1);

    /**
     * Get the grid refinement criteria. dom must be greater than
     * or equal to zero (that is, the domain must be specified).
     * @see Refiner::getCriteria
     */
    vector<double> getRefineCriteria(int dom);

    /**
     * Set the maximum number of grid points in the domain. If dom >= 0,
     * then the settings apply only to the specified domain. If dom < 0,
     * the settings are applied to each domain.  @see Refiner::setMaxPoints.
     */
    void setMaxGridPoints(int dom, int npoints);

    /**
     * Get the maximum number of grid points in this domain. @see Refiner::maxPoints
     *
     * @param dom  domain number, beginning with 0 for the leftmost domain.
     */
    size_t maxGridPoints(size_t dom);

    //! Set the minimum grid spacing in the specified domain(s).
    /*!
     * @param dom  Domain index. If dom == -1, the specified spacing is applied
     *             to all domains.
     * @param gridmin  The minimum allowable grid spacing [m]
     */
    void setGridMin(int dom, double gridmin);

    //! Set the current solution vector to the last successful time-stepping
    //! solution. This can be used to examine the solver progress after a failed
    //! integration.
    void restoreTimeSteppingSolution();

    //! Set the current solution vector and grid to the last successful steady-
    //! state solution. This can be used to examine the solver progress after a
    //! failure during grid refinement.
    void restoreSteadySolution();

    //! Get the initial value of the system state from each domain in the simulation.
    void getInitialSoln();

    //! Get the Jacobian element @f$ J_{ij} = \partial f_i / \partial x_j \f$
    double jacobian(int i, int j);

    //! Evaluate the Jacobian in steady-state mode.
    void evalSSJacobian();

    //! Solve the equation @f$ J^T \lambda = b @f$.
    /**
     * Here, @f$ J = \partial f/\partial x @f$ is the Jacobian matrix of the
     * system of equations @f$ f(x,p)=0 @f$. This can be used to efficiently
     * solve for the sensitivities of a scalar objective function @f$ g(x,p) @f$
     * to a vector of parameters @f$ p @f$ by solving:
     * @f[ J^T \lambda = \left( \frac{\partial g}{\partial x} \right)^T @f]
     * for @f$ \lambda @f$ and then computing:
     * @f[
     *     \left.\frac{dg}{dp}\right|_{f=0} = \frac{\partial g}{\partial p}
     *         - \lambda^T \frac{\partial f}{\partial p}
     * @f]
     */
    void solveAdjoint(const double* b, double* lambda);

    void resize() override;

    //! Set a function that will be called after each successful steady-state
    //! solve, before regridding. Intended to be used for observing solver
    //! progress for debugging purposes.
    void setSteadyCallback(Func1* callback) {
        m_steady_callback = callback;
    }

protected:
    //! the solution vector after the last successful timestepping
    vector<double> m_xlast_ts;

    //! the solution vector after the last successful steady-state solve (stored
    //! before grid refinement)
    vector<double> m_xlast_ss;

    //! the grids for each domain after the last successful steady-state solve
    //! (stored before grid refinement)
    vector<vector<double>> m_grid_last_ss;

    //! a work array used to hold the residual or the new solution
    vector<double> m_xnew;

    //! timestep
    double m_tstep;

    //! array of number of steps to take before re-attempting the steady-state
    //! solution
    vector<int> m_steps;

    //! User-supplied function called after a successful steady-state solve.
    Func1* m_steady_callback;

private:
    //! Calls method _finalize in each domain.
    void finalize();

    //! Wrapper around the Newton solver
    /*!
     * @return 0 if successful, -1 on failure
     */
    int newtonSolve(int loglevel);
};

}
#endif
