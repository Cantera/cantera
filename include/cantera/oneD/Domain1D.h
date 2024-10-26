 //! @file Domain1D.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DOMAIN1D_H
#define CT_DOMAIN1D_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

class MultiJac;
class OneDim;
class Refiner;
class AnyMap;
class Kinetics;
class Transport;
class Solution;
class SolutionArray;

/**
 * Base class for one-dimensional domains.
 * @ingroup flowGroup
 */
class Domain1D
{
public:
    /**
     * Constructor.
     * @param nv      Number of variables at each grid point.
     * @param points  Number of grid points.
     * @param time    (unused)
     */
    Domain1D(size_t nv=1, size_t points=1, double time=0.0);

    virtual ~Domain1D();
    Domain1D(const Domain1D&) = delete;
    Domain1D& operator=(const Domain1D&) = delete;

    //! Domain type flag.
    //! @since Starting in %Cantera 3.1, the return type is a `string`.
    virtual string domainType() const { return "domain"; }

    //! String indicating the domain implemented.
    //! @since New in %Cantera 3.0.
    //! @deprecated Transitional method. Use domainType() instead.
    string type() const { return domainType(); }

    //! The left-to-right location of this domain.
    size_t domainIndex() {
        return m_index;
    }

    //! True if the domain is a connector domain.
    virtual bool isConnector() {
        return false;
    }

    //! Set the solution manager.
    //! @since New in %Cantera 3.0.
    void setSolution(shared_ptr<Solution> sol);

    //! Set the kinetics manager.
    //! @since New in %Cantera 3.0.
    virtual void setKinetics(shared_ptr<Kinetics> kin) {
        throw NotImplementedError("Domain1D::setKinetics");
    }

    //! Set transport model to existing instance
    //! @since New in %Cantera 3.0.
    virtual void setTransport(shared_ptr<Transport> trans) {
        throw NotImplementedError("Domain1D::setTransport");
    }

    //! The container holding this domain.
    const OneDim& container() const {
        return *m_container;
    }

    //! Specify the container object for this domain, and the position of this
    //! domain in the list.
    void setContainer(OneDim* c, size_t index) {
        m_container = c;
        m_index = index;
    }

    //! Set the Jacobian bandwidth. See the discussion of method bandwidth().
    void setBandwidth(int bw = -1) {
        m_bw = bw;
    }

    //! Set the Jacobian bandwidth for this domain.
    /**
     * When class OneDim computes the bandwidth of the overall multi-domain
     * problem (in OneDim::resize()), it calls this method for the bandwidth
     * of each domain. If setBandwidth has not been called, then a negative
     * bandwidth is returned, in which case OneDim assumes that this domain is
     * dense -- that is, at each point, all components depend on the value of
     * all other components at that point. In this case, the bandwidth is bw =
     * 2*nComponents() - 1. However, if this domain contains some components
     * that are uncoupled from other components at the same point, then this
     * default bandwidth may greatly overestimate the true bandwidth, with a
     * substantial penalty in performance. For such domains, use method
     * setBandwidth to specify the bandwidth before passing this domain to the
     * Sim1D or OneDim constructor.
     */
    size_t bandwidth() {
        return m_bw;
    }

    /**
     * Initialize. This method is called by OneDim::init() for each domain once
     * at the beginning of a simulation. Base class method does nothing, but may
     * be overloaded.
     */
    virtual void init() {  }

    //! @deprecated Unused. To be removed after Cantera 3.1.
    virtual void setInitialState(double* xlocal = 0) {
        warn_deprecated("Domain1D::setInitialState", "To be removed after Cantera 3.1.");
    }

    //! @deprecated Unused. To be removed after Cantera 3.1.
    virtual void setState(size_t point, const double* state, double* x) {
        warn_deprecated("Domain1D::setState", "To be removed after Cantera 3.1.");
    }

    /**
     * When called, this function should reset "bad" values in the state vector
     * such as negative species concentrations. This function may be called
     * after a failed solution attempt.
     */
    virtual void resetBadValues(double* xg) {}

    /**
     * Resize the domain to have nv components and np grid points. This method
     * is virtual so that subclasses can perform other actions required to
     * resize the domain.
     */
    virtual void resize(size_t nv, size_t np);

    //! Return a reference to the grid refiner.
    Refiner& refiner() {
        return *m_refiner;
    }

    //! Number of components at each grid point.
    size_t nComponents() const {
        return m_nv;
    }

    //! Check that the specified component index is in range.
    //! Throws an exception if n is greater than nComponents()-1
    void checkComponentIndex(size_t n) const {
        if (n >= m_nv) {
            throw IndexError("Domain1D::checkComponentIndex", "points", n, m_nv-1);
        }
    }

    //! Check that an array size is at least nComponents().
    //! Throws an exception if nn is less than nComponents(). Used before calls
    //! which take an array pointer.
    void checkComponentArraySize(size_t nn) const {
        if (m_nv > nn) {
            throw ArraySizeError("Domain1D::checkComponentArraySize", nn, m_nv);
        }
    }

    //! Number of grid points in this domain.
    size_t nPoints() const {
        return m_points;
    }

    //! Check that the specified point index is in range.
    //! Throws an exception if n is greater than nPoints()-1
    void checkPointIndex(size_t n) const {
        if (n >= m_points) {
            throw IndexError("Domain1D::checkPointIndex", "points", n, m_points-1);
        }
    }

    //! Check that an array size is at least nPoints().
    //! Throws an exception if nn is less than nPoints(). Used before calls
    //! which take an array pointer.
    void checkPointArraySize(size_t nn) const {
        if (m_points > nn) {
            throw ArraySizeError("Domain1D::checkPointArraySize", nn, m_points);
        }
    }

    //! Name of component `n`. May be overloaded.
    virtual string componentName(size_t n) const;

    //! Set the name of the component `n` to `name`.
    void setComponentName(size_t n, const string& name) {
        m_name[n] = name;
    }

    //! index of component with name `name`.
    virtual size_t componentIndex(const string& name) const;

    /**
     * Set the upper and lower bounds for a solution component, n.
     *
     * @param n  solution component index
     * @param lower  lower bound on component n
     * @param upper  upper bound on component n
     */
    void setBounds(size_t n, double lower, double upper) {
        m_min[n] = lower;
        m_max[n] = upper;
    }

    //! Set tolerances for time-stepping mode
    /*!
     * @param rtol  Relative tolerance
     * @param atol  Absolute tolerance
     * @param n  component index these tolerances apply to. If set to -1 (the
     *      default), these tolerances will be applied to all solution
     *      components.
     */
    void setTransientTolerances(double rtol, double atol, size_t n=npos);

    //! Set tolerances for steady-state mode
    /*!
     * @param rtol  Relative tolerance
     * @param atol  Absolute tolerance
     * @param n  component index these tolerances apply to. If set to -1 (the
     *     default), these tolerances will be applied to all solution
     *     components.
     */
    void setSteadyTolerances(double rtol, double atol, size_t n=npos);

    //! Relative tolerance of the nth component.
    double rtol(size_t n) {
        return (m_rdt == 0.0 ? m_rtol_ss[n] : m_rtol_ts[n]);
    }

    //! Absolute tolerance of the nth component.
    double atol(size_t n) {
        return (m_rdt == 0.0 ? m_atol_ss[n] : m_atol_ts[n]);
    }

    //! Steady relative tolerance of the nth component
    double steady_rtol(size_t n) {
        return m_rtol_ss[n];
    }

    //! Steady absolute tolerance of the nth component
    double steady_atol(size_t n) {
        return m_atol_ss[n];
    }

    //! Transient relative tolerance of the nth component
    double transient_rtol(size_t n) {
        return m_rtol_ts[n];
    }

    //! Transient absolute tolerance of the nth component
    double transient_atol(size_t n) {
        return m_atol_ts[n];
    }

    //! Upper bound on the nth component.
    double upperBound(size_t n) const {
        return m_max[n];
    }

    //! Lower bound on the nth component
    double lowerBound(size_t n) const {
        return m_min[n];
    }

    /**
     * Performs the setup required before starting a time-stepping solution.
     * Stores the solution provided in `x0` to the internal storage, and sets
     * the reciprocal of the time step to `1/dt`.
     *
     * @param[in] dt  Time step
     * @param[in] x0  Array to store the solution at the last time step
     */
    void initTimeInteg(double dt, const double* x0) {
        std::copy(x0 + loc(), x0 + loc() + size(), m_slast.begin());
        m_rdt = 1.0/dt;
    }

    /**
     * Set the internally-stored reciprocal of the time step to 0.0, which is
     * used to indicate that the problem is in steady-state mode.
     */
    void setSteadyMode() {
        m_rdt = 0.0;
    }

    //! True if in steady-state mode
    bool steady() {
        return (m_rdt == 0.0);
    }

    //! True if not in steady-state mode
    bool transient() {
        return (m_rdt != 0.0);
    }

    /**
     * Set this if something has changed in the governing
     * equations (for example, the value of a constant has been changed,
     * so that the last-computed Jacobian is no longer valid.
     */
    void needJacUpdate();

    //! Evaluate the residual function at point j. If j == npos,
    //! evaluate the residual function at all points.
    /*!
     *  This function must be implemented in classes derived from Domain1D.
     *
     *  @param[in] j  Grid point at which to update the residual
     *  @param[in] x  State vector
     *  @param[out] r  residual vector
     *  @param[out] mask  Boolean mask indicating whether each solution
     *      component has a time derivative (1) or not (0).
     *  @param[in] rdt  Reciprocal of the timestep (`rdt=0` implies steady-state.)
     */
    virtual void eval(size_t j, double* x, double* r, integer* mask, double rdt=0.0) {
        throw NotImplementedError("Domain1D::eval");
    }

    /**
     * Returns the index of the solution vector, which corresponds to component
     * n at grid point j.
     *
     * @param n  component index
     * @param j  grid point index
     */
    size_t index(size_t n, size_t j) const {
        return m_nv*j + n;
    }

    /**
     * Returns the value of solution component n at grid point j of the solution
     * vector x.
     *
     * @param x  solution vector
     * @param n  component index
     * @param j  grid point index
     */
    double value(const double* x, size_t n, size_t j) const {
        return x[index(n,j)];
    }

    //! @deprecated To be removed after Cantera 3.1.
    virtual void setJac(MultiJac* jac) {
        warn_deprecated("Domain1D::setJac", "To be removed after Cantera 3.1.");
    }

    //! Save the state of this domain as a SolutionArray.
    /*!
     * @param soln  local solution vector for this domain
     * @todo  Despite the method's name, data are copied; the intent is to access data
     *      directly in future revisions, where a non-const version will be implemented.
     *
     * @since New in %Cantera 3.0.
     */
    virtual shared_ptr<SolutionArray> asArray(const double* soln) const {
        throw NotImplementedError("Domain1D::asArray", "Needs to be overloaded.");
    }

    //! Save the state of this domain to a SolutionArray.
    /*!
     * This method serves as an external interface for high-level API's; it does not
     * provide direct access to memory.
     * @param normalize  If true, normalize concentrations (default=false)
     *
     * @since New in %Cantera 3.0.
     */
    shared_ptr<SolutionArray> toArray(bool normalize=false) const;

    //! Restore the solution for this domain from a SolutionArray
    /*!
     * @param[in] arr  SolutionArray defining the state of this domain
     * @param[out] soln  Value of the solution vector, local to this domain
     *
     * @since New in %Cantera 3.0.
     */
    virtual void fromArray(SolutionArray& arr, double* soln) {
        throw NotImplementedError("Domain1D::fromArray", "Needs to be overloaded.");
    }

    //! Restore the solution for this domain from a SolutionArray.
    /*!
     * This method serves as an external interface for high-level API's.
     * @param arr  SolutionArray defining the state of this domain
     * @since New in %Cantera 3.0.
     */
    void fromArray(const shared_ptr<SolutionArray>& arr);

    //! Return thermo/kinetics/transport manager used in the domain
    //! @since New in %Cantera 3.0.
    shared_ptr<Solution> solution() const {
        return m_solution;
    }

    //! Return the size of the solution vector (the product of #m_nv and #m_points).
    size_t size() const {
        return m_nv*m_points;
    }

    /**
     * Find the index of the first grid point in this domain, and
     * the start of its variables in the global solution vector.
     */
    void locate();

    //! Location of the start of the local solution vector in the global solution vector
    virtual size_t loc(size_t j = 0) const {
        return m_iloc;
    }

    //! The index of the first (that is, left-most) grid point belonging to this domain.
    size_t firstPoint() const {
        return m_jstart;
    }

    //! The index of the last (that is, right-most) grid point belonging to this domain.
    size_t lastPoint() const {
        return m_jstart + m_points - 1;
    }

    /**
     * Set the left neighbor to domain 'left.' Method 'locate' is called to
     * update the global positions of this domain and all those to its right.
     */
    void linkLeft(Domain1D* left) {
        m_left = left;
        if (!m_solution && left && left->solution()) {
            setSolution(left->solution());
        }
        locate();
    }

    //! Set the right neighbor to domain 'right.'
    void linkRight(Domain1D* right) {
        m_right = right;
        if (!m_solution && right && right->solution()) {
            setSolution(right->solution());
        }
    }

    //! Append domain 'right' to this one, and update all links.
    void append(Domain1D* right) {
        linkRight(right);
        right->linkLeft(this);
    }

    //! Return a pointer to the left neighbor.
    Domain1D* left() const {
        return m_left;
    }

    //! Return a pointer to the right neighbor.
    Domain1D* right() const {
        return m_right;
    }

    //! Value of component n at point j in the previous solution.
    double prevSoln(size_t n, size_t j) const {
        return m_slast[m_nv*j + n];
    }

    //! Specify an identifying tag for this domain.
    void setID(const string& s) {
        m_id = s;
    }

    //! Returns the identifying tag for this domain.
    string id() const {
        if (m_id != "") {
            return m_id;
        } else {
            return fmt::format("domain {}", m_index);
        }
    }

    //! Print the solution.
    //! @deprecated Not implemented. To be removed after Cantera 3.1.
    virtual void show(std::ostream& s, const double* x) {
        warn_deprecated("Domain1D::show(std::ostream, double*)",
                        "Not implemented. To be removed after Cantera 3.1.");
    }

    //! Print the solution.
    //! @param x  Pointer to the local portion of the system state vector
    virtual void show(const double* x);

    //! Get the coordinate [m] of the point with local index `jlocal`
    double z(size_t jlocal) const {
        return m_z[jlocal];
    }

    //! Get the coordinate [m] of the first (leftmost) grid point in this domain
    double zmin() const {
        return m_z[0];
    }

    //! Get the coordinate [m] of the last (rightmost) grid point in this domain
    double zmax() const {
        return m_z[m_points - 1];
    }

    //! Set initial values for a component at each grid point
    //! @param name  Name of the component
    //! @param values  Array of length nPoints() containing the initial values
    //! @param soln  Pointer to the local portion of the system state vector
    void setProfile(const string& name, double* values, double* soln);

    //! Access the array of grid coordinates [m]
    vector<double>& grid() {
        return m_z;
    }

    //! Access the array of grid coordinates [m]
    const vector<double>& grid() const {
        return m_z;
    }

    //! @deprecated To be removed after Cantera 3.1. Use z() instead.
    double grid(size_t point) const {
        warn_deprecated("Domain1D::grid",
            "To be removed after Cantera 3.1. Use z() instead.");
        return m_z[point];
    }

    //! called to set up initial grid, and after grid refinement
    virtual void setupGrid(size_t n, const double* z);

    /**
     * Writes some or all initial solution values into the global solution
     * array, beginning at the location pointed to by x. This method is called
     * by the Sim1D constructor, and allows default values or ones that have
     * been set locally prior to installing this domain into the container to be
     * written to the global solution vector.
     */
    virtual void _getInitialSoln(double* x);

    //! Initial value of solution component @e n at grid point @e j.
    virtual double initialValue(size_t n, size_t j);

    /**
     * In some cases, a domain may need to set parameters that depend on the
     * initial solution estimate. In such cases, the parameters may be set in
     * method _finalize. This method is called just before the Newton solver is
     * called, and the x array is guaranteed to be the local solution vector for
     * this domain that will be used as the initial guess. If no such parameters
     * need to be set, then method _finalize does not need to be overloaded.
     */
    virtual void _finalize(const double* x) {}

    /**
     * In some cases, for computational efficiency some properties (such as
     * transport coefficients) may not be updated during Jacobian evaluations.
     * Set this to `true` to force these properties to be updated even while
     * calculating Jacobian elements.
     */
    void forceFullUpdate(bool update) {
        m_force_full_update = update;
    }

    //! Set shared data pointer
    void setData(shared_ptr<vector<double>>& data) {
        m_state = data;
    }

protected:
    //! Retrieve meta data
    virtual AnyMap getMeta() const;

    //! Retrieve meta data
    virtual void setMeta(const AnyMap& meta);

    shared_ptr<vector<double>> m_state; //!< data pointer shared from OneDim

    double m_rdt = 0.0; //!< Reciprocal of the time step
    size_t m_nv = 0; //!< Number of solution components
    size_t m_points; //!< Number of grid points
    vector<double> m_slast; //!< Solution vector at the last time step
    vector<double> m_max; //!< Upper bounds on solution components
    vector<double> m_min; //!< Lower bounds on solution components
    vector<double> m_rtol_ss; //!< Relative tolerances for steady mode
    vector<double> m_rtol_ts; //!< Relative tolerances for transient mode
    vector<double> m_atol_ss; //!< Absolute tolerances for steady mode
    vector<double> m_atol_ts; //!< Absolute tolerances for transient mode
    vector<double> m_z; //!< 1D spatial grid coordinates

    //! Parent OneDim simulation containing this and adjacent domains
    OneDim* m_container = nullptr;

    size_t m_index; //!< Left-to-right location of this domain

    //! Starting location within the solution vector for unknowns that
    //! correspond to this domain
    /*!
     * Remember there may be multiple domains associated with this problem
     */
    size_t m_iloc = 0;

    //! Index of the first point in this domain in the global point list.
    //! @see firstPoint(), lastPoint()
    size_t m_jstart = 0;

    Domain1D* m_left = nullptr; //!< Pointer to the domain to the left
    Domain1D* m_right = nullptr; //!< Pointer to the domain to the right

    //! Identity tag for the domain
    string m_id;
    unique_ptr<Refiner> m_refiner; //!< Refiner object used for placing grid points
    vector<string> m_name; //!< Names of solution components
    int m_bw = -1; //!< See bandwidth()
    bool m_force_full_update = false; //!< see forceFullUpdate()

    //! Composite thermo/kinetics/transport handler
    shared_ptr<Solution> m_solution;
};
}

#endif
