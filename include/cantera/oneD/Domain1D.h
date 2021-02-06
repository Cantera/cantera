 //! @file Domain1D.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DOMAIN1D_H
#define CT_DOMAIN1D_H

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include "refine.h"

namespace Cantera
{

// domain types
const int cFlowType = 50;
const int cFreeFlow = 51;
const int cAxisymmetricStagnationFlow = 52;
const int cConnectorType = 100;
const int cSurfType = 102;
const int cInletType = 104;
const int cSymmType = 105;
const int cOutletType = 106;
const int cEmptyType = 107;
const int cOutletResType = 108;
const int cPorousType = 109;

class MultiJac;
class OneDim;
class XML_Node;

/**
 * Base class for one-dimensional domains.
 * @ingroup onedim
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

    virtual ~Domain1D() {}
    Domain1D(const Domain1D&) = delete;
    Domain1D& operator=(const Domain1D&) = delete;

    //! Domain type flag.
    int domainType() {
        return m_type;
    }

    //! The left-to-right location of this domain.
    size_t domainIndex() {
        return m_index;
    }

    //! True if the domain is a connector domain.
    bool isConnector() {
        return (m_type >= cConnectorType);
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

    /*!
     * Initialize. This method is called by OneDim::init() for each domain once
     * at the beginning of a simulation. Base class method does nothing, but may
     * be overloaded.
     */
    virtual void init() {  }

    virtual void setInitialState(doublereal* xlocal = 0) {}
    virtual void setState(size_t point, const doublereal* state, doublereal* x) {}

    /*!
     * When called, this function should reset "bad" values in the state vector
     * such as negative species concentrations. This function may be called
     * after a failed solution attempt.
     */
    virtual void resetBadValues(double* xg) {}

    /*!
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

    //! Name of the nth component. May be overloaded.
    virtual std::string componentName(size_t n) const;

    void setComponentName(size_t n, const std::string& name) {
        m_name[n] = name;
    }

    //! index of component with name \a name.
    virtual size_t componentIndex(const std::string& name) const;

    void setBounds(size_t n, doublereal lower, doublereal upper) {
        m_min[n] = lower;
        m_max[n] = upper;
    }

    //! Set tolerances for time-stepping mode
    /*!
     * @param rtol Relative tolerance
     * @param atol Absolute tolerance
     * @param n    component index these tolerances apply to. If set to -1 (the
     *      default), these tolerances will be applied to all solution
     *      components.
     */
    void setTransientTolerances(doublereal rtol, doublereal atol, size_t n=npos);

    //! Set tolerances for steady-state mode
    /*!
     * @param rtol Relative tolerance
     * @param atol Absolute tolerance
     * @param n    component index these tolerances apply to. If set to -1 (the
     *     default), these tolerances will be applied to all solution
     *     components.
     */
    void setSteadyTolerances(doublereal rtol, doublereal atol, size_t n=npos);

    //! Relative tolerance of the nth component.
    doublereal rtol(size_t n) {
        return (m_rdt == 0.0 ? m_rtol_ss[n] : m_rtol_ts[n]);
    }

    //! Absolute tolerance of the nth component.
    doublereal atol(size_t n) {
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
    doublereal upperBound(size_t n) const {
        return m_max[n];
    }

    //! Lower bound on the nth component
    doublereal lowerBound(size_t n) const {
        return m_min[n];
    }

    //! Prepare to do time stepping with time step dt
    /*!
     * Copy the internally-stored solution at the last time step to array x0.
     */
    void initTimeInteg(doublereal dt, const doublereal* x0) {
        std::copy(x0 + loc(), x0 + loc() + size(), m_slast.begin());
        m_rdt = 1.0/dt;
    }

    //! Prepare to solve the steady-state problem
    /*!
     * Set the internally-stored reciprocal of the time step to 0.0
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

    /*!
     * Set this if something has changed in the governing
     * equations (e.g. the value of a constant has been changed,
     * so that the last-computed Jacobian is no longer valid.
     */
    void needJacUpdate();

    //! Evaluate the residual function at point j. If j == npos,
    //! evaluate the residual function at all points.
    /*!
     *  This function must be implemented in classes derived from Domain1D.
     *
     *  @param j  Grid point at which to update the residual
     *  @param[in] x  State vector
     *  @param[out] r  residual vector
     *  @param[out] mask  Boolean mask indicating whether each solution
     *      component has a time derivative (1) or not (0).
     *  @param[in] rdt Reciprocal of the timestep (`rdt=0` implies steady-
     *  state.)
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt=0.0) {
        throw NotImplementedError("Domain1D::eval");
    }

    size_t index(size_t n, size_t j) const {
        return m_nv*j + n;
    }
    doublereal value(const doublereal* x, size_t n, size_t j) const {
        return x[index(n,j)];
    }

    virtual void setJac(MultiJac* jac) {}

    //! Save the current solution for this domain into an XML_Node
    /*!
     * Base class version of the general domain1D save function. Derived classes
     * should call the base class method in addition to saving their own data.
     *
     * @param o    XML_Node to save the solution to.
     * @param sol  Current value of the solution vector. The object will pick
     *             out which part of the solution vector pertains to this
     *             object.
     * @return     XML_Node created to represent this domain
     *
     * @deprecated The XML output format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    //! Restore the solution for this domain from an XML_Node
    /*!
     * Base class version of the general Domain1D restore function. Derived
     * classes should call the base class method in addition to restoring
     * their own data.
     *
     * @param dom XML_Node for this domain
     * @param soln Current value of the solution vector, local to this object.
     * @param loglevel 0 to suppress all output; 1 to show warnings; 2 for
     *      verbose output
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

    size_t size() const {
        return m_nv*m_points;
    }

    /**
     * Find the index of the first grid point in this domain, and
     * the start of its variables in the global solution vector.
     */
    void locate();

    /**
     * Location of the start of the local solution vector in the global
     * solution vector,
     */
    virtual size_t loc(size_t j = 0) const {
        return m_iloc;
    }

    /**
     * The index of the first (i.e., left-most) grid point belonging to this
     * domain.
     */
    size_t firstPoint() const {
        return m_jstart;
    }

    /**
     * The index of the last (i.e., right-most) grid point belonging to this
     * domain.
     */
    size_t lastPoint() const {
        return m_jstart + m_points - 1;
    }

    /**
     * Set the left neighbor to domain 'left.' Method 'locate' is called to
     * update the global positions of this domain and all those to its right.
     */
    void linkLeft(Domain1D* left) {
        m_left = left;
        locate();
    }

    //! Set the right neighbor to domain 'right.'
    void linkRight(Domain1D* right) {
        m_right = right;
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
    void setID(const std::string& s) {
        m_id = s;
    }

    std::string id() const {
        if (m_id != "") {
            return m_id;
        } else {
            return fmt::format("domain {}", m_index);
        }
    }

    virtual void showSolution_s(std::ostream& s, const doublereal* x) {}

    //! Print the solution.
    virtual void showSolution(const doublereal* x);

    doublereal z(size_t jlocal) const {
        return m_z[jlocal];
    }
    doublereal zmin() const {
        return m_z[0];
    }
    doublereal zmax() const {
        return m_z[m_points - 1];
    }

    void setProfile(const std::string& name, double* values, double* soln);

    vector_fp& grid() {
        return m_z;
    }
    const vector_fp& grid() const {
        return m_z;
    }
    doublereal grid(size_t point) const {
        return m_z[point];
    }

    //! called to set up initial grid, and after grid refinement
    virtual void setupGrid(size_t n, const doublereal* z);

    /**
     * Writes some or all initial solution values into the global solution
     * array, beginning at the location pointed to by x. This method is called
     * by the Sim1D constructor, and allows default values or ones that have
     * been set locally prior to installing this domain into the container to be
     * written to the global solution vector.
     */
    virtual void _getInitialSoln(doublereal* x);

    //! Initial value of solution component \a n at grid point \a j.
    virtual doublereal initialValue(size_t n, size_t j);

    /**
     * In some cases, a domain may need to set parameters that depend on the
     * initial solution estimate. In such cases, the parameters may be set in
     * method _finalize. This method is called just before the Newton solver is
     * called, and the x array is guaranteed to be the local solution vector for
     * this domain that will be used as the initial guess. If no such parameters
     * need to be set, then method _finalize does not need to be overloaded.
     */
    virtual void _finalize(const doublereal* x) {}

    /**
     * In some cases, for computational efficiency some properties (e.g.
     * transport coefficients) may not be updated during Jacobian evaluations.
     * Set this to `true` to force these properties to be udpated even while
     * calculating Jacobian elements.
     */
    void forceFullUpdate(bool update) {
        m_force_full_update = update;
    }

protected:
    doublereal m_rdt;
    size_t m_nv;
    size_t m_points;
    vector_fp m_slast;
    vector_fp m_max;
    vector_fp m_min;
    vector_fp m_rtol_ss, m_rtol_ts;
    vector_fp m_atol_ss, m_atol_ts;
    vector_fp m_z;
    OneDim* m_container;
    size_t m_index;
    int m_type;

    //! Starting location within the solution vector for unknowns that
    //! correspond to this domain
    /*!
     * Remember there may be multiple domains associated with this problem
     */
    size_t m_iloc;

    size_t m_jstart;

    Domain1D* m_left, *m_right;

    //! Identity tag for the domain
    std::string m_id;
    std::string m_desc;
    std::unique_ptr<Refiner> m_refiner;
    std::vector<std::string> m_name;
    int m_bw;
    bool m_force_full_update;
};
}

#endif
