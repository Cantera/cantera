/// @file BoundaryValueProblem.h
/// Simplified interface to the capabilities provided by Cantera to
/// solve boundary value problems.

#ifndef BVP_H
#define BVP_H

#include "cantera/onedim.h"
#include <fstream>

/// Namespace for the boundary value problem package.
namespace BVP
{

// default grid refinement parameters
const double max_grid_ratio = 4.0; ///< max ratio of neighboring grid intervals
const double max_delta = 0.01; ///< max difference in function values
const double max_delta_slope = 0.02; ///< max difference in slopes
const double prune = 0.000; ///< don't remove grid points

/**
 * Used to specify component-specific options for method
 * setComponent of method BoundaryValueProblem. An instance of
 * class Component should be created for each solution component,
 * and its values set appropriately.
 */
class Component
{
public:
    double lower; ///< lower bound
    double upper; ///< upper bound
    double rtol; ///< relative error tolerance
    double atol; ///< absolute error tolerance
    bool refine; ///< make this component active for grid refinement
    std::string name; ///< component name

    /**
     * Constructor. Sets default values.
     */
    Component() : lower(0.0), upper(1.0), rtol(1.0e-9), atol(1.0e-12),
        refine(true) {}
};

/**
 * Base class for boundary value problems. This class is designed
 * to provide a simplified interface to the capabilities Cantera
 * provides to solve boundary value problems. Classes for specific
 * boundary value problems should be derived from this one.
 *
 * Class BoundaryValueProblem derives from Cantera's Domain1D
 * class.
 */
class BoundaryValueProblem : public Cantera::Domain1D
{

public:

    /**
     * Constructor. This constructor begins with a uniform grid of
     * np points starting at zmin, and ending at zmax.
     *
     * @param nv Number of solution components
     * @param np Number of grid points in initial grid
     * @param zmin Location of left-hand side of domain
     * @param zmax Location of right-hand side of domain
     */
    BoundaryValueProblem(int nv, int np,
                         doublereal zmin, doublereal zmax) :
        m_left(0), m_right(0), m_sim(0)
    {
        // Create the initial uniform grid
        Cantera::vector_fp z(np);
        int iz;
        for (iz = 0; iz < np; iz++) {
            z[iz] = zmin + iz*(zmax - zmin)/(np-1);
        }
        setupGrid(np, z.data());
        resize(nv, np);
    }

    /**
     * Constructor. This alternate constructor starts with a
     * specified grid, unlike the above that uses a uniform grid
     * to start. The array z must contain the z coordinates of np
     * grid points.
     */
    BoundaryValueProblem(int nv, int np,
                         doublereal* z) :
        m_left(0), m_right(0), m_sim(0)
    {
        setupGrid(np, z);
        resize(nv, np);
    }

    /**
     * Destructor. Deletes the dummy terminator domains, and the
     * solver.
     */
    virtual ~BoundaryValueProblem() {
        delete m_left;
        delete m_right;
        delete m_sim;
    }

    /**
     *  Set parameters and options for solution component \a n.
     *  This method should be invoked for each solution component
     *  before calling 'solve'. The parameter values should first
     *  be set by creating an instance of class Component, and
     *  setting its member data appropriately.
     *
     *  @param n Component number.
     *  @param c Component parameter values
     */
    void setComponent(size_t n, Component& c) {
        if (m_sim == 0) {
            start();
        }
        if (n >= m_nv) {
            throw Cantera::CanteraError("BoundaryValueProblem::setComponent",
                                        "Illegal solution component number");
        }
        // set the upper and lower bounds for this component
        setBounds(n, c.lower, c.upper);
        // set the error tolerances
        setSteadyTolerances(c.rtol, c.atol, n);
        setTransientTolerances(c.rtol, c.atol, n);
        // specify whether this component should be considered in
        // refining the grid
        m_refiner->setActive(n, c.refine);
        // set a default name if one has not been entered
        if (c.name == "") {
            c.name = fmt::format("Component {}", n);
        }
        setComponentName(n, c.name);
    }

    /**
     * Solve the boundary value problem.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel=0) {
        if (m_sim == 0) {
            start();
        }
        bool refine = true;
        m_sim->solve(loglevel, refine);
    }

    /**
     * Write the solution to a CSV file.
     * @param filename CSV file name.
     * @param ztitle Title for 'z' column.
     * @param dotitles If true, begin with a row of column titles.
     */
    void writeCSV(std::string filename = "output.csv",
                  bool dotitles = true, std::string ztitle = "z") const {
        std::ofstream f(filename);
        int np = nPoints();
        int nc = nComponents();
        int n, m;
        if (dotitles) {
            f << ztitle << ", ";
            for (m = 0; m < nc; m++) {
                f << componentName(m);
                if (m != nc - 1) {
                    f << ", ";
                }
            }
            f << std::endl;
        }
        for (n = 0; n < np; n++) {
            f << z(n) << ", ";
            for (m = 0; m < nc; m++) {
                f << m_sim->value(1, m, n);
                if (m != nc - 1) {
                    f << ", ";
                }
            }
            f << std::endl;
        }
    }

    /**
     * Initial value of solution component \a n at initial grid
     * point \a j. The default is zero for all components at all
     * grid points. Overload in derived classes to specify other
     * choices for initial values.
     */
    virtual doublereal initialValue(size_t n, size_t j) {
        return 0.0;
    }

protected:
    Cantera::Domain1D* m_left; ///< dummy terminator
    Cantera::Domain1D* m_right; ///< dummy terminator
    Cantera::Sim1D* m_sim; ///< controller for solution

    /**
     * Set up the problem. Creates the solver instance, and sets
     * default grid refinement parameters. This method is called
     * internally, and does not need to be invoked explicitly in
     * derived classes.
     */
    void start() {
        // Add dummy terminator domains on either side of this one.
        m_left = new Cantera::Empty1D;
        m_right = new Cantera::Empty1D;
        std::vector<Cantera::Domain1D*> domains { m_left, this, m_right };

        // create the Sim1D instance that will control the
        // solution process
        m_sim = new Cantera::Sim1D(domains);

        // set default grid refinement parameters
        m_sim->setRefineCriteria(1, max_grid_ratio, max_delta,
                                 max_delta_slope, prune);
    }

    /**
     * @name Trial Solution Derivatives
     * These methods return
     * derivatives of individual components at specified grid
     * points, using a given trial solution.  They are designed
     * for use in writing overloaded versions of method 'residual'
     * in derived classes.
     */

    //@{

    /**
     * This method is provided for use in method residual when
     * central-differenced second derivatives are needed.
     * @param x The current trial solution vector.
     * @param n Component index.
     * @param j Grid point number.
     */
    doublereal cdif2(const doublereal* x, int n, int j) const {
        doublereal c1 = value(x,n,j) - value(x,n,j-1);
        doublereal c2 = value(x,n,j+1) - value(x,n,j);
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/
               (z(j+1) - z(j-1));
    }

    /**
     * The first derivative of solution component n at point j.
     * If type is -1, the first derivative is computed using the
     * value to the left of point j, if it is +1 then the
     * value to the right is used, and if it is zero (default) a
     * central-differenced first derivative is computed.
    */
    doublereal firstDeriv(const doublereal* x, int n, int j,
                          int type = 0) const {
        switch (type) {
        case -1:
            return leftFirstDeriv(x, n, j);
        case 1:
            return rightFirstDeriv(x, n, j);
        default:
            return centralFirstDeriv(x, n, j);
        }
    }

    /**
     * First derivative of component \a n at point \a j. The derivative
     * is formed to the right of point j, using values at point j
     * and point j + 1.
     */
    doublereal rightFirstDeriv(const doublereal* x, int n, int j) const {
        return (value(x,n,j+1) - value(x,n,j))/(z(j+1) - z(j));
    }

    /**
     * First derivative of component \a n at point \a j. The derivative
     * is formed to the left of point j, using values at point j
     * and point j - 1.
     */

    doublereal leftFirstDeriv(const doublereal* x, int n, int j) const {
        return (value(x,n,j) - value(x,n,j-1))/(z(j) - z(j-1));
    }

    /**
     * This method is provided for use in method residual when
     * central-differenced first derivatives are needed.
     * @param x The current trial solution vector.
     * @param n Component index.
     * @param j Grid point number.
     */
    doublereal centralFirstDeriv(const doublereal* x, int n, int j) const {
        doublereal c1 = value(x,n,j+1) - value(x,n,j-1);
        return c1/(z(j+1) - z(j-1));
    }

    //@}
};
}
#endif
