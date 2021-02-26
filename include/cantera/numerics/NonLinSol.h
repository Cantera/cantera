/**
 * @file NonLinSol.h
 * A nonlinear algebraic system solver, built upon Cantera's 1D multi-domain damped newton solver.
 *
 * NonLinSol is a simplified interface to the Cantera solver, designed for solving standard
 * systems of nonlinear algebraic equations. Systems are solved by Cantera as single-domain,
 * single-point problems.
 */

#include "cantera/onedim.h"

class NonLinSol : public Cantera::Domain1D
{
public:
/* TO USE THIS SOLVER:
 *     1. Include this header file in your code.
 *          #include "path/to/NonLinSol.h"
 *     2. Create a NonLinSol child class, where you'll set up your problem.
 *          class YourClass : public NonLinSol
 *     3. In the child class, provide implementations for problem-specific functions.
 *          void nonlinsol_residFunction(args) { ... }
 *          doublereal nonlinsol_initialValue(size_t i, size_t j) { ... }
 *          size_t nonlinsol_neq() { ... }
 *     4. Initialize your class and call its inherited solve() function.
 *          YourClassObject.solve();
 *     5. (Optional) Reconfigure solver settings. 
 *          YourClassObject.reconfigure(neq, lowerBound, upperBound, rtol, atol);
 */

/// IMPLEMENT THESE FUNCTIONS:

    // Specify the residual function for the system
    //  sol - iteration solution vector (input)
    //  rsd - residual vector (output)
    virtual void nonlinsol_residFunction(double *sol, double *rsd) = 0;

    // Specify guesses for the initial values.
    //  Note: called during Sim1D initialization
    virtual doublereal nonlinsol_initialValue(size_t i) = 0;

    // Number of equations (state variables) for this reactor
    virtual size_t nonlinsol_nEqs() = 0;

/// CALLABLE FUNCTIONS:

//TODO: need per component reconfigure here
    void reconfigure(int neq, double lowerBound = -1.0e-3, double upperBound = 1.01,
                     double rtol = 1.0e-4, double atol = 1.0e-9)
    {
        Domain1D::resize(neq, 1);
        for (int i = 0; i < neq; i++)
        {
            Domain1D::setBounds(i, lowerBound, upperBound);
            Domain1D::setSteadyTolerances(rtol, atol, i);
            Domain1D::setTransientTolerances(rtol, atol, i);
            Domain1D::setComponentName(i, std::to_string(i));
        }
        //TODO: per component reconfigure
        Domain1D::setBounds(0, 0, 100);
        Domain1D::setBounds(1, 0, 5);
        Domain1D::setBounds(2, -10000000, 10000000);
    }

    /**
     * Solve the nonlinear algebraic system.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel = 0)
    {
        if (Domain1D::nComponents() != nonlinsol_nEqs())
            reconfigure(nonlinsol_nEqs());

        std::vector<Cantera::Domain1D *> domains{this};
        Cantera::Sim1D(domains).solve(loglevel);
    }

private:
/// INTERNAL FUNCTIONS:

    // Implementing the residual function for the Cantera solver to use. Handles time integration, calls subclass residFunction for residuals.
    void eval(size_t jg, double *sol, double *rsd, int *timeintMask, double rdt)
    {
        nonlinsol_residFunction(sol, rsd); // call subclass residFunction() implementation to update rsd

        if (rdt == 0)
            return; // rdt is the reciprocal of the time step "dt"; rdt != 0 for time integration...

        // -------------------- TIME INTEGRATION --------------------------
        for (int i = 0; i < nonlinsol_nEqs(); i++)
        {
            rsd[i] -= rdt * (sol[i] - prevSoln(i, 0)); // backward euler method (result will be appropriately extracted from the residual)
            timeintMask[i] = 1;                        // enable time stepping for this solution component (automatically resets each iteration)
        }
    }

    // Implementing the initial value function for the Cantera solver. Grid point j is always 0 (the only point in the single-point simulation), and thus unneeded
    doublereal initialValue(size_t n, size_t j)
    {
        return nonlinsol_initialValue(n);
    }
};
