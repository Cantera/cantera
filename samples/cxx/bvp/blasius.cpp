/// @file blasius.cpp
/// The Blasius boundary layer

#include "BoundaryValueProblem.h"

/**
 * This class solves the Blasius boundary value problem on the domain (0,L):
 * \f[
 *             \frac{d\zeta}{dz} = u.
 * \f]
 * \f[
 *             \frac{d^2u}{dz^2} + 0.5\zeta \frac{du}{dz} = 0.
 * \f]
 * with boundary conditions
 * \f[
 * \zeta(0) = 0, u(0) = 0, u(L) = 1.
 * \f]
 * Note that this is formulated as a system of two equations, with maximum
 * order of 2, rather than as a single third-order boundary value problem.
 * For reasons having to do with the band structure of the Jacobian, no
 * equation in the system should have order greater than 2.
 */
class Blasius : public BVP::BoundaryValueProblem
{
public:
    // This problem has two components (zeta and u)
    Blasius(int np, double L) : BVP::BoundaryValueProblem(2, np, 0.0, L) {
        // specify the component bounds, error tolerances, and names.
        BVP::Component A;
        A.lower = -200.0;
        A.upper = 200.0;
        A.rtol = 1.0e-12;
        A.atol = 1.0e-15;
        A.name = "zeta";
        setComponent(0, A); // zeta will be component 0

        BVP::Component B;
        B.lower = -200.0;
        B.upper = 200.0;
        B.rtol = 1.0e-12;
        B.atol = 1.0e-15;
        B.name = "u";
        setComponent(1, B); // u will be component 1
    }

    // destructor
    virtual ~Blasius() {}

    // specify guesses for the initial values. These can be anything
    // that leads to a converged solution.
    virtual doublereal initialValue(size_t n, size_t j) {
        switch (n) {
        case 0:
            return 0.1*z(j);
        case 1:
            return 0.5*z(j);
        default:
            return 0.0;
        }
    }

    // Specify the residual. This is where the ODE system and boundary
    // conditions are specified. The solver will attempt to find a solution
    // x so that this function returns 0 for all n and j.
    virtual doublereal residual(doublereal* x, size_t n, size_t j) {
        // if n = 0, return the residual for the first ODE
        if (n == 0) {
            if (isLeft(j)) { // here we specify zeta(0) = 0
                return zeta(x,j);
            } else {
                // this implements d(zeta)/dz = u
                return (zeta(x,j) - zeta(x,j-1))/(z(j)-z(j-1)) - u(x,j);
            }
        } else {
            // if n = 1, then return the residual for the second ODE
            if (isLeft(j)) { // here we specify u(0) = 0
                return u(x,j);
            } else if (isRight(j)) { // and here we specify u(L) = 1
                return u(x,j) - 1.0;
            } else {
                // this implements the 2nd ODE
                return cdif2(x,1,j) + 0.5*zeta(x,j)*centralFirstDeriv(x,1,j);
            }
        }
    }

private:
    // for convenience only. Note that the compiler will inline these.
    double zeta(double* x, int j) {
        return value(x,0,j);
    }
    double u(double* x, int j) {
        return value(x,1,j);
    }
};

int main()
{
    try {
        // Specify a problem on (0,10), with an initial uniform grid of
        // 6 points.
        Blasius eqs(6, 10.0);
        // Solve the equations, refining the grid as needed, and print lots of diagnostic output (loglevel = 4)
        eqs.solve(4);
        // write the solution to a CSV file.
        eqs.writeCSV();
        return 0;
    } catch (Cantera::CanteraError& err) {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
