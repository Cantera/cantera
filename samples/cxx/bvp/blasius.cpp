/// @file blasius.cpp
/// The Blasius boundary layer

#include "BoundaryValueProblem.h"

using Cantera::npos;

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

    // Specify the residual function. This is where the ODE system and boundary
    // conditions are specified. The solver will attempt to find a solution
    // x so that rsd is zero.
    void eval(size_t jg, double* x, double* rsd, int* diag, double rdt) {
        size_t jpt = jg - firstPoint();
        size_t jmin, jmax;
        if (jg == npos) {  // evaluate all points
            jmin = 0;
            jmax = m_points - 1;
        } else {  // evaluate points for Jacobian
            jmin = std::max<size_t>(jpt, 1) - 1;
            jmax = std::min(jpt+1,m_points-1);
        }

        for (size_t j = jmin; j <= jmax; j++) {
            if (j == 0) {
                rsd[index(0,j)] = zeta(x,j);
                rsd[index(1,j)] = u(x,j);
            } else if (j == m_points - 1) {
                rsd[index(0,j)] = leftFirstDeriv(x,0,j) - u(x,j);
                rsd[index(1,j)] = u(x,j) - 1.0;
            } else {
                rsd[index(0,j)] = leftFirstDeriv(x,0,j) - u(x,j);
                rsd[index(1,j)] = cdif2(x,1,j) + 0.5*zeta(x,j)*centralFirstDeriv(x,1,j)
                                  - rdt*(value(x,1,j) - prevSoln(1,j));
                diag[index(1,j)] = 1;
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
        // Solve the equations, refining the grid as needed
        eqs.solve(1);
        // write the solution to a CSV file.
        eqs.writeCSV("blasius.csv");
        return 0;
    } catch (Cantera::CanteraError& err) {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
