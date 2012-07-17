/// @file AxiStagnBVP.h

#include "cantera/Cantera.h"
#include "BoundaryValueProblem.h"


/**
 * This class solves
 */
class AxiStagnBVP : public BVP::BoundaryValueProblem
{

public:

    AxiStagnBVP(int nsp, int np, double L) : BVP::BoundaryValueProblem(nsp+4,
                np, 0.0, L) {

        // specify the component bounds, error tolerances, and names.
        BVP::Component u;
        u.lower = -200.0;
        u.upper = 200.0;
        u.rtol = 1.0e-8;
        u.atol = 1.0e-15;
        u.name = "u";
        setComponent(0, u);  // the axial velocity will be component 0

        BVP::Component V;
        V.lower = -1.0e8;
        V.upper = 1.0e8;
        V.rtol = 1.0e-8;
        V.atol = 1.0e-15;
        V.name = "V";
        setComponent(1, V);  // the radial velocity will be component 1

        BVP::Component T;
        T.lower = 200.0;
        T.upper = 1.0e9;
        T.rtol = 1.0e-8;
        T.atol = 1.0e-15;
        T.name = "T";
        setComponent(2, T);  // the temperature will be component 2

        BVP::Component lambda;
        lambda.lower = -1.0e20;
        lambda.upper = 1.0e20;
        lambda.rtol = 1.0e-8;
        lambda.atol = 1.0e-15;
        lambda.name = "Lambda";
        setComponent(3, lambda);  // the pressure-gradient eigenvalue will be
        //component 3
        BVP::Component Y;
        Y.lower = -1.0e-5;
        Y.upper = 1.0e2;
        Y.rtol = 1.0e-8;
        Y.atol = 1.0e-15;
        for (k = 0; k < nsp; k++) {
            Y.name = thermo->speciesName(k);
            setComponent(k+4, Y);
        }
    }


    // destructor
    virtual ~AxiStagnBVP() {}


    // specify guesses for the initial values. These can be anything
    // that leads to a converged solution.
    doublereal AxiStagnBVP::initialValue(int n, int j) {
        switch (n) {
        case 0:
            return m;
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
            } else
                // this implements d(zeta)/dz = u
            {
                return (zeta(x,j) - zeta(x,j-1))/(z(j)-z(j-1)) - u(x,j);
            }
        }
        // if n = 1, then return the residual for the second ODE
        else {
            if (isLeft(j)) { // here we specify u(0) = 0
                return u(x,j);
            } else if (isRight(j)) { // and here we specify u(L) = 1
                return u(x,j) - 1.0;
            } else
                // this implements the 2nd ODE
            {
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
        AxiStagnBVP eqs(6, 10.0);
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


