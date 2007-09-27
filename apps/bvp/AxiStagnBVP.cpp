/// @file AxiStagnBVP.cpp

#include <cantera/Cantera.h>
#include "AxiStagnBVP.h"

AxiStagnBVP::AxiStagnBVP(int nsp, int np, double L) : 
    BVP::BoundaryValueProblem(nsp+4, 
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
AxiStagnBVP::~AxiStagnBVP() {}

// specify guesses for the initial values. These can be anything
// that leads to a converged solution.
doublereal initialValue(int n, int j) { 
    switch (n) {
    case 0:
        return m_u0;
    case 1:
        return m_u0/m_L;
    case 2:
        return m_Tinf;
    case 4:
        return 1.0;
    default:
        return 0.0;
    }
}

// Specify the residual. This is where the ODE system and boundary 
// conditions are specified. The solver will attempt to find a solution
// x so that this function returns 0 for all n and j.
doublereal AxiStagnFlow::residual(doublereal* x, int n, int j) {
    
    // if n = 0, return the residual for the continuity equation
    if (n == 0) {
        if (isRight(j)) {
            return -rho_u(x,j);  // force u to zero at the right
        }
        else {
            return -(rho_u(x, j+1) - rho_u(x,j))/m_dz[j]
                -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
        }
    }

    else if (n == 1) {

        // if n = 1, then return the residual for radial momentum
        if (isLeft(j)) {
            return V(x,j);
        }
        else if (isRight(j)) {
            return V(x,j);             // force V to zero at the wall
        }
        else {
            return (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j) 
                - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
                - rdt*(V(x,j) - V_prev(j));
        }

    }

    else if (n == 2) {
        if (isLeft(j)) {
            return T(x,j) - m_Tinlet;
        }
        else if (isRight(j)) {
            return T(x,j) - m_Tsurf;
        }
        else {
            setGas(x,j);

            // heat release term
            const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
            const vector_fp& cp_R = m_thermo->cp_R_ref();

            sum = 0.0;
            sum2 = 0.0;
            doublereal flxk;
            for (k = 0; k < m_nsp; k++) {
                flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                sum += wdot(k,j)*h_RT[k];
                sum2 += flxk*cp_R[k]/m_wt[k];
            }
            sum *= GasConstant * T(x,j);
            dtdzj = dTdz(x,j);
            sum2 *= GasConstant * dtdzj;

            rsd =    - m_cp[j]*rho_u(x,j)*dtdzj 
                     - divHeatFlux(x,j) - sum - sum2;
            rsd /= (m_rho[j]*m_cp[j]);

            rsd -= rdt*(T(x,j) - T_prev(j));
        }            
    }

}


