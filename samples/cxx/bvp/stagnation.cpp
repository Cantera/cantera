/**
 * @file AxiStagnBVP.cpp
 */

// Copyright 2002  California Institute of Technology

#include <stdlib.h>
#include <time.h>

#include "AxiStagnBVP.h"
#include "cantera/base/ctml.h"
#include "cantera/oneD/MultiJac.h"

using namespace ctml;
using namespace Cantera;
using namespace std;

static void st_drawline()
{
    writelog("\n-------------------------------------"
             "------------------------------------------");
}

AxiStagnBVP::AxiStagnBVP(IdealGasPhase* ph, int nsp, int points) :
    Domain1D(nsp+4, points),
    m_inlet_u(0.0),
    m_inlet_V(0.0),
    m_inlet_T(-1.0),
    m_surface_T(-1.0),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_jac(0),
    m_ok(false),
    m_do_soret(false),
    m_transport_option(-1),
    m_efctr(0.0)
{
    m_type = cFlowType;

    m_points = points;
    m_thermo = ph;

    if (ph == 0) {
        return;    // used to create a dummy object
    }

    int nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+4, points);
    }


    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = c_offset_Y + m_nsp;

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_surfdot.resize(m_nsp, 0.0);
    m_ybar.resize(m_nsp);


    //-------------- default solution bounds --------------------

    vector_fp vmin(m_nv), vmax(m_nv);

    // no bounds on u
    vmin[0] = -1.e20;
    vmax[0] = 1.e20;

    // V
    vmin[1] = -1.e20;
    vmax[1] = 1.e20;

    // temperature bounds
    vmin[2] = 200.0;
    vmax[2]= 1.e9;

    // lamda should be negative
    vmin[3] = -1.e20;
    vmax[3] = 1.e20;

    // mass fraction bounds
    int k;
    for (k = 0; k < m_nsp; k++) {
        vmin[4+k] = -1.0e-5;
        vmax[4+k] = 1.0e5;
    }
    setBounds(vmin.size(), DATA_PTR(vmin), vmax.size(), DATA_PTR(vmax));


    //-------------------- default error tolerances ----------------
    vector_fp rtol(m_nv, 1.0e-8);
    vector_fp atol(m_nv, 1.0e-15);
    setTolerances(rtol.size(), DATA_PTR(rtol), atol.size(), DATA_PTR(atol),false);
    setTolerances(rtol.size(), DATA_PTR(rtol), atol.size(), DATA_PTR(atol),true);

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    m_refiner->setActive(3, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, DATA_PTR(gr));
    setID("stagnation flow");
}


/**
 * Change the grid size. Called after grid refinement.
 */
void AxiStagnBVP::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_enth.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);

    if (m_transport_option ==  c_Mixav_Transport) {
        m_diff.resize(m_nsp*m_points);
    } else {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_diff.resize(m_nsp*m_points);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);

    m_fixedy.resize(m_nsp, m_points);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}


void AxiStagnBVP::setupGrid(int n, const doublereal* z)
{
    resize(m_nv, n);
    int j;

    m_z[0] = z[0];
    for (j = 1; j < m_points; j++) {
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }
}


/**
 * Install a transport manager.
 */
void AxiStagnBVP::setTransport(Transport& trans, bool withSoret)
{
    m_trans = &trans;
    m_do_soret = withSoret;

    if (m_trans->model() == cMulticomponent) {
        m_transport_option = c_Multi_Transport;
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_diff.resize(m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    } else if (m_trans->model() == cMixtureAveraged) {
        m_transport_option = c_Mixav_Transport;
        m_diff.resize(m_nsp*m_points);
        if (withSoret)
            throw CanteraError("setTransport",
                               "Thermal diffusion (the Soret effect) "
                               "requires using a multicomponent transport model.");
    } else {
        throw CanteraError("setTransport","unknown transport model.");
    }
}

void AxiStagnBVP::enableSoret(bool withSoret)
{
    if (m_transport_option == c_Multi_Transport) {
        m_do_soret = withSoret;
    } else {
        throw CanteraError("setTransport",
                           "Thermal diffusion (the Soret effect) "
                           "requires using a multicomponent transport model.");
    }
}


/**
 * Set the gas object state to be consistent with the solution at
 * point j.
 */
void AxiStagnBVP::setGas(const doublereal* x,int j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}


/**
 * Set the gas state to be consistent with the solution at the
 * midpoint between j and j + 1.
 */
void AxiStagnBVP::setGasAtMidpoint(const doublereal* x,int j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(DATA_PTR(m_ybar));
    m_thermo->setPressure(m_press);
}


void AxiStagnBVP::_finalize(const doublereal* x)
{
    int k, j;
    doublereal zz, tt;
    int nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            setTemperature(j, T(x, j));
        } else {
            zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            tt = linearInterp(zz, m_zfix, m_tfix);
            setTemperature(j, tt);
        }
        for (k = 0; k < m_nsp; k++) {
            setMassFraction(j, k, Y(x, k, j));
        }
    }
    if (e) {
        solveEnergyEqn();
    }
}


//------------------------------------------------------

/**
 *  Evaluate the residual function for axisymmetric stagnation
 *  flow. If jpt is less than zero, the residual function is
 *  evaluated at all grid points. If jpt >= 0, then the residual
 *  function is only evaluated at grid points jpt-1, jpt, and
 *  jpt+1. This option is used to efficiently evaluate the
 *  Jacobian numerically.
 *
 */

void AxiStagnFlowBVP::prepare(doublereal* x)
{
    int j;

    // update thermo properties
    for (j = 0; j < m_points; j++) {
        setGas(x,j);
        m_rho[j] = m_thermo->density();
        m_wtm[j] = m_thermo->meanMolecularWeight();
        m_cp[j]  = m_thermo->cp_mass();
    }

    // update transport properties
    int k,m;

    if (m_transport_option == c_Mixav_Transport) {
        for (j = 0; j < m_points; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = m_trans->viscosity();
            m_trans->getMixDiffCoeffs(DATA_PTR(m_diff) + j*m_nsp);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    } else if (m_transport_option == c_Multi_Transport) {
        doublereal sum, sumx, wtm, dz;
        doublereal eps = 1.0e-12;
        for (m = 0; m < m_points-1; m++) {
            setGasAtMidpoint(x,m);
            dz = m_z[m+1] - m_z[m];
            wtm = m_thermo->meanMolecularWeight();

            m_visc[m] = m_trans->viscosity();
            m_trans->getMultiDiffCoeffs(m_nsp,
                                        DATA_PTR(m_multidiff) + mindex(0,0,m));

            for (k = 0; k < m_nsp; k++) {
                sum = 0.0;
                sumx = 0.0;
                for (j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum += m_wt[j]*m_multidiff[mindex(k,j,m)]*
                               ((X(x,j,m+1) - X(x,j,m))/dz + eps);
                        sumx += (X(x,j,m+1) - X(x,j,m))/dz;
                    }
                }
                m_diff[k + m*m_nsp] = sum/(wtm*(sumx+eps));
            }

            m_tcon[m] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + m*m_nsp);
            }
        }
    }
}


void AxiStagnFlowBVP::residual(doublereal* x,
                               int n, int j)
{

    int j, k;

    updateDiffFluxes(x, j0, j1);


    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    doublereal sum, sum2, dtdzj;


    //----------------------------------------------
    //         left boundary
    //----------------------------------------------

    if (j == 0) {

        // Continuity. This propagates information right-to-left,
        // since rho_u at point 0 is dependent on rho_u at point 1,
        // but not on mdot from the inlet.
        if (n == c_offset_U) {
            return -(rho_u(x,1) - rho_u(x,0))/m_dz[0]  - (density(1)*V(x,1) + density(0)*V(x,0));
        }

        // the inlet (or other) object connected to this one
        // will modify these equations by subtracting its values
        // for V, T, and mdot. As a result, these residual equations
        // will force the solution variables to the values for
        // the boundary object
        else if (n == c_offset_V) {
            return V(x,0);
        } else if (n == c_offset_T) {
            return m_Tsurf - T(x,0);
        }
        rsd[index(c_offset_T,0)] = T(x,0);
        rsd[index(c_offset_L,0)] = -rho_u(x,0);

        // The default boundary condition for species is zero
        // flux. However, the boundary object may modify
        // this.
        sum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            sum += Y(x,k,0);
            rsd[index(c_offset_Y + k, 0)] =
                -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
        }
        rsd[index(c_offset_Y, 0)] = 1.0 - sum;
    }


    //----------------------------------------------
    //
    //         right boundary
    //
    //----------------------------------------------

    else if (j == m_points - 1) {

        // the boundary object connected to the right of this
        // one may modify or replace these equations. The
        // default boundary conditions are zero u, V, and T,
        // and zero diffusive flux for all species.

        rsd[index(0,j)] = rho_u(x,j);
        rsd[index(1,j)] = V(x,j);
        rsd[index(2,j)] = T(x,j);
        rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
        diag[index(c_offset_L, j)] = 0;
        doublereal sum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            sum += Y(x,k,j);
            rsd[index(k+4,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
        }
        rsd[index(4,j)] = 1.0 - sum;
        diag[index(4,j)] = 0;

    }


    //------------------------------------------
    //     interior points
    //------------------------------------------

    else {

        //----------------------------------------------
        //    Continuity equation
        //
        //    Note that this propagates the mass flow rate
        //    information to the left (j+1 -> j) from the
        //    value specified at the right boundary. The
        //    lambda information propagates in the opposite
        //    direction.
        //
        //    d(\rho u)/dz + 2\rho V = 0
        //
        //------------------------------------------------

        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));

        //algebraic constraint
        diag[index(c_offset_U, j)] = 0;


        //------------------------------------------------
        //    Radial momentum equation
        //
        //    \rho u dV/dz + \rho V^2 = d(\mu dV/dz)/dz - lambda
        //
        //-------------------------------------------------
        rsd[index(c_offset_V,j)]
        = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
           - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
          - rdt*(V(x,j) - V_prev(j));
        diag[index(c_offset_V, j)] = 1;


        //-------------------------------------------------
        //    Species equations
        //
        //   \rho u dY_k/dz + dJ_k/dz + M_k\omega_k
        //
        //-------------------------------------------------
        getWdot(x,j);

        doublereal convec, diffus;
        for (k = 0; k < m_nsp; k++) {
            convec = rho_u(x,j)*dYdz(x,k,j);
            diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                     /(z(j+1) - z(j-1));
            rsd[index(c_offset_Y + k, j)]
            = (m_wt[k]*(wdot(k,j))
               - convec - diffus)/m_rho[j]
              - rdt*(Y(x,k,j) - Y_prev(k,j));
            diag[index(c_offset_Y + k, j)] = 1;
        }


        //-----------------------------------------------
        //    energy equation
        //-----------------------------------------------

        if (m_do_energy[j]) {

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

            rsd[index(c_offset_T, j)]   =
                - m_cp[j]*rho_u(x,j)*dtdzj
                - divHeatFlux(x,j) - sum - sum2;
            rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);

            rsd[index(c_offset_T, j)] =
                rsd[index(c_offset_T, j)] + m_efctr*(T_fixed(j) - T(x,j));

            rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
            diag[index(c_offset_T, j)] = 1;
        }

        // residual equations if the energy equation is disabled

        if (!m_do_energy[j]) {
            rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
            diag[index(c_offset_T, j)] = 0;
        }

        rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
        diag[index(c_offset_L, j)] = 0;
    }
}
}



/**
 * Update the transport properties at grid points in the range
 * from j0 to j1, based on solution x.
 */
void AxiStagnBVP::updateTransport(doublereal* x,int j0, int j1)
{
    int j,k,m;

    if (m_transport_option == c_Mixav_Transport) {
        for (j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(DATA_PTR(m_diff) + j*m_nsp);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    } else if (m_transport_option == c_Multi_Transport) {
        doublereal sum, sumx, wtm, dz;
        doublereal eps = 1.0e-12;
        for (m = j0; m < j1; m++) {
            setGasAtMidpoint(x,m);
            dz = m_z[m+1] - m_z[m];
            wtm = m_thermo->meanMolecularWeight();

            m_visc[m] = (m_dovisc ? m_trans->viscosity() : 0.0);

            m_trans->getMultiDiffCoeffs(m_nsp,
                                        DATA_PTR(m_multidiff) + mindex(0,0,m));

            for (k = 0; k < m_nsp; k++) {
                sum = 0.0;
                sumx = 0.0;
                for (j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum += m_wt[j]*m_multidiff[mindex(k,j,m)]*
                               ((X(x,j,m+1) - X(x,j,m))/dz + eps);
                        sumx += (X(x,j,m+1) - X(x,j,m))/dz;
                    }
                }
                m_diff[k + m*m_nsp] = sum/(wtm*(sumx+eps));
            }

            m_tcon[m] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + m*m_nsp);
            }
        }
    }
}


//------------------------------------------------------

/**
 *  Evaluate the residual function for axisymmetric stagnation
 *  flow. If jpt is less than zero, the residual function is
 *  evaluated at all grid points. If jpt >= 0, then the residual
 *  function is only evaluated at grid points jpt-1, jpt, and
 *  jpt+1. This option is used to efficiently evaluate the
 *  Jacobian numerically.
 *
 */

void FreeFlame::eval(int jg, doublereal* xg,
                     doublereal* rg, integer* diagg, doublereal rdt)
{

    // if evaluating a Jacobian, and the global point is outside
    // the domain of influence for this domain, then skip
    // evaluating the residual
    if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) {
        return;
    }

    // if evaluating a Jacobian, compute the steady-state residual
    if (jg >= 0) {
        rdt = 0.0;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    int jmin, jmax, jpt;
    jpt = jg - firstPoint();

    if (jg < 0) {      // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else {          // evaluate points for Jacobian
        jmin = std::max(jpt-1, 0);
        jmax = std::min(jpt+1,m_points-1);
    }

    // properties are computed for grid points from j0 to j1
    int j0 = std::max(jmin-1,0);
    int j1 = std::min(jmax+1,m_points-1);


    int j, k;


    //-----------------------------------------------------
    //              update properties
    //-----------------------------------------------------

    // update thermodynamic properties only if a Jacobian is not
    // being evaluated
    if (jpt < 0) {
        updateThermo(x, j0, j1);
        updateTransport(x, j0, j1);
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);


    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    doublereal sum, sum2, dtdzj;

    for (j = jmin; j <= jmax; j++) {


        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {

            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left,
            // since rho_u at point 0 is dependent on rho_u at point 1,
            // but not on mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one
            // will modify these equations by subtracting its values
            // for V, T, and mdot. As a result, these residual equations
            // will force the solution variables to the values for
            // the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            rsd[index(c_offset_T,0)] = T(x,0);
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero
            // flux
            sum = 0.0;
            for (k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y, 0)] = 1.0 - sum;
        }


        //----------------------------------------------
        //
        //         right boundary
        //
        //----------------------------------------------

        else if (j == m_points - 1) {

            // the boundary object connected to the right of this
            // one may modify or replace these equations. The
            // default boundary conditions are zero u, V, and T,
            // and zero diffusive flux for all species.

            // zero gradient
            rsd[index(0,j)] = rho_u(x,j) - rho_u(x,j-1);
            rsd[index(1,j)] = V(x,j);
            rsd[index(2,j)] = T(x,j) - T(x,j-1);
            doublereal sum = 0.0;
            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
            for (k = 0; k < m_nsp; k++) {
                sum += Y(x,k,j);
                rsd[index(k+4,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
            }
            rsd[index(4,j)] = 1.0 - sum;
            diag[index(4,j)] = 0;
        }

        //------------------------------------------
        //     interior points
        //------------------------------------------

        else {

            //----------------------------------------------
            //    Continuity equation
            //----------------------------------------------

            if (grid(j) > m_zfixed) {
                rsd[index(c_offset_U,j)] =
                    - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1]
                    - (density(j-1)*V(x,j-1) + density(j)*V(x,j));
            }

            else if (grid(j) == m_zfixed) {
                if (m_do_energy[j]) {
                    rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
                } else {
                    rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                                - m_rho[0]*0.3);
                }
            } else if (grid(j) < m_zfixed) {
                rsd[index(c_offset_U,j)] =
                    - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
                    - (density(j+1)*V(x,j+1) + density(j)*V(x,j));
            }
            //algebraic constraint
            diag[index(c_offset_U, j)] = 0;

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho u dV/dz + \rho V^2 = d(\mu dV/dz)/dz - lambda
            //
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;


            //-------------------------------------------------
            //    Species equations
            //
            //   \rho u dY_k/dz + dJ_k/dz + M_k\omega_k
            //
            //-------------------------------------------------
            getWdot(x,j);

            doublereal convec, diffus;
            for (k = 0; k < m_nsp; k++) {
                convec = rho_u(x,j)*dYdz(x,k,j);
                diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                         /(z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j))
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }


            //-----------------------------------------------
            //    energy equation
            //-----------------------------------------------

            if (m_do_energy[j]) {

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

                rsd[index(c_offset_T, j)]   =
                    - m_cp[j]*rho_u(x,j)*dtdzj
                    - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);

                rsd[index(c_offset_T, j)] =
                    rsd[index(c_offset_T, j)] + m_efctr*(T_fixed(j) - T(x,j));

                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                diag[index(c_offset_T, j)] = 1;
            }
            // residual equations if the energy equation is disabled
            else {
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}


/**
 * Print the solution.
 */
void AxiStagnBVP::showSolution(const doublereal* x)
{
    int nn = m_nv/5;
    int i, j, n;
    //char* buf = new char[100];
    char buf[100];

    // The mean molecular weight is needed to convert
    updateThermo(x, 0, m_points-1);

    sprintf(buf, "    Pressure:  %10.4g Pa \n", m_press);
    writelog(buf);
    for (i = 0; i < nn; i++) {
        st_drawline();
        sprintf(buf, "\n        z   ");
        writelog(buf);
        for (n = 0; n < 5; n++) {
            sprintf(buf, " %10s ",componentName(i*5 + n).c_str());
            writelog(buf);
        }
        st_drawline();
        for (j = 0; j < m_points; j++) {
            sprintf(buf, "\n %10.4g ",m_z[j]);
            writelog(buf);
            for (n = 0; n < 5; n++) {
                sprintf(buf, " %10.4g ",component(x, i*5+n,j));
                writelog(buf);
            }
        }
        writelog("\n");
    }
    int nrem = m_nv - 5*nn;
    st_drawline();
    sprintf(buf, "\n        z   ");
    writelog(buf);
    for (n = 0; n < nrem; n++) {
        sprintf(buf, " %10s ", componentName(nn*5 + n).c_str());
        writelog(buf);
    }
    st_drawline();
    for (j = 0; j < m_points; j++) {
        sprintf(buf, "\n %10.4g ",m_z[j]);
        writelog(buf);
        for (n = 0; n < nrem; n++) {
            sprintf(buf, " %10.4g ",component(x, nn*5+n,j));
            writelog(buf);
        }
    }
    writelog("\n");
}


/**
 * Update the diffusive mass fluxes.
 */
void AxiStagnBVP::updateDiffFluxes(const doublereal* x, int j0, int j1)
{
    int j, k, m;
    doublereal sum, wtm, rho, dz, gradlogT;

    switch (m_transport_option) {

    case c_Mixav_Transport:
    case c_Multi_Transport:
        for (j = j0; j < j1; j++) {
            sum = 0.0;
            wtm = m_wtm[j];
            rho = density(j);
            dz = z(j+1) - z(j);

            for (k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
            for (k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum*Y(x,k,j);
            }
        }
        break;

    default:
        throw CanteraError("updateDiffFluxes","unknown transport model");
    }

    if (m_do_soret) {
        for (m = j0; m < j1; m++) {
            gradlogT = 2.0*(T(x,m+1) - T(x,m))/(T(x,m+1) + T(x,m));
            for (k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}


string AxiStagnBVP::componentName(int n) const
{
    switch (n) {
    case 0:
        return "u";
    case 1:
        return "V";
    case 2:
        return "T";
    case 3:
        return "lambda";
    default:
        if (n >= (int) c_offset_Y && n < (int)(c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}


int AxiStagnBVP::componentIndex(string name) const
{


    if (name=="u") {
        return 0;
    } else if (name=="V") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else {
        for (int n=4; n<m_nsp+4; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
    }

    return -1;
}


void AxiStagnBVP::restore(const XML_Node& dom, doublereal* soln)
{

    vector<string> ignored;
    int nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str;
    dom.getChildren("string",str);
    int nstr = static_cast<int>(str.size());
    for (int istr = 0; istr < nstr; istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    //map<string, double> params;
    double pp = -1.0;
    pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);


    vector<XML_Node*> d;
    dom.child("grid_data").getChildren("floatArray",d);
    int nd = static_cast<int>(d.size());

    vector_fp x;
    int n, np = 0, j, ks, k;
    string nm;
    bool readgrid = false, wrote_header = false;
    for (n = 0; n < nd; n++) {
        const XML_Node& fa = *d[n];
        nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            writelog("Grid contains "+int2str(np)+
                     " points.\n");
            readgrid = true;
            setupGrid(np, DATA_PTR(x));
        }
    }
    if (!readgrid) {
        throw CanteraError("AxiStagnBVP::restore",
                           "domain contains no grid points.");
    }

    writelog("Importing datasets:\n");
    for (n = 0; n < nd; n++) {
        const XML_Node& fa = *d[n];
        nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            writelog("axial velocity   ");
            if ((int) x.size() == np) {
                for (j = 0; j < np; j++) {
                    soln[index(0,j)] = x[j];
                }
            } else {
                goto error;
            }
        } else if (nm == "z") {
            ;   // already read grid
        } else if (nm == "V") {
            writelog("radial velocity   ");
            if ((int) x.size() == np) {
                for (j = 0; j < np; j++) {
                    soln[index(1,j)] = x[j];
                }
            } else {
                goto error;
            }
        } else if (nm == "T") {
            writelog("temperature   ");
            if ((int) x.size() == np) {
                for (j = 0; j < np; j++) {
                    soln[index(2,j)] = x[j];
                }

                // For fixed-temperature simulations, use the
                // imported temperature profile by default.  If
                // this is not desired, call setFixedTempProfile
                // *after* restoring the solution.

                vector_fp zz(np);
                for (int jj = 0; jj < np; jj++) {
                    zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
                }
                setFixedTempProfile(zz, x);
            } else {
                goto error;
            }
        } else if (nm == "L") {
            writelog("lambda   ");
            if ((int) x.size() == np) {
                for (j = 0; j < np; j++) {
                    soln[index(3,j)] = x[j];
                }
            } else {
                goto error;
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            writelog(nm+"   ");
            if ((int) x.size() == np) {
                k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (j = 0; j < np; j++) {
                    soln[index(k+4,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (ignored.size() != 0) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        int nn = static_cast<int>(ignored.size());
        for (int n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    for (ks = 0; ks < nsp; ks++) {
        if (did_species[ks] == 0) {
            if (!wrote_header) {
                writelog("Missing data for species:\n");
                wrote_header = true;
            }
            writelog(m_thermo->speciesName(ks)+" ");
        }
    }

    return;
error:
    throw CanteraError("AxiStagnBVP::restore","Data size error");
}



void AxiStagnBVP::save(XML_Node& o, doublereal* sol)
{
    int k;

    Array2D soln(m_nv, m_points, sol + loc());

    XML_Node& flow = (XML_Node&)o.addChild("domain");
    flow.addAttribute("type",flowType());
    flow.addAttribute("id",m_id);
    flow.addAttribute("points",m_points);
    flow.addAttribute("components",m_nv);

    if (m_desc != "") {
        addString(flow,"description",m_desc);
    }
    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");
    addFloatArray(gv,"z",m_z.size(),DATA_PTR(m_z),
                  "m","length");
    vector_fp x(static_cast<size_t>(soln.nColumns()));

    soln.getRow(0,DATA_PTR(x));
    addFloatArray(gv,"u",x.size(),DATA_PTR(x),"m/s","velocity");

    soln.getRow(1,DATA_PTR(x));
    addFloatArray(gv,"V",
                  x.size(),DATA_PTR(x),"1/s","rate");

    soln.getRow(2,DATA_PTR(x));
    addFloatArray(gv,"T",x.size(),DATA_PTR(x),"K","temperature",0.0);

    soln.getRow(3,DATA_PTR(x));
    addFloatArray(gv,"L",x.size(),DATA_PTR(x),"N/m^4");

    for (k = 0; k < m_nsp; k++) {
        soln.getRow(4+k,DATA_PTR(x));
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),DATA_PTR(x),"","massFraction",0.0,1.0);
    }
}


void AxiStagnBVP::setJac(MultiJac* jac)
{
    m_jac = jac;
}


}  // namespace
