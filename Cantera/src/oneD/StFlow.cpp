/**
 * @file StFlow.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <stdlib.h>
#include <time.h>

#include "StFlow.h"
#include "../ArrayViewer.h"
#include "ctml.h"
#include "MultiJac.h"

using namespace ctml;

namespace Cantera {


    //-------------------  importSolution ------------------------

    /**
     * Import a previous solution to use as an initial estimate. The
     * previous solution may have been computed using a different
     * reaction mechanism. Species in the old and new mechanisms are
     * matched by name, and any species in the new mechanism that were
     * not in the old one are set to zero. The new solution is created
     * with the same number of grid points as in the old solution.
     */ 
    void importSolution(int points,  
        doublereal* oldSoln, igthermo_t& oldmech,
        int size_new, doublereal* newSoln, igthermo_t& newmech) {
        
        // Number of components in old and new solutions
        int nv_old = oldmech.nSpecies() + 4;
        int nv_new = newmech.nSpecies() + 4;

        if (size_new < nv_new*points) {
            throw CanteraError("importSolution",
                "new solution array must have length "+
                int2str(nv_new*points));
        }

        int n, j, knew;
        string nm;

        // copy u,V,T,lambda
        for (j = 0; j < points; j++) 
            for (n = 0; n < 4; n++) 
                newSoln[nv_new*j + n] = oldSoln[nv_old*j + n];

        // copy mass fractions        
        int nsp0 = oldmech.nSpecies();
        int nsp1 = newmech.nSpecies();

        // loop over the species in the old mechanism
        for (int k = 0; k < nsp0; k++) {
            nm = oldmech.speciesName(k);      // name

            // location of this species in the new mechanism.
            // If < 0, then the species is not in the new mechanism.
            knew = newmech.speciesIndex(nm); 

            // copy this species from the old to the new solution vectors
            if (knew >= 0) {
                for (j = 0; j < points; j++) {
                    newSoln[nv_new*j + 4 + knew] = oldSoln[nv_old*j + 4 + k];
                }
            }
        }


        // normalize mass fractions
        for (j = 0; j < points; j++) {
            newmech.setMassFractions(&newSoln[nv_new*j + 4]);
            newmech.getMassFractions(nsp1,&newSoln[nv_new*j + 4]);
        }
    }


    //---------------------- drawline ----------------------------------

    static void drawline(ostream& s) {
        s << "\n-------------------------------------"
          <<  "------------------------------------------";
    }

    static void drawline() {
        writelog("\n-------------------------------------"
            "------------------------------------------");
    }


    //--------------------- linear interp ------------------------------

    StFlow::StFlow(igthermo_t* ph, int nsp, int points) : 
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

        m_boundary.resize(2,0);

        m_points = points;
        m_thermo = ph;

        if (ph == 0) return; // used to create a dummy object

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
        setBounds(vmin.size(), vmin.begin(), vmax.size(), vmax.begin());


        //-------------------- default error tolerances ----------------
        vector_fp rtol(m_nv, 1.0e-8);
        vector_fp atol(m_nv, 1.0e-15);
        setTolerances(rtol.size(), rtol.begin(), atol.size(), atol.begin(),false);
        setTolerances(rtol.size(), rtol.begin(), atol.size(), atol.begin(),true);

        //-------------------- grid refinement -------------------------
        m_refiner->setActive(0, false);
        m_refiner->setActive(1, false);
        m_refiner->setActive(2, false);
        m_refiner->setActive(3, false);

        vector_fp gr;
        for (int ng = 0; ng < m_points; ng++) gr.push_back(1.0*ng/m_points);
        setupGrid(m_points, gr.begin());
        setID("stagnation flow");
    }


    /**
     * Change the grid size. Called after grid refinement.
     */
    void StFlow::resize(int points) {
        Domain1D::resize(m_nv, points);

        m_rho.resize(m_points, 0.0);
        m_wtm.resize(m_points, 0.0);
        m_cp.resize(m_points, 0.0);
        m_enth.resize(m_points, 0.0);
        m_visc.resize(m_points, 0.0);
        m_tcon.resize(m_points, 0.0);

        if (m_transport_option ==  c_Mixav_Transport) {
            m_diff.resize(m_nsp*m_points);
        }
        else {
            m_diff.resize(m_nsp*m_nsp*m_points);
        }
        m_flux.resize(m_nsp,m_points);
        m_wdot.resize(m_nsp,m_points, 0.0);
        m_do_energy.resize(m_points,false);

        m_fixedy.resize(m_nsp, m_points);
        m_fixedtemp.resize(m_points);        

        m_dz.resize(m_points-1);
        m_z.resize(m_points);
    }        
        


    void StFlow::setupGrid(int n, const doublereal* z) {        
        resize(n);
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
    void StFlow::setTransport(Transport& trans, bool withSoret) {
        m_trans = &trans;
        m_do_soret = withSoret;

        if (m_trans->model() == cMulticomponent) {
            m_transport_option = c_Multi_Transport;
        }
        else if (m_trans->model() == cMixtureAveraged) {
            m_transport_option = c_Mixav_Transport;
            if (withSoret) 
                throw CanteraError("setTransport",
                    "Thermal diffusion (the Soret effect) "
                    "requires using a multicomponent transport model.");
        }
        else 
            throw CanteraError("setTransport","unknown transport model.");
    }


    /**
     * Set the gas object state to be consistent with the solution at
     * point j.
     */
    void StFlow::setGas(const doublereal* x,int j) {
        m_thermo->setTemperature(T(x,j));
        const doublereal* yy = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions_NoNorm(yy);
        m_thermo->setPressure(m_press);
    }


    /**
     * Set the gas state to be consistent with the solution at the
     * midpoint between j and j + 1.
     */
    void StFlow::setGasAtMidpoint(const doublereal* x,int j) {
        m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
        const doublereal* yyj = x + m_nv*j + c_offset_Y;
        const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
        for (int k = 0; k < m_nsp; k++)
            m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
        m_thermo->setMassFractions_NoNorm(m_ybar.begin());
        m_thermo->setPressure(m_press);
    }



    /**
     *  Evaluate the residual function for axisymmetric stagnation
     *  flow. If jpt is less than zero, the residual function is
     *  evaluated at all grid points. If jpt >= 0, then the residual
     *  function is only evaluated at grid points jpt-1, jpt, and
     *  jpt+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     *
     */

    void AxiStagnFlow::eval(int jg, doublereal* xg, 
        doublereal* rg, integer* diagg, doublereal rdt) {

        // if evaluating a Jacobian, and the global point is outside
        // the domain of influence for this domain, then skip
        // evaluating the residual
        if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) return;

        // if evaluating a Jacobian, compute the steady-state residual
        if (jg >= 0) rdt = 0.0;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* rsd = rg + loc();
        integer* diag = diagg + loc();
        
        int jmin, jmax, jpt;
        jpt = jg - firstPoint();

        if (jg < 0) {      // evaluate all points
            jmin = 0;
            jmax = m_points - 1;
        }
        else {            // evaluate points for Jacobian
            jmin = max(jpt-1, 0);
            jmax = min(jpt+1,m_points-1);
        }

        // properties are computed for grid points from j0 to j1
        int j0 = max(jmin-1,0);
        int j1 = min(jmax+1,m_points-1);


        int j, k;


        //----------------------------------------------------- 
        //              update properties
        //-----------------------------------------------------

        // thermodynamic properties only if a Jacobian is
        // not being evaluated
        if (jpt < 0) 
            updateThermo(x, j0, j1);

        // update transport properties only if a Jacobian is
        // not being evaluated
        if (jpt < 0) 
            updateTransport(x, j0, j1);

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
                 // flux. However, the boundary object may modify
                 // this.
                 for (k = 0; k < m_nsp; k++) {
                     rsd[index(c_offset_Y + k, 0)] =  
                         -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
                 }
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
                    //if (m_do_species[k]) {
                    convec = rho_u(x,j)*dYdz(x,k,j);
                    diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                             /(z(j+1) - z(j-1));
                    rsd[index(c_offset_Y + k, j)]   
                        = (m_wt[k]*(wdot(k,j) ) 
                            - convec - diffus)/m_rho[j]
                        - rdt*(Y(x,k,j) - Y_prev(k,j));
                    diag[index(c_offset_Y + k, j)] = 1;
                    //}
                }


                //-----------------------------------------------
                //    energy equation
                //-----------------------------------------------

                if (m_do_energy[j]) {

                    setGas(x,j);

                    // heat release term
                    const vector_fp& h_RT = m_thermo->enthalpy_RT();
                    const vector_fp& cp_R = m_thermo->cp_R();
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
            }

            // residual equations if the energy or species equations
            // are disabled

            for (k = 0; k < m_nsp; k++) {
                if (!m_do_species[k]) {
                    rsd[index(c_offset_Y+k,j)] = Y(x,k,j) - Y_fixed(k,j);
                    diag[index(c_offset_Y+k, j)] = 0;
                }
            }

            if (!m_do_energy[j]) {
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            // Propagate lambda from left to right
            if (j > 0) {
                rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
                diag[index(c_offset_L, j)] = 0;
            }
        }
    }



    /**
     * Update the transport properties at grid points in the range
     * from j0 to j1, based on solution x.
     */
    void AxiStagnFlow::updateTransport(doublereal* x,int j0, int j1) {
        int j;
        if (m_transport_option == c_Mixav_Transport) {
            for (j = j0; j < j1; j++) {
                setGasAtMidpoint(x,j);
                m_visc[j] = m_trans->viscosity();
                m_trans->getMixDiffCoeffs(m_diff.begin() + j*m_nsp);
                m_tcon[j] = m_trans->thermalConductivity();
            }
        }
        else if (m_transport_option == c_Multi_Transport) {
            for (j = j0; j < j1; j++) {
                setGasAtMidpoint(x,j);
                m_visc[j] = m_trans->viscosity();
                m_trans->getMultiDiffCoeffs(m_nsp, 
                    m_diff.begin() + mindex(0,0,j));
                m_tcon[j] = m_trans->thermalConductivity();
            }
        }
    }


    void OneDFlow::eval(int jg, doublereal* xg, doublereal* rg, integer* diagg,
        doublereal rdt) {

        static double elapsed;
        //        doublereal rtau = 1.e5;

        clock_t t0 = clock();

        //        doublereal rdt_save = rdt;
        if (jg >= 0) rdt = 0.0;

        if (jg >= 0 && (jg < firstPoint() || jg > lastPoint())) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* rsd = rg + loc();
        integer* diag = diagg + loc();

        int jmin, jmax, jpt;
        jpt = jg - firstPoint();

        for (int jj = 0; jj < m_points*m_nv; jj++) {
            if (x[jj] < -1.e20 || x[jj] > 1.e20) {
                showSolution(cout, x);
                throw CanteraError("tlt","tlt");
            }
        }

        // the residual function is evaluated for jmin <= j <= jmax, and
        // properties and evaluated for j0 <= j <= j1. 

        if (jg < 0) {
            jmin = 0; 
            jmax = m_points - 1;
        }
        else {                   
            jmin = max(jpt-1,0);
            jmax = min(jpt+1,m_points-1);
        }
        int j0 = max(jmin-1,0);
        int j1 = min(jmax+1,m_points-1);

        int j, k;

        //----------------------------------------------------- 
        // compute properties needed in the residual equations
        //-----------------------------------------------------

        // for each point, synchronize the state of the fluid object
        // with the current solution values, and then use this object
        // to compute the density, mean molecular weight, and mean
        // specific heat at constant pressure.
        if (jpt < 0) updateThermo(x, j0, j1);

        // skip updating transport properties if a Jacobian is
        // being evaluated
        if (jpt < 0) updateTransport(x, j0, j1);

        // update the species diffusive mass fluxes
        updateDiffFluxes(x, j0, j1);


        //----------------------------------------------------
        // evaluate the residual equations at all required
        // grid points
        //----------------------------------------------------

        doublereal sum, sum2, deltaz, dtdzj;


        for (j = jmin; j <= jmax; j++) {


            //----------------------------------------------
            //         boundaries
            //----------------------------------------------

            if (j == 0) {
                setGas(x,0);
                m_boundary[0]->eval(x, m_rho[0], m_flux.begin(), 
                    rsd);
            }

            else if (j == m_points - 1) {
                m_boundary[1]->eval(x + index(0, j), m_rho[j], 
                    m_flux.begin() + m_nsp*(j-1), 
                    rsd + index(0, j));
            }


            //------------------------------------------
            //     interior points                 
            //------------------------------------------
            
            else {
                
                //    continuity
                rsd[index(c_offset_U,j)] = (rho_u(x,j-1) - rho_u(x,j));
                
                //    radial velocity = 0
                rsd[index(c_offset_V,j)] = V(x,j);
                
                //    species equations
                getWdot(x,j);

                doublereal convec, diffus;
                for (k = 0; k < m_nsp; k++) {
                    if (m_do_species[k]) {

                        convec = rho_u(x,j) * dYdz(x,k,j);               
                        diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))/(z(j+1) - z(j-1));
                        rsd[index(c_offset_Y + k, j)] = 
                            (m_wt[k]*wdot(k,j) - convec - diffus)/m_rho[j]
                            - rdt*(Y(x,k,j) - Y_prev(k,j));
                        diag[index(c_offset_Y + k, j)] = 1;
                    }
                    else 
                        rsd[index(c_offset_Y+k,j)] = (Y(x,k,j) - Y_fixed(k,j));
                }


                //    energy equation

                if (m_do_energy[j]) {
                    setGas(x,j);

                    // heat release term
                    const vector_fp& h_RT = m_thermo->enthalpy_RT();
                    const vector_fp& cp_R = m_thermo->cp_R();
                    sum = 0.0;
                    sum2 = 0.0;
                    deltaz = (z(j+1) - z(j-1));
                    doublereal flxk;
                    for (k = 0; k < m_nsp; k++) {
                        flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));

                        sum += wdot(k,j)*h_RT[k];
                        sum2 += flxk*cp_R[k]/m_wt[k];
                    }
                    sum *= GasConstant * T(x,j);
                    dtdzj = (T(x,j+1) - T(x,j-1))/deltaz; // dTdz(x,j) + (m_dz[j-1]/deltaz)*(dTdz(x,j+1) - dTdz(x,j));
                    sum2 *= GasConstant * dtdzj;
                    rsd[index(c_offset_T, j)]   = - m_cp[j]*rho_u(x,j)*dtdzj 
                                           - divHeatFlux(x,j) - sum - sum2;
                    rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                    rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                    diag[index(c_offset_T, j)] = 1;
                }
                
                //    lambda = 0
                rsd[index(c_offset_L, j)] = lambda(x,j);

            }
            for (k = 0; k < m_nsp; k++) {
                if (!m_do_species[k]) {
                    rsd[index(c_offset_Y+k,j)] = 
                        (Y(x,k,j) - Y_fixed(k,j));
                    diag[index(c_offset_Y+k, j)] = 0;
                }
            }
            if (!m_do_energy[j]) {
                rsd[index(c_offset_T, j)] = (T(x,j) - T_fixed(j));
                diag[index(c_offset_T, j)] = 0;
            }

        }
        clock_t t1 = clock();
        elapsed += double(t1 - t0)/CLOCKS_PER_SEC;
    }


    /**
     * Update the transport properties at grid points in the range
     * from j0 to j1, based on solution x.
     */
    void OneDFlow::updateTransport(doublereal* x,int j0, int j1) {
        int j;
        if (m_transport_option == c_Mixav_Transport) {
            for (j = j0; j < j1; j++) {
                setGasAtMidpoint(x,j);
                m_trans->getMixDiffCoeffs(m_diff.begin() + j*m_nsp);
                m_tcon[j] = m_trans->thermalConductivity();
            }
        }
        else if (m_transport_option == c_Multi_Transport) {
            for (j = j0; j < j1; j++) {
                setGasAtMidpoint(x,j);
                m_trans->getMultiDiffCoeffs(m_nsp, m_diff.begin() + mindex(0,0,j));
                m_tcon[j] = m_trans->thermalConductivity();
            }
        }
    }


    /**
     * Print the solution.
     */    
    void StFlow::showSolution(ostream& s, const doublereal* x) {
        int nn = m_nv/5;
        int i, j, n;
        char* buf = new char[100];

        // The mean molecular weight is needed to convert
        updateThermo(x, 0, m_points-1);

        for (i = 0; i < nn; i++) {
            drawline(s);
            sprintf(buf, "\n        z   ");
            s << buf;
            for (n = 0; n < 5; n++) { 
                sprintf(buf, " %10s ",componentName(i*5 + n).c_str());
                s << buf;
            }
            drawline(s);
            for (j = 0; j < m_points; j++) {
                sprintf(buf, "\n %10.4g ",m_z[j]);
                s << buf;
                for (n = 0; n < 5; n++) { 
                    sprintf(buf, " %10.4g ",component(x, i*5+n,j));
                    s << buf;
                }
            }
            s << endl;
        }
        int nrem = m_nv - 5*nn;
        drawline(s);
        sprintf(buf, "\n        z   ");
        s << buf;
        for (n = 0; n < nrem; n++) {
            sprintf(buf, "  %10s  ", componentName(nn*5 + n).c_str());
            s << buf;
        }
        drawline(s);
        for (j = 0; j < m_points; j++) {
            sprintf(buf, "\n  %10.4g ",m_z[j]);
            s << buf;
            for (n = 0; n < nrem; n++) { 
                sprintf(buf, "  %10.4g  ",component(x, nn*5+n,j));
                s << buf;
            }
        }
        s << endl;
    }


    /**
     * Print the solution.
     */    
    void StFlow::showSolution(const doublereal* x) {
        int nn = m_nv/5;
        int i, j, n;
        char* buf = new char[100];

        // The mean molecular weight is needed to convert
        updateThermo(x, 0, m_points-1);

        sprintf(buf, "    Pressure:  %10.4g Pa \n", m_press);
        writelog(buf);
        for (i = 0; i < nn; i++) {
            drawline();
            sprintf(buf, "\n        z   ");
            writelog(buf);
            for (n = 0; n < 5; n++) { 
                sprintf(buf, " %10s ",componentName(i*5 + n).c_str());
                writelog(buf);
            }
            drawline();
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
        drawline();
        sprintf(buf, "\n        z   ");
        writelog(buf);
        for (n = 0; n < nrem; n++) {
            sprintf(buf, " %10s ", componentName(nn*5 + n).c_str());
            writelog(buf);
        }
        drawline();
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
    void StFlow::updateDiffFluxes(const doublereal* x, int j0, int j1) {
        int j, k, m;
        doublereal sum, wtm, rho, dz, gradlogT, s;

        switch (m_transport_option) {

        case c_Mixav_Transport:
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
                for (k = 0; k < m_nsp; k++) m_flux(k,j) += sum*Y(x,k,j);
            } 
            break;

        case c_Multi_Transport:
            for (m = j0; m < j1; m++) {
                wtm = m_wtm[m];
                rho = density(m);
                dz = z(m+1) - z(m);
                
                for (k = 0; k < m_nsp; k++) {
                    sum = 0.0;
                    for (j = 0; j < m_nsp; j++) {
                        s = m_wt[j]*m_diff[mindex(k,j,m)];
                        s *= (X(x,k,m+1) - X(x,k,m))/dz;
                        sum += s;
                    }
                    m_flux(k,m) = sum*rho*m_wt[k]/(wtm*wtm);
                }
            } 
            break;
        default:
            throw CanteraError("updateDiffFluxes","unknown transport model");
        }

        if (m_do_soret) {
            throw CanteraError("updateDiffFluxes","not yet");
            for (m = j0; m < j1; m++) {
                gradlogT = 2.0*(T(x,m+1) - T(x,m))/(T(x,m+1) + T(x,m));
                for (k = 0; k < m_nsp; k++) {
                    m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
                }
            }
        }
    }



    void StFlow::outputTEC(ostream &s, const doublereal* x, 
        string title, int zone) {
        int j,k;
        s << "TITLE     = \"" + title + "\"" << endl;
        s << "VARIABLES = \"Z (m)\"" << endl;
        s << "\"u (m/s)\"" << endl;
        s << "\"V (1/s)\"" << endl;
        s << "\"T (K)\"" << endl;
        s << "\"lambda\"" << endl;

        for (k = 0; k < m_nsp; k++) {
            s << "\"" << m_thermo->speciesName(k) << "\"" << endl;
        }
        s << "ZONE T=\"c" << zone << "\"" << endl;
        s << " I=" << m_points << ",J=1,K=1,F=POINT" << endl;
        s << "DT=(SINGLE SINGLE SINGLE SINGLE";
        for (k = 0; k < m_nsp; k++) s << " SINGLE";
        s << " )" << endl;
        for (j = 0; j < m_points; j++) {
            s << z(j) << " ";
            for (k = 0; k < m_nv; k++) {
                s << component(x, k, j) << " ";
            }
            s << endl;
        }
    }


    string StFlow::componentName(int n) const {
        switch(n) {
        case 0: return "u";
        case 1: return "V";
        case 2: return "T";
        case 3: return "lambda";
        default:
            if (n >= (int) c_offset_Y && n < (int) (c_offset_Y + m_nsp)) {
                    return m_thermo->speciesName(n - c_offset_Y);
            }
                 //    if (m_do_species[n - c_offset_Y]) 
                 //   return m_thermo->speciesName(n - c_offset_Y)+"  ";
                 // else
                 //   return m_thermo->speciesName(n - c_offset_Y)+" *";
                 //}
            else 
                return "<unknown>";
        }
    }


    /**
     * Returns true if all necessary parameters have been set; otherwise it 
     * throws an exception.
     */
    bool StFlow::ready() {
        if (m_press < 0.0) { 
            throw CanteraError("StFlow::ready",
                "pressure not specified - call setPressure");
            return false;
        }
        if (m_points == 0) {
            throw CanteraError("StFlow::ready",
                "grid not specified - call setupGrid");
            return false;
        }
        if (m_nsp < 0) {
            throw CanteraError("StFlow::ready",
                "fluid not specified - call specifyFluid");
            return false;
        }
        if (m_boundary[0] == 0 || m_boundary[1] == 0) {
            throw CanteraError("StFlow::ready",
                "boundaries not specified - call setBoundary");
            return false;
        }
        m_ok = true;
        return m_ok;
    }

    
    void StFlow::restore(int job, 
        string fname, string id, int& size_z, doublereal* z, 
        int& size_soln, doublereal* soln) {

        vector<string> ignored;
        int nsp = m_thermo->nSpecies();
        vector_int did_species(nsp, 0);

        ifstream s(fname.c_str());
        if (!s) 
            throw CanteraError("StFlow::restore",
                "could not open input file "+fname);

        XML_Node root;
        root.build(s);
        s.close();
        int k;

        const XML_Node* f = root.findID(id);
        if (!f) {
            throw CanteraError("StFlow::restore","No solution with id = "+id);
        }
        
        const XML_Node& flow = f->child("flowfield");
        f = &flow;

        //if (f->name() != "flowfield") {
        //    throw CanteraError("StFlow::restore","The element with id "
        //        +id+" does not contain flowfield data.");
        //}

        vector<XML_Node*> str;
        f->getChildren("string",str);
        int nstr = str.size();
        for (int istr = 0; istr < nstr; istr++) {
            const XML_Node& nd = *str[istr];
            writelog(nd["title"]+": "+nd.value()+"\n");
        }

        vector<XML_Node*> d;
        f->child("grid_data").getChildren("floatArray",d);
        int nd = d.size();

        vector_fp x;
        int n, np = 0, j, ks;
        string nm;
        bool readgrid = false, wrote_header = false;
        for (n = 0; n < nd; n++) {
            const XML_Node& fa = *d[n];
            nm = fa["title"];
            if (nm == "z") {
                getFloatArray(fa,x,false);
                np = x.size();
                if (job == -1) {
                    size_z = np;
                    //size_soln = (nd - 1)*np;
                    size_soln = (m_nsp + 4)*np;
                    return;
                }
                writelog("Grid contains "+int2str(np)+
                    " points.\n");
                if (size_z < np) {
                    throw CanteraError("restore",
                        "grid array must be have length at least "
                        +int2str(np));
                } 
//                 if (size_soln < (m_nsp + 4)*np) {
//                     throw CanteraError("restore",
//                         "solution array must have length at least "
//                         +int2str((m_nsp + 4)*np));
//                 } 
                copy(x.begin(), x.end(), z);
                readgrid = true;
            }
        }
        if (!readgrid) {
            throw CanteraError("StFlow::restore",
                "solution contains no grid points.");
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
                }
                else {
                    goto error;
                }
            }
            else if (nm == "z") {
                ;
            }
            else if (nm == "V") {
                writelog("radial velocity   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(1,j)] = x[j];
                }
                else goto error;
            }
            else if (nm == "T") {
                writelog("temperature   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(2,j)] = x[j];
                }
                else goto error;
            }
            else if (nm == "L") {
                writelog("lambda   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(3,j)] = x[j];
                }
                else goto error;
            }
            else if (m_thermo->speciesIndex(nm) >= 0) {
                writelog(nm+"   ");
                if ((int) x.size() == np) {
                    k = m_thermo->speciesIndex(nm);
                    did_species[k] = 1;
                    for (j = 0; j < np; j++) 
                        soln[index(k+4,j)] = x[j];
                }
            }
            else
                ignored.push_back(nm);
        }

        if (ignored.size() != 0) {
            writelog("\n\n");
            writelog("Ignoring datasets:\n");
            int nn = ignored.size();
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

        writelog("\n\nFinished importing solution.\n\n");
        return;
 error:
        throw CanteraError("StFlow::restore","Data size error");
    }



    void StFlow::restore(const XML_Node& dom, doublereal* soln) {

        vector<string> ignored;
        int nsp = m_thermo->nSpecies();
        vector_int did_species(nsp, 0);

        vector<XML_Node*> str;
        dom.getChildren("string",str);
        int nstr = str.size();
        for (int istr = 0; istr < nstr; istr++) {
            const XML_Node& nd = *str[istr];
            writelog(nd["title"]+": "+nd.value()+"\n");
        }

        map<string, double> params;
        getFloats(dom, params);
        setPressure(params["pressure"]);


        vector<XML_Node*> d;
        dom.child("grid_data").getChildren("floatArray",d);
        int nd = d.size();

        vector_fp x;
        int n, np, j, ks, k;
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

                // note that setupGrid also resizes the domain.
                setupGrid(np, x.begin());
            }
        }
        if (!readgrid) {
            throw CanteraError("StFlow::restore",
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
                }
                else {
                    goto error;
                }
            }
            else if (nm == "z") {
                ;   // already read grid
            }
            else if (nm == "V") {
                writelog("radial velocity   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(1,j)] = x[j];
                }
                else goto error;
            }
            else if (nm == "T") {
                writelog("temperature   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(2,j)] = x[j];

                    // For fixed-temperature simulations, use the imported temperature profile by default. 
                    // If this is not desired, call setFixedTempProfile *after* restoring the solution.
                    vector_fp zz(np);
                    for (int jj = 0; jj < np; jj++) zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
                    setFixedTempProfile(zz, x);
                }
                else goto error;
            }
            else if (nm == "L") {
                writelog("lambda   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(3,j)] = x[j];
                }
                else goto error;
            }
            else if (m_thermo->speciesIndex(nm) >= 0) {
                writelog(nm+"   ");
                if ((int) x.size() == np) {
                    k = m_thermo->speciesIndex(nm);
                    did_species[k] = 1;
                    for (j = 0; j < np; j++) 
                        soln[index(k+4,j)] = x[j];
                }
            }
            else
                ignored.push_back(nm);
        }

        if (ignored.size() != 0) {
            writelog("\n\n");
            writelog("Ignoring datasets:\n");
            int nn = ignored.size();
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
        throw CanteraError("StFlow::restore","Data size error");
    }



    void StFlow::save(XML_Node& o, doublereal* sol) {
        int k;

        ArrayViewer soln(m_nv, m_points, sol + loc());

        XML_Node& flow = (XML_Node&)o.addChild("domain");
        flow.addAttribute("type",flowType());
        flow.addAttribute("id",m_id);
        flow.addAttribute("points",m_points);
        flow.addAttribute("components",m_nv);

        if (m_desc != "") addString(flow,"description",m_desc);
        XML_Node& gv = flow.addChild("grid_data");
        addFloat(flow, "pressure", m_press, "Pa", "pressure"); 
        addFloatArray(gv,"z",m_z.size(),m_z.begin(),
            "m","length");
        vector_fp x(soln.nColumns());

        soln.getRow(0,x.begin());
        addFloatArray(gv,"u",x.size(),x.begin(),"m/s","velocity");

        soln.getRow(1,x.begin());
        addFloatArray(gv,"V",
            x.size(),x.begin(),"1/s","rate");

        soln.getRow(2,x.begin());
        addFloatArray(gv,"T",x.size(),x.begin(),"K","temperature",0.0);

        soln.getRow(3,x.begin());
        addFloatArray(gv,"L",x.size(),x.begin(),"N/m^4");

        for (k = 0; k < m_nsp; k++) {
            soln.getRow(4+k,x.begin());
            addFloatArray(gv,m_thermo->speciesName(k),
                x.size(),x.begin(),"","massFraction",0.0,1.0);
        }
    }


    void StFlow::setJac(MultiJac* jac) {
        m_jac = jac;
    }

    //void StFlow::requestJacUpdate() {
    //    if (m_jac) m_jac->setAge(10000);
    //}

    void StFlow::setEnergyFactor(doublereal efctr) {
        doublereal de = efctr - m_efctr;
        m_efctr = efctr;
        int strt = loc();
        int jg;
        for (int j = 1; j < m_points - 1; j++) {
            jg = strt + index(c_offset_T, j);
            m_jac->incrementDiagonal(jg, -de);
        }
    }


}
