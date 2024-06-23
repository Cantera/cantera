//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/refine.h"
#include "cantera/transport/Transport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_nsp(nsp)
{
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }
    m_thermo = ph;

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // set pressure based on associated thermo object
    setPressure(m_thermo->pressure());

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
    m_hk.resize(m_nsp, m_points, 0.0);
    m_dhk_dz.resize(m_nsp, m_points - 1, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(c_offset_U, -1e20, 1e20); // no bounds on u
    setBounds(c_offset_V, -1e20, 1e20); // V
    setBounds(c_offset_T, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    setBounds(c_offset_L, -1e20, 1e20); // lambda should be negative
    setBounds(c_offset_E, -1e20, 1e20); // no bounds for inactive component

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_L, false);

    vector<double> gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");
}

StFlow::StFlow(shared_ptr<ThermoPhase> th, size_t nsp, size_t points)
    : StFlow(th.get(), nsp, points)
{
    m_solution = Solution::create();
    m_solution->setThermo(th);
}

StFlow::StFlow(shared_ptr<Solution> sol, const string& id, size_t points)
    : StFlow(sol->thermo().get(), sol->thermo()->nSpecies(), points)
{
    m_solution = sol;
    m_id = id;
    m_kin = m_solution->kinetics().get();
    m_trans = m_solution->transport().get();
    
    m_solution->registerChangedCallback(this, [this]() {
        setKinetics(m_solution->kinetics());
        setTransport(m_solution->transport());
    });
    
    if (m_trans->transportModel() == "none") {
        // @deprecated
        warn_deprecated("StFlow",
            "An appropriate transport model\nshould be set when instantiating the "
            "Solution ('gas') object.\nImplicit setting of the transport model "
            "is deprecated and\nwill be removed after Cantera 3.0.");
        setTransportModel("mixture-averaged");
    }

}

StFlow::~StFlow()
{
    if (m_solution) {
        m_solution->removeChangedCallback(this);
    }
}

string StFlow::type() const {
    if (m_isFree) {
        return "free-flow";
    }
    if (m_usesLambda) {
        return "axisymmetric-flow";
    }
    return "unstrained-flow";
}

void StFlow::setThermo(ThermoPhase& th) {
    warn_deprecated("StFlow::setThermo", "To be removed after Cantera 3.0.");
    m_thermo = &th;
}

void StFlow::setKinetics(shared_ptr<Kinetics> kin)
{
    if (!m_solution) {
        // @todo remove after Cantera 3.0
        throw CanteraError("StFlow::setKinetics",
            "Unable to update object that was not constructed from smart pointers.");
    }
    m_kin = kin.get();
    m_solution->setKinetics(kin);
}

void StFlow::setKinetics(Kinetics& kin)
{
    warn_deprecated("StFlow::setKinetics(Kinetics&)", "To be removed after Cantera 3.0."
        " Replaced by setKinetics(shared_ptr<Kinetics>).");
    m_kin = &kin;
}

void StFlow::setTransport(shared_ptr<Transport> trans)
{
    if (!m_solution) {
        // @todo remove after Cantera 3.0
        throw CanteraError("StFlow::setTransport",
            "Unable to update object that was not constructed from smart pointers.");
    }
    if (!trans) {
        throw CanteraError("StFlow::setTransport", "Unable to set empty transport.");
    }
    m_trans = trans.get();
    if (m_trans->transportModel() == "none") {
        throw CanteraError("StFlow::setTransport", "Invalid Transport model 'none'.");
    }
    m_do_multicomponent = (m_trans->transportModel() == "multicomponent" ||
        m_trans->transportModel() == "multicomponent-CK");

    m_diff.resize(m_nsp * m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp * m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_solution->setTransport(trans);
}

void StFlow::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);
    m_enthp.resize(m_nsp,m_points,0.0);

    Eneg_RHS_term_1.resize(m_points, 0.0);
    Eneg_RHS_term_2.resize(m_points, 0.0);
    Eneg_RHS_term_3.resize(m_points, 0.0);
    Eneg_RHS_term_4.resize(m_points, 0.0);
    
    Sp_RHS_term_1.resize(m_points, m_nsp, 0.0);
    Sp_RHS_term_2.resize(m_points, m_nsp, 0.0);
    Sp_RHS_term_3.resize(m_points, m_nsp, 0.0);
    
    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_hk.resize(m_nsp, m_points, 0.0);
    m_dhk_dz.resize(m_nsp, m_points - 1, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}

void StFlow::setupGrid(size_t n, const double* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("StFlow::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }
}

void StFlow::resetBadValues(double* xg)
{
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void StFlow::setTransportModel(const string& trans)
{
    if (!m_solution) {
        // @todo remove after Cantera 3.0
        throw CanteraError("StFlow::setTransportModel",
            "Unable to set Transport manager by name as object was not initialized\n"
            "from a Solution manager: set Transport object directly instead.");
    }
    m_solution->setTransportModel(trans);
}

string StFlow::transportModel() const {
    return m_trans->transportModel();
}

void StFlow::setTransport(Transport& trans)
{
    warn_deprecated("StFlow::setTransport(Transport&)", "To be removed after"
        " Cantera 3.0. Replaced by setTransport(shared_ptr<Transport>).");
    m_trans = &trans;
    if (m_trans->transportModel() == "none") {
        throw CanteraError("StFlow::setTransport",
            "Invalid Transport model 'none'.");
    }
    m_do_multicomponent = (m_trans->transportModel() == "multicomponent" ||
                           m_trans->transportModel() == "multicomponent-CK");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void StFlow::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
        m_rho[j] = m_thermo->density();
    }
}

void StFlow::setGas(const double* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const double* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const double* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const double* yyj = x + m_nv*j + c_offset_Y;
    const double* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

bool StFlow::fixed_mdot() {
    warn_deprecated("StFlow::fixed_mdot", "To be removed after"
        " Cantera 3.0. Replaced by isFree().");
    return !m_isFree;
}

void StFlow::_finalize(const double* x)
{
    if (!m_do_multicomponent && m_do_soret) {
        throw CanteraError("StFlow::_finalize",
            "Thermal diffusion (the Soret effect) is enabled, and requires "
            "using a multicomponent transport model.");
    }

    size_t nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (size_t j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            m_fixedtemp[j] = T(x, j);
        } else {
            double zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            double tt = linearInterp(zz, m_zfix, m_tfix);
            m_fixedtemp[j] = tt;
        }
    }
    if (e) {
        solveEnergyEqn();
    }

    if (m_isFree) {
        // If the domain contains the temperature fixed point, make sure that it
        // is correctly set. This may be necessary when the grid has been modified
        // externally.
        if (m_tfixed != Undef) {
            for (size_t j = 0; j < m_points; j++) {
                if (z(j) == m_zfixed) {
                    return; // fixed point is already set correctly
                }
            }

            for (size_t j = 0; j < m_points - 1; j++) {
                // Find where the temperature profile crosses the current
                // fixed temperature.
                if ((T(x, j) - m_tfixed) * (T(x, j+1) - m_tfixed) <= 0.0) {
                    m_tfixed = T(x, j+1);
                    m_zfixed = z(j+1);
                    return;
                }
            }
        }
    }
}

void StFlow::eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jg, x, jmin, jmax);
    evalResidual(x, rsd, diag, rdt, jmin, jmax);
}

void StFlow::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
{
    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
    }
    if (jg == npos) {
        double* Yleft = x + index(c_offset_Y, jmin);
        m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
        double* Yright = x + index(c_offset_Y, jmax);
        m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);
}

void StFlow::evalResidual(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation (see docstring of enableRadiation)
    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        double k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const double c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const double c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            rsd[index(c_offset_T,0)] = T(x,0);
            if (m_usesLambda) {
                rsd[index(c_offset_L, 0)] = -rho_u(x, 0);
            } else {
                rsd[index(c_offset_L, 0)] = lambda(x, 0);
                diag[index(c_offset_L, 0)] = 0;
            }

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            if (m_usesLambda) {
                rsd[index(c_offset_V,j)] =
                    (shear(x, j) - lambda(x, j) - rho_u(x, j) * dVdz(x, j)
                    - m_rho[j] * V(x, j) * V(x, j)) / m_rho[j]
                    - rdt * (V(x, j) - V_prev(j));
                diag[index(c_offset_V, j)] = 1;
            } else {
                rsd[index(c_offset_V, j)] = V(x, j);
                diag[index(c_offset_V, j)] = 0;
            }

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j))
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
                
                //output terms in specis equations                
                Sp_RHS_term_1(j,k)=-convec;
                Sp_RHS_term_2(j,k)=-diffus;
                Sp_RHS_term_3(j,k)=m_wt[k]*(wdot(k,j));                
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {

                setGas(x,j);
                double dtdzj = dTdz(x,j);
                double sum = 0.0;
                double sum1=0;//used for heat release term
                double sum2=0;//used for heat convection from species diffusion

                grad_hk(x, j);
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*m_hk(k,j);
                    sum += flxk * m_dhk_dz(k,j) / m_wt[k];
                    sum1 += wdot(k,j)*m_hk(k,j);
                    sum2 += flxk * m_dhk_dz(k,j) / m_wt[k];;
                }
                
                Eneg_RHS_term_1[j]=m_cp[j]*rho_u(x,j)*dtdzj;//heat convection from species convection
                Eneg_RHS_term_2[j]=divHeatFlux(x,j);//heat conduction
                Eneg_RHS_term_3[j]=sum1;//heat generation
                Eneg_RHS_term_4[j]=sum2;//heat convection from species diffusion                

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            if (m_usesLambda) {
                rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j - 1);
            } else {
                rsd[index(c_offset_L, j)] = lambda(x, j);
            }
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::updateTransport(double* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            double wtm = m_thermo->meanMolecularWeight();
            double rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void StFlow::show(const double* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::show(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

void StFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                double sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            double wtm = m_wtm[j];
            double rho = density(j);
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum*Y(x,k,j);
            }
        }
    }

    if (m_do_soret) {
        for (size_t m = j0; m < j1; m++) {
            double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                              ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}

string StFlow::componentName(size_t n) const
{
    switch (n) {
    case c_offset_U:
        return "velocity";
    case c_offset_V:
        return "spread_rate";
    case c_offset_T:
        return "T";
    case c_offset_L:
        return "lambda";
    case c_offset_E:
        return "eField";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const string& name) const
{
    if (name=="velocity") {
        return c_offset_U;
    } else if (name=="spread_rate") {
        return c_offset_V;
    } else if (name=="T") {
        return c_offset_T;
    } else if (name=="lambda") {
        return c_offset_L;
    } else if (name == "eField") {
        return c_offset_E;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
        throw CanteraError("StFlow1D::componentIndex",
                           "no component named " + name);
    }
}

bool StFlow::componentActive(size_t n) const
{
    switch (n) {
    case c_offset_V: // spread_rate
        return m_usesLambda;
    case c_offset_L: // lambda
        return m_usesLambda;
    case c_offset_E: // eField
        return false;
    default:
        return true;
    }
}

AnyMap StFlow::getMeta() const
{
    AnyMap state = Domain1D::getMeta();
    state["transport-model"] = m_trans->transportModel();

    state["phase"]["name"] = m_thermo->name();
    AnyValue source = m_thermo->input().getMetadata("filename");
    state["phase"]["source"] = source.empty() ? "<unknown>" : source.asString();

    state["radiation-enabled"] = m_do_radiation;
    if (m_do_radiation) {
        state["emissivity-left"] = m_epsilon_left;
        state["emissivity-right"] = m_epsilon_right;
    }

    set<bool> energy_flags(m_do_energy.begin(), m_do_energy.end());
    if (energy_flags.size() == 1) {
        state["energy-enabled"] = m_do_energy[0];
    } else {
        state["energy-enabled"] = m_do_energy;
    }

    state["Soret-enabled"] = m_do_soret;

    set<bool> species_flags(m_do_species.begin(), m_do_species.end());
    if (species_flags.size() == 1) {
        state["species-enabled"] = m_do_species[0];
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            state["species-enabled"][m_thermo->speciesName(k)] = m_do_species[k];
        }
    }

    state["refine-criteria"]["ratio"] = m_refiner->maxRatio();
    state["refine-criteria"]["slope"] = m_refiner->maxDelta();
    state["refine-criteria"]["curve"] = m_refiner->maxSlope();
    state["refine-criteria"]["prune"] = m_refiner->prune();
    state["refine-criteria"]["grid-min"] = m_refiner->gridMin();
    state["refine-criteria"]["max-points"] =
        static_cast<long int>(m_refiner->maxPoints());

    if (m_zfixed != Undef) {
        state["fixed-point"]["location"] = m_zfixed;
        state["fixed-point"]["temperature"] = m_tfixed;
    }

    return state;
}

shared_ptr<SolutionArray> StFlow::asArray(const double* soln) const
{
    auto arr = SolutionArray::create(
        m_solution, static_cast<int>(nPoints()), getMeta());
    arr->addExtra("grid", false); // leading entry
    AnyValue value;
    value = m_z;
    arr->setComponent("grid", value);
    vector<double> data(nPoints());
    for (size_t i = 0; i < nComponents(); i++) {
        if (componentActive(i)) {
            auto name = componentName(i);
            for (size_t j = 0; j < nPoints(); j++) {
                data[j] = soln[index(i, j)];
            }
            if (!arr->hasComponent(name)) {
                arr->addExtra(name, componentIndex(name) > c_offset_Y);
            }
            value = data;
            arr->setComponent(name, value);
        }
    }
    value = m_rho;
    arr->setComponent("D", value); // use density rather than pressure

    if (m_do_radiation) {
        arr->addExtra("radiative-heat-loss", true); // add at end
        value = m_qdotRadiation;
        arr->setComponent("radiative-heat-loss", value);
    }

    return arr;
}

void StFlow::fromArray(SolutionArray& arr, double* soln)
{
    Domain1D::setMeta(arr.meta());
    arr.setLoc(0);
    auto phase = arr.thermo();
    m_press = phase->pressure();

    const auto grid = arr.getComponent("grid").as<vector<double>>();
    setupGrid(nPoints(), &grid[0]);

    for (size_t i = 0; i < nComponents(); i++) {
        if (!componentActive(i)) {
            continue;
        }
        string name = componentName(i);
        if (arr.hasComponent(name)) {
            const vector<double> data = arr.getComponent(name).as<vector<double>>();
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(i,j)] = data[j];
            }
        } else {
            warn_user("StFlow::fromArray", "Saved state does not contain values for "
                "component '{}' in domain '{}'.", name, id());
        }
    }

    updateProperties(npos, soln + loc(), 0, m_points - 1);
    setMeta(arr.meta());
}

string StFlow::flowType() const {
    warn_deprecated("StFlow::flowType",
        "To be removed after Cantera 3.0; superseded by 'type'.");
    if (m_type == cFreeFlow) {
        return "Free Flame";
    } else if (m_type == cAxisymmetricStagnationFlow) {
        return "Axisymmetric Stagnation";
    } else {
        throw CanteraError("StFlow::flowType", "Unknown value for 'm_type'");
    }
}

void StFlow::setMeta(const AnyMap& state)
{
    if (state.hasKey("energy-enabled")) {
        const AnyValue& ee = state["energy-enabled"];
        if (ee.isScalar()) {
            m_do_energy.assign(nPoints(), ee.asBool());
        } else {
            m_do_energy = ee.asVector<bool>(nPoints());
        }
    }

    setTransportModel(state.getString("transport-model", "mixture-averaged"));

    if (state.hasKey("Soret-enabled")) {
        m_do_soret = state["Soret-enabled"].asBool();
    }

    if (state.hasKey("species-enabled")) {
        const AnyValue& se = state["species-enabled"];
        if (se.isScalar()) {
            m_do_species.assign(m_thermo->nSpecies(), se.asBool());
        } else {
            m_do_species = se.asVector<bool>(m_thermo->nSpecies());
        }
    }

    if (state.hasKey("radiation-enabled")) {
        m_do_radiation = state["radiation-enabled"].asBool();
        if (m_do_radiation) {
            m_epsilon_left = state["emissivity-left"].asDouble();
            m_epsilon_right = state["emissivity-right"].asDouble();
        }
    }

    if (state.hasKey("refine-criteria")) {
        const AnyMap& criteria = state["refine-criteria"].as<AnyMap>();
        double ratio = criteria.getDouble("ratio", m_refiner->maxRatio());
        double slope = criteria.getDouble("slope", m_refiner->maxDelta());
        double curve = criteria.getDouble("curve", m_refiner->maxSlope());
        double prune = criteria.getDouble("prune", m_refiner->prune());
        m_refiner->setCriteria(ratio, slope, curve, prune);

        if (criteria.hasKey("grid-min")) {
            m_refiner->setGridMin(criteria["grid-min"].asDouble());
        }
        if (criteria.hasKey("max-points")) {
            m_refiner->setMaxPoints(criteria["max-points"].asInt());
        }
    }

    if (state.hasKey("fixed-point")) {
        m_zfixed = state["fixed-point"]["location"].asDouble();
        m_tfixed = state["fixed-point"]["temperature"].asDouble();
    }
}

void StFlow::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    if (changed) {
        needJacUpdate();
    }
}

size_t StFlow::getSolvingStage() const
{
    throw NotImplementedError("StFlow::getSolvingStage",
        "Not used by '{}' objects.", type());
}

void StFlow::setSolvingStage(const size_t stage)
{
    throw NotImplementedError("StFlow::setSolvingStage",
        "Not used by '{}' objects.", type());
}

void StFlow::solveElectricField(size_t j)
{
    throw NotImplementedError("StFlow::solveElectricField",
        "Not used by '{}' objects.", type());
}

void StFlow::fixElectricField(size_t j)
{
    throw NotImplementedError("StFlow::fixElectricField",
        "Not used by '{}' objects.", type());
}

bool StFlow::doElectricField(size_t j) const
{
    throw NotImplementedError("StFlow::doElectricField",
        "Not used by '{}' objects.", type());
}

void StFlow::setBoundaryEmissivities(double e_left, double e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void StFlow::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    if (changed) {
        needJacUpdate();
    }
}

void StFlow::evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
{
    size_t j = m_points - 1;

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    rsd[index(c_offset_V,j)] = V(x,j);
    diag[index(c_offset_V,j)] = 0;
    double sum = 0.0;
    // set residual of poisson's equ to zero
    rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
    if (m_usesLambda) {
        rsd[index(c_offset_U, j)] = rho_u(x, j);
    } else {
        rsd[index(c_offset_U, j)] = rho_u(x, j) - rho_u(x, j-1);
    }

    rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j-1);
    diag[index(c_offset_L, j)] = 0;
    rsd[index(c_offset_T, j)] = T(x, j);
}

void StFlow::evalContinuity(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (m_usesLambda) {
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (m_isFree) {
        // terms involving V are zero as V=0 by definition
        if (grid(j) > m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1];
        } else if (grid(j) == m_zfixed) {
            if (m_do_energy[j]) {
                rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
            } else {
                rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                            - m_rho[0]*0.3); // why 0.3?
            }
        } else if (grid(j) < m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j];
        }
    } else {
        // unstrained with fixed mass flow rate
        rsd[index(c_offset_U, j)] = rho_u(x, j) - rho_u(x, j - 1);
    }
}

void StFlow::grad_hk(const double* x, size_t j)
{
    for(size_t k = 0; k < m_nsp; k++) {
        if (u(x, j) > 0.0) {
            m_dhk_dz(k,j) = (m_hk(k,j) - m_hk(k,j-1)) / m_dz[j - 1];
        }
        else {
            m_dhk_dz(k,j) = (m_hk(k,j+1) - m_hk(k,j)) / m_dz[j];
        }
    }
}

} // namespace
