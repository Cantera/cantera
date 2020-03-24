//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"

using namespace std;

namespace Cantera
{

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_epsilon_left(0.0),
    m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_kExcessLeft(0),
    m_kExcessRight(0),
    m_zfixed(Undef),
    m_tfixed(-1.)
{
    if (ph->type() == "IdealGas") {
        m_thermo = static_cast<IdealGasPhase*>(ph);
    } else {
        throw CanteraError("StFlow::StFlow",
                           "Unsupported phase type: need 'IdealGasPhase'");
    }
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
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
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(0, -1e20, 1e20); // no bounds on u
    setBounds(1, -1e20, 1e20); // V
    setBounds(2, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    setBounds(3, -1e20, 1e20); // lambda should be negative

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_L, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");
}

void StFlow::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}

void StFlow::setupGrid(size_t n, const doublereal* z)
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

void StFlow::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi");

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
    }
}

void StFlow::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void StFlow::_finalize(const doublereal* x)
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

    if (domainType() == cFreeFlow) {
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

void StFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
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

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const doublereal c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const doublereal c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
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
            if (doEnergy(0)) {
                rsd[index(c_offset_T,0)] = T(x,0);
            } else {
                rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
            }
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

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
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
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
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

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

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                double dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::updateTransport(doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
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

void StFlow::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

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

void StFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                doublereal sum = 0.0;
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
    case 0:
        return "velocity";
    case 1:
        return "spread_rate";
    case 2:
        return "T";
    case 3:
        return "lambda";
    case 4:
        return "eField";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const std::string& name) const
{
    if (name=="u") {
        warn_deprecated("StFlow::componentIndex",
                        "To be changed after Cantera 2.5. "
                        "Solution component 'u' renamed to 'velocity'");
        return 0;
    } else if (name=="velocity") {
        return 0;
    } else if (name=="V") {
        warn_deprecated("StFlow::componentIndex",
                        "To be changed after Cantera 2.5. "
                        "Solution component 'V' renamed to 'spread_rate'");
        return 1;
    } else if (name=="spread_rate") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else if (name == "eField") {
        return 4;
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

void StFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    vector<string> ignored;
    size_t nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);
    vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_fp x;
    size_t np = 0;
    bool readgrid = false, wrote_header = false;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            readgrid = true;
            setupGrid(np, x.data());
        }
    }
    if (!readgrid) {
        throw CanteraError("StFlow::restore",
                           "domain contains no grid points.");
    }

    debuglog("Importing datasets:\n", loglevel >= 2);
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_U,j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V") {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_V,j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_T,j)] = x[j];
            }

            // For fixed-temperature simulations, use the imported temperature
            // profile by default.  If this is not desired, call
            // setFixedTempProfile *after* restoring the solution.
            vector_fp zz(np);
            for (size_t jj = 0; jj < np; jj++) {
                zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
            }
            setFixedTempProfile(zz, x);
        } else if (nm == "L") {
            debuglog("lambda   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "lambda arary size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_L,j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (loglevel >=2 && !ignored.empty()) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        size_t nn = ignored.size();
        for (size_t n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    if (loglevel >= 1) {
        for (size_t ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }
    }

    if (dom.hasChild("energy_enabled")) {
        getFloatArray(dom, x, false, "", "energy_enabled");
        if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_energy[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            throw CanteraError("StFlow::restore", "energy_enabled is length {}"
                               "but should be length {}", x.size(), nPoints());
        }
    }

    if (dom.hasChild("species_enabled")) {
        getFloatArray(dom, x, false, "", "species_enabled");
        if (x.size() == m_nsp) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_species[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            // This may occur when restoring from a mechanism with a different
            // number of species.
            if (loglevel > 0) {
                warn_user("StFlow::restore", "species_enabled is "
                    "length {} but should be length {}. Enabling all species "
                    "equations by default.", x.size(), m_nsp);
            }
            m_do_species.assign(m_nsp, true);
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }

    if (domainType() == cFreeFlow) {
        getOptionalFloat(dom, "t_fixed", m_tfixed);
        getOptionalFloat(dom, "z_fixed", m_zfixed);
    }
}

XML_Node& StFlow::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    flow.addAttribute("type",flowType());

    if (m_desc != "") {
        addString(flow,"description",m_desc);
    }
    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");

    addFloatArray(gv,"z",m_z.size(), m_z.data(),
                  "m","length");
    vector_fp x(soln.nColumns());

    soln.getRow(c_offset_U, x.data());
    addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

    soln.getRow(c_offset_V, x.data());
    addFloatArray(gv,"V",x.size(),x.data(),"1/s","rate");

    soln.getRow(c_offset_T, x.data());
    addFloatArray(gv,"T",x.size(),x.data(),"K","temperature");

    soln.getRow(c_offset_L, x.data());
    addFloatArray(gv,"L",x.size(),x.data(),"N/m^4");

    for (size_t k = 0; k < m_nsp; k++) {
        soln.getRow(c_offset_Y+k, x.data());
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),x.data(),"","massFraction");
    }
    if (m_do_radiation) {
        addFloatArray(gv, "radiative_heat_loss", m_z.size(),
            m_qdotRadiation.data(), "W/m^3", "specificPower");
    }
    vector_fp values(nPoints());
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = m_do_energy[i];
    }
    addNamedFloatArray(flow, "energy_enabled", nPoints(), &values[0]);

    values.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        values[i] = m_do_species[i];
    }
    addNamedFloatArray(flow, "species_enabled", m_nsp, &values[0]);

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    if (m_zfixed != Undef) {
        addFloat(flow, "z_fixed", m_zfixed, "m");
        addFloat(flow, "t_fixed", m_tfixed, "K");
    }
    return flow;
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

void StFlow::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
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
    doublereal sum = 0.0;
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
    if (domainType() == cAxisymmetricStagnationFlow) {
        rsd[index(c_offset_U,j)] = rho_u(x,j);
        if (m_do_energy[j]) {
            rsd[index(c_offset_T,j)] = T(x,j);
        } else {
            rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
        }
    } else if (domainType() == cFreeFlow) {
        rsd[index(c_offset_U,j)] = rho_u(x,j) - rho_u(x,j-1);
        rsd[index(c_offset_T,j)] = T(x,j) - T(x,j-1);
    }
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
    if (domainType() == cAxisymmetricStagnationFlow) {
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (domainType() == cFreeFlow) {
        if (grid(j) > m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1]
                - (density(j-1)*V(x,j-1) + density(j)*V(x,j));
        } else if (grid(j) == m_zfixed) {
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
    }
}

} // namespace
