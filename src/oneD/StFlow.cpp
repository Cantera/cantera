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
#include "cantera/base/Parser.h"
#include <iostream>

using namespace std;

namespace Cantera
{

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t nsoot, size_t neq, size_t points) :
    Domain1D(nsp, points, 0.0, nsoot, neq),
    m_nsp(nsp),
    m_neq(neq),
    m_nsoot(nsoot)
{
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }
    m_thermo = ph;

    size_t nsp2 = m_thermo->nSpecies();
    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.


    m_nv = m_neq + m_nsp + m_nsoot;
    if (nsp2 + m_nsoot != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nv, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // set pressure based on associated thermo object
    setPressure(m_thermo->pressure());

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    // m_nv = c_offset_Y + m_nsp;

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
    if (m_neq == 1) {
        setBounds(c_offset_Tflamelet, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    } else {
        setBounds(c_offset_U, -1e20, 1e20); // no bounds on u
        setBounds(c_offset_V, -1e20, 1e20); // V
        setBounds(c_offset_T, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
        setBounds(c_offset_L, -1e20, 1e20); // lambda should be negative
        setBounds(c_offset_E, -1e20, 1e20); // no bounds for inactive component
    }
    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(m_neq+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    if (m_neq == 1) {
        m_refiner->setActive(c_offset_Tflamelet, false);
    } else {
        m_refiner->setActive(c_offset_U, false);
        m_refiner->setActive(c_offset_V, false);
        m_refiner->setActive(c_offset_T, false);
        m_refiner->setActive(c_offset_L, false);
    }

    vector<double> gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");

    // P. Wolf
    //---------------------AVBP PEA initialization-------------------
    AVBPReadInputPea();

    // B. Franzelli
    //---------------------AVBP thickening initialisation------------
    AVBPReadInputChem();

    // E. Lameloise
    //---------------------Soot computation initialisation-----------
    if (m_nsoot > 0){
        initSoot();
        }
}

StFlow::StFlow(shared_ptr<ThermoPhase> th, size_t nsp, size_t nsoot, size_t neq, size_t points)
    : StFlow(th.get(), nsp, nsoot, neq, points)
{
    m_solution = Solution::create();
    m_solution->setThermo(th);
}

StFlow::StFlow(shared_ptr<Solution> sol, const string& id, size_t nsoot, size_t neq, size_t points)
    : StFlow(sol->thermo().get(), sol->thermo()->nSpecies(), nsoot, neq, points)
{
    m_solution = sol;
    m_id = id;
    m_kin = m_solution->kinetics().get();
    m_trans = m_solution->transport().get();
    if (m_trans->transportModel() == "none") {
        // @deprecated
        warn_deprecated("StFlow",
            "An appropriate transport model\nshould be set when instantiating the "
            "Solution ('gas') object.\nImplicit setting of the transport model "
            "is deprecated and\nwill be removed after Cantera 3.0.");
        setTransportModel("mixture-averaged");
    }
    m_solution->registerChangedCallback(this, [this]() {
        setKinetics(m_solution->kinetics());
        setTransport(m_solution->transport());
    });
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
    m_h.resize(m_points, 0.0);
    m_hr.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);
    avbp_thick.resize(m_points,0.0);

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
    
    // P. Wolf
    m_zmixfrac.resize(m_points);
    // Soot
    if (m_nsoot > 0){
        sootConsumption.resize(m_nsp, 0.0);
        m_soot_diff.resize(m_nsoot, m_points, 0.0);
        m_soot_soret.resize(m_nsoot, m_points, 0.0);
        m_qdotNucleation.resize(m_points, 0.0);
        m_qdotCondensation.resize(m_nsoot, m_points, 0.0);
        m_qdotCoagulation.resize(m_nsoot, m_points, 0.0);
        m_qdotSg.resize(m_nsoot, m_points, 0.0);
        m_qdotOxidation.resize(m_nsoot, m_points, 0.0);
    }
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

void StFlow::setupGrid(size_t n)
{
    resize(m_nv, n);

    m_z[0] = 0.0;
    for (size_t j = 1; j < m_points; j++) {
	m_z[j] = (double) j / ( (double) (n-1) );
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
    // P. Wolf
    if ((avbp_ipea==1)||(avbp_ipea==2)||(avbp_ipea==3)){
        calcMixFrac(x,j);
    }
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
    // P. Wolf
    if ((avbp_ipea==1)||(avbp_ipea==2)||(avbp_ipea==3)){
        calcMixFrac(x,j);
    }
}

    /// Calculates and updates the mixture fraction at point j
    /// WARNING: Only correct if no Carbon atoms in the oxidizer line
    /// and pure fuel in the fuel line
    /// Added by P. Wolf - March 2010
    /// Thanks to Alireza Najafiyazdi

void StFlow::calcMixFrac(const double* x, size_t j) {
    	size_t id_C = this->m_thermo->elementIndex("C"); // Index of C atom
    	size_t id_fuel = this->m_thermo->speciesIndex(avbp_fuel);  // Index of fuel
    	double Z_C = 0.0; // passive scalar for C atom
    	double Z_C_F;     // value of the passive scalar in fuel line
    	double W_C = this->m_thermo->atomicWeight(id_C);  // Molecular weight of C atom
    	double W_fuel = this->m_thermo->molecularWeight(id_fuel);
    	double n_C;     // Number of C atoms in each species

    	for (size_t k=0; k<this->m_nsp; k++) {
    		n_C = this->m_thermo->nAtoms(k,id_C);
    		Z_C += n_C*x[index(c_offset_Y + k, j)]/this->m_thermo->molecularWeight(k);
    	}
    	Z_C = Z_C*W_C;
    	Z_C_F = W_C*this->m_thermo->nAtoms(id_fuel,id_C)/W_fuel;
		m_zmixfrac[j] = Z_C/Z_C_F;    // ie Z_C=0 in oxidizer line
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
    if (m_nsoot > 0){
        updateSootDiffFluxes(x, j0, j1);
    }
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
            // absorption coefficient for soot
            if (m_nsoot != 0 && m_do_soot_radiation) {
                double k_P_soot = 0.0;
                k_P_soot = 3.83 * CRad0 * fv(x,j) * T(x,j) / CRad2;
                k_P += k_P_soot;
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
        if (m_nsoot > 0){
            getDistributionOrdre0(x,j);
        }
        //avbp_thick.resize(m_points,avbp_fthick); ///////////////////////////// NEW to handle pytest
        AVBPcompute_local_thick(x,j);
        double thick_j = avbp_thick[j];

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
// no thick         -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
                    -(m_flux(k,0) * thick_j + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;
            if (m_nsoot > 0){
                for (size_t k = 0; k < m_nsoot; k++){
                    rsd[index(c_offset_S + k, 0)] = - (m_soot_diff(k,0)
                                                    + rho_u(x,0)* Ys(x,k,0));
                }
            }

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
            AVBPcompute_local_thick(x,j);
            double thick_j = avbp_thick[j];
            double inv_thick = 1.0 / thick_j;
            double thick_prev = avbp_thick[j-1];

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
            //    Soot equation
            //
            //   \rho u dYs/dz + d(-0.554 \mu 1/T dT/dx Ys)/dz
            //   = d(\rho Ds dYs/dz)/dz + \rhos \rho qdots
            //
            //   convection - Thermophoresis = Diffusion + Source
            //-------------------------------------------------
            if (m_nsoot > 0) {
                sootSource(x,j);
                for (size_t k = 0; k < m_nsoot; k++){
                    // Convection
                    double soot_convec = rho_u(x,j)*dYsdz(x,k,j);
                    // Diffusion
                    double soot_diffus = 2.0 * (m_soot_diff(k,j) - m_soot_diff(k,j-1))
                                         / (z(j+1) - z(j-1));
                    // Thermophoresis
                    double soot_soret = 2.0 * (m_soot_soret(k,j) - m_soot_soret(k,j-1))
                                        / (z(j+1) - z(j-1));
                    // Source terms [m3/kg/s]
                    double soot_source = 0;
                    soot_source += (k == 0 ? m_qdotNucleation[j] : 0.0);
                    soot_source += (m_do_condensation ? m_qdotCondensation(k,j) : 0.0);
                    soot_source += (m_do_coagulation ? m_qdotCoagulation(k,j) : 0.0);
                    soot_source += (m_do_sg ? m_qdotSg(k,j) : 0.0);
                    soot_source += (m_do_oxidation ? m_qdotOxidation(k,j) : 0.0);
                    // kg/s/m^3
                    soot_source *= (rho_soot * m_rho[j]);
                    if(k == m_nsoot - 1 && m_trash_section){soot_source = 0.0;}
                    // Residual
                    rsd[index(c_offset_S + k, j)] =
                      (soot_source - soot_convec + soot_diffus + soot_soret) / m_rho[j]
                      - rdt * (Ys(x,k,j) - Ys_prev(k,j));
                    diag[index(c_offset_S + k, j)] = 1;
                }
            }
            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            if (m_nsoot > 0 && m_do_retroaction){
                for (size_t k=0; k < m_nsp; k++){
                    // (kmol.m^(-3).s^(-1))
                    m_wdot(k,j) -= sootConsumption[k];
                }
            }
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
// no thick                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                double diffus = 2.0*(m_flux(k,j)*thick_j - m_flux(k,j-1)*thick_prev)
                                / (z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j)*inv_thick)
// no thick     = (m_wt[k]*(wdot(k,j))
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
                double dtdzj = dTdz(x,j);
                double sum = 0.0;
                double sum2 = 0.0;

                grad_hk(x, j);
                for (size_t k = 0; k < m_nsp; k++) {
// no thick         double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    double flxk = 0.5*(m_flux(k,j-1)*thick_prev + m_flux(k,j)*thick_j);
                    sum += wdot(k,j)*m_hk(k,j);
                    sum2 += flxk * m_dhk_dz(k,j) / m_wt[k];
                }

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
// no thick                               //- divHeatFlux(x,j) - sum;
                                            - AVBPdivHeatFlux(x,j) - sum*inv_thick - sum2;

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
        } else if (n >= c_offset_S && n < (c_offset_S + m_nsoot)) {
            return sectionName(n-c_offset_S);
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
    } else if (find(m_section_name.begin(), m_section_name.end(), name) != m_section_name.end()) {
        return stoi(name.substr(2))+c_offset_S; // offset not required in previous version, why?
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
    vector<bool> did_section(m_nsoot,false);
    bool wrote_header = false;

    const auto grid = arr.getComponent("grid").as<vector<double>>();
    setupGrid(nPoints(), &grid[0]);

    for (size_t i = 0; i < nComponents(); i++) {
        if (!componentActive(i)) {
            continue;
        }
        string name = componentName(i);
        if (arr.hasComponent(name)) {
            const vector<double> data = arr.getComponent(name).as<vector<double>>();
            if (find(m_section_name.begin(), m_section_name.end(), name) != m_section_name.end()){
                    size_t n = stoi(name.substr(2));
                    did_section[n]=true;
                    writelog("Adding soot section to soln "+name+"\n");
                }
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(i,j)] = data[j];
            }
        } else if (find(m_section_name.begin(), m_section_name.end(), name) != m_section_name.end()){
            size_t n = stoi(name.substr(2));    
            did_section[n]=true;
            writelog("Creating soot section in soln "+name+"\n");
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(n+c_offset_S,j)] = 0.0;
            }
        } else {
            warn_user("StFlow::fromArray", "Saved state does not contain values for "
                "component '{}' in domain '{}'.", name, id());
        }
    }

    if (soot_loglevel >= 1 && m_nsoot != 0) {
        wrote_header = false;
        for (size_t ks = 0; ks < m_nsoot; ks++) {
            if (did_section[ks] == false) {
                if (!wrote_header) {
                    writelog("\n\n");
                    writelog("Missing data for soot sections:\n");
                    wrote_header = true;
                }
                writelog(m_section_name[ks]+" ");
            }
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
    } else if (m_type == cFlameletFlow){
        return "Flamelet";
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
/* ------------------------------------------------------------------------
	Reads the mixture_database.dat file and sets accordingly the PEA function(s)
    P. Wolf modif by J. Wirtz
--------------------------------------------------------------------------*/
void StFlow::AVBPReadInputPea() {
  Parser parser;
  Param* param;

  ifstream avbp_mixturedb("mixture_database.dat");
  if (avbp_mixturedb.good()){
    parser.parseFile("mixture_database.dat");
    // In case of pea, the name of the gas should be the same than the one in the mixture_database.
    // Get number of mixtures in database
    size_t n_mixtures = parser.nbParamOccurence("mixture_name");
    size_t idx_beg = std::string::npos;
    size_t idx_end = std::string::npos;
    // Loop through all occurences of keyword "mixture name"
    for (size_t i=0; i< n_mixtures; ++i) {
      param = parser.getParam("mixture_name",i);
      std::string current_str = param->readString(0);
      // Check if mixture_name is the one requested in the cti file
      if ( current_str == m_thermo->name() ) {
        // Store bounding idx to get important info
        idx_beg = parser.getParamNumber("mixture_name",i);
	      idx_end = parser.getParamNumber("mixture_name",i+1);

        // Get the fuel of the 2-step chemistry
        param = parser.getParam("species_name", idx_beg, idx_end);
        for (size_t i=0; i<m_nsp; ++i) {
          std::string avbp_species = param->readString(i);
          if (avbp_species.find('C')!=string::npos && avbp_species.find('H')!=string::npos){
            avbp_fuel = avbp_species;
          }
          if (avbp_species == "H2"){
            avbp_fuel = avbp_species;
          }

        }

        // Check if PEA method and search for the coefficient
        if (parser.checkParam("pea_method", idx_beg, idx_end) == true){
          param = parser.getParam("pea_method", idx_beg, idx_end);
          std::string avbp_pea = param->readString(0);
          if (avbp_pea == "pea1"){
            avbp_ipea = 1;
            cout<<"ipea set to 1, PEA_1 will be applied"<<endl;
				    avbp_pea_coeffs.resize(7,0.0);
            param = parser.getParam("pea_coeff", idx_beg, idx_end);
            for (size_t i=0; i<7; i++) {
				      avbp_pea_coeffs[i] = param->readDouble(i);
			      }
          }
          else if (avbp_pea == "pea2") {
            avbp_ipea = 2;
            cout<<"ipea set to 2, PEA_2 will be applied"<<endl;
				    avbp_pea_coeffs.resize(18,0.0);
            param = parser.getParam("pea_coeff", idx_beg, idx_end, 0);
            for (size_t i=0; i<8; i++) {
				      avbp_pea_coeffs[i] = param->readDouble(i);
			      }
            param = parser.getParam("pea_coeff", idx_beg, idx_end, 1);
            for (size_t i=0; i<10; i++) {
				      avbp_pea_coeffs[i+8] = param->readDouble(i);
			      }

          }
          else{
					  throw CanteraError("readAvbpInputPEA","ipea is set to a weird value");
          }
        }
        else {
          avbp_ipea = 0;
        }
      }
      else {
        // pea not considered as mixture_name (mixture_database.dat) does not correspond to gas name (cti)
        avbp_ipea = 0;
      }

      // calculate the global equivalence ratio phi_cst
      if ((avbp_ipea==1)||(avbp_ipea==2)){
        size_t id_C = this->m_thermo->elementIndex("C"); // Index of C atom
        size_t id_H = this->m_thermo->elementIndex("H"); // Index of H atom
        size_t id_fuel = this->m_thermo->speciesIndex(avbp_fuel);  // Index of fuel
        size_t id_ox = this->m_thermo->speciesIndex("O2");    // Index of O2
        size_t n_C = this->m_thermo->nAtoms(id_fuel,id_C);
        size_t n_H = this->m_thermo->nAtoms(id_fuel,id_H);
        //Calculate in mass units phi = s_m*Yf/Yo
        phi_cst = (n_C + n_H/4)*this->m_thermo->molecularWeight(id_ox)/this->m_thermo->molecularWeight(id_fuel);
        phi_cst = phi_cst*1./0.2333;			// WARNING: valid only for Fuel/AIR mixture
      }

      if (avbp_ipea!=0){
        break;
      }

    }
  }
  else {
    avbp_ipea = 0;
  }
}
/**********
    Reading the AVBP input_chem.dat
  B. Franzelli
*************/
void StFlow::AVBPReadInputChem() {
  avbp_fthick = m_thick;
  if (avbp_fthick > 1.0) {
    cout << "INFO: Thickening applied with constant value." << endl;
  }
  else {
    avbp_fthick = 1.0;
  }
}
void StFlow::AVBPcompute_local_thick(double* x,size_t j){
  avbp_thick[j] = avbp_fthick;
}

//----------------------//
// SOOT RELATED MTEHODS //
//----------------------//
void StFlow::initSoot(){
    //----------------------------------------------------
    // Initializes soot computation
    //----------------------------------------------------

    // Sets offset of soot sections in solution array
    // Resizes necessay arrays
    // Sets sections names "YsXX", disables refinement on soot and sets clipping

    // Offset definition
    c_offset_S = c_offset_Y + m_nsp;
    // Array allocation
    m_section_name.resize(m_nsoot);
    sootConsumption.resize(m_nsp, 0.0);
    q.resize(m_nsoot,0.0);
    vSectMin.resize(m_nsoot, 0.0);
    vSectMax.resize(m_nsoot, 0.0);
    vSectMean.resize(m_nsoot, 0.0);
    sSectMean.resize(m_nsoot, 0.0);
    dSectMean.resize(m_nsoot, 0.0);
    dSectCol.resize(m_nsoot, 0.0);
    kfractal.resize(m_nsoot, 1.0);
    Dfractal.resize(m_nsoot, 3);
    nNucMean.resize(m_nsoot, 0.0);
    rNucMean.resize(m_nsoot, 0.0);
    rSectCol.resize(m_nsoot, 0.0);
    aSectCol.resize(m_nsoot, 0.0);
    theta_surf.resize(m_nsoot,2.0);
    m_soot_diff.resize(m_nsoot, m_points, 0.0);
    m_soot_soret.resize(m_nsoot, m_points, 0.0);
    m_qdotNucleation.resize(m_points, 0.0);
    m_qdotCondensation.resize(m_nsoot, m_points, 0.0);
    m_qdotCoagulation.resize(m_nsoot, m_points, 0.0);
    collision_mat.resize(m_nsoot, m_nsoot, 0.0);
    m_qdotSg.resize(m_nsoot, m_points, 0.0);
    m_qdotOxidation.resize(m_nsoot, m_points, 0.0);
    // Sections setup
    for (size_t k = 0; k < m_nsoot; k++) {
        std::string string_k = std::to_string(k);
        m_section_name[k] = "Ys" + string_k;
        m_refiner->setActive(c_offset_S+k, false);
        setBounds(c_offset_S + k, -1.0e-5, 1e5);
    }
    
}

void StFlow::setPrecursors(vector<size_t> id_precursors){
    //----------------------------------------------------
    // Recovers precursor data for further computation
    //----------------------------------------------------

    if (m_nsoot <= 0){
        writelog("\n--> CANTERA SOOT : no sections :(\n");
        return;
    }
    size_t id_C = m_thermo->elementIndex("C");
    W_C = m_thermo->atomicWeight(id_C);
    n_PAH = id_precursors.size();
    m_precursors.resize(n_PAH);
    W_PAH.resize(n_PAH,0.0);
    n_C.resize(n_PAH,0.0);
    for (size_t k = 0; k != n_PAH; k++){
        m_precursors[k] = id_precursors[k] - c_offset_Y;
        W_PAH[k] = m_thermo->molecularWeight(m_precursors[k]); 
        n_C[k] = m_thermo->nAtoms(m_precursors[k],id_C);
    }
}

void StFlow::finalizeSoot(){
    //----------------------------------------------------
    // Soot setup finalization
    //----------------------------------------------------
    // Creates soot volumes
    sootCreationVol();
    // Computes collision matrix
    sootCollisionInverse();
    // Load surface reactions
    if (m_do_sg || m_do_oxidation){
        loadHaca(m_haca_model);
        loadSurface();}

    if (soot_loglevel >= 1){showSootInfo();}
    if (soot_loglevel >= 2){showSootSections();}
}

void StFlow::sootCreationVol(){
    // Sections creations
    // PAH list boundaries
    size_t n_C_min = *min_element(n_C.begin(),n_C.end());
    size_t n_C_max = *max_element(n_C.begin(),n_C.end());
    // First section lower boundary : smallest possible dimer
    vSectMin[0] = 2.0 * (double)n_C_min * V_C2;
    // First section upper boundary : biggest possible dimer 
    vSectMax[0] = V_C2 * max(2.0 * n_C_min + n_C_max, 2.0 * n_C_max);
    // Following sections : geometrical progression
    vSectMin[1] = vSectMax[0];
    if (m_trash_section){
        double vTrashLowLim = Pi * pow(dTrashLowLim,3.0) / 6;
        for (size_t k = 1; k < m_nsoot - 2; k++){
            double exponent = double(k) / double(m_nsoot-2);
            vSectMax[k] = vSectMin[1] * pow(vTrashLowLim / vSectMin[1],exponent);
            vSectMin[k+1] = vSectMax[k];
        }
        vSectMax[m_nsoot - 2] = vTrashLowLim;
        vSectMin[m_nsoot - 1] = vTrashLowLim;
    }else{
        for (size_t k = 1; k < m_nsoot - 1; k++){
            double exponent = double(k) / double(m_nsoot-1);
            vSectMax[k] = vSectMin[1] * pow(bigSoot / vSectMin[1],exponent);
            vSectMin[k+1] = vSectMax[k];
        }
    }

    vSectMax[m_nsoot-1] = bigSoot;
    // Building mean values
    for (size_t k = 0; k < m_nsoot; k++){
        vSectMean[k] =  (vSectMax[k] + vSectMin[k]) / 2.0 ; //(vSectMax[k] - vSectMin[k]) / log (vSectMax[k] / vSectMin[k]);
        dSectMean[k] = pow(6 / Pi * vSectMean[k], 1.0 / 3.0);
    }

    // Sections morphology
    for (size_t k = 0; k < m_nsoot; k++){
        if (vSectMean[k] >= vlim1){
            // Fractal dimension
            Dfractal[k] = Df;
            // Fractal prefactor
            kfractal[k] = 4.46 * pow(Dfractal[k], -2.08);
        }
    }

    if (m_soot_morphology == "sphere" || m_soot_morphology == "rodrigues"){
        // Soot morphological properties as defined by P. Rodrigues (PhD Thesis, 2018)
        // For spheres, theta = 2.0
        for (size_t k = 0; k < m_nsoot; k++){
            if (m_soot_morphology == "rodrigues" && vSectMin[k] >= vlim1){
                theta_surf[k] = (3.0 * (log(vSectMean[k] / vlim1))
                                + 2.0 * (log(vlim1 / V_C2)))
                                / (log(vSectMean[k] / V_C2));
                if (k == m_nsoot-1 && m_trash_section){theta_surf[k] = theta_surf[k-1];} 
            }
            // Collision diameter
            dSectCol[k] = 6 / pow(36 * Pi, 1.0 /Dfractal[k])
                        * pow(vSectMean[k], (Dfractal[k]-2.0) / Dfractal[k])
                        * pow(S_C2 *
                        pow(vSectMean[k]/V_C2,theta_surf[k] / 3.0),(3.0-Dfractal[k]) / Dfractal[k]);
            // Collision radius
            rSectCol[k] = dSectCol[k] / 2.0;
            // Collision cross-section
            aSectCol[k] = Pi * pow(rSectCol[k], 2.0);
            // Surface area
            sSectMean[k] = S_C2 * pow(vSectMean[k] / V_C2, theta_surf[k] / 3.0);
            // Nuclei number
            nNucMean[k] = pow(sSectMean[k], 2) / (36. * Pi * pow(vSectMean[k], 2.));
            // Nuclei radius[m]
            rNucMean[k] = 3. * vSectMean[k] / sSectMean[k];
        }
        // // From P. Rodrigues, trick to help numerical stability (not tested its impact)
        // dSectCol[m_nsoot-1] = dSectMean[m_nsoot-1];
        // rSectCol[m_nsoot-1] = dSectCol[m_nsoot-1] / 2.0;
        // aSectCol[m_nsoot-1] = Pi * pow(rSectCol[m_nsoot-1], 2.0);
    }else if (m_soot_morphology == "thajudeen"){
        // Soot morphological properties as defined by Thajudeen (2012)
        for (size_t k = 0; k < m_nsoot; k++){
            // Nuclei number
            nNucMean[k] = max(1., vSectMean[k] / vlim1);
            // Nuclei radius
            rNucMean[k] = min(pow(3.0 * vSectMean[k] / (4.0 * Pi), 1./3.), pow(3.0 * vlim1 / (4.0 * Pi), 1.0/3.0)); 
            if (vSectMean[k] >= vlim1){
                // Collision radius (Smoluchowski radius)
                double alpha1 = 0.253 * pow(Dfractal[k], 2.0) - 1.209 * Dfractal[k] + 1.433;
                double alpha2 = -0.218 * pow(Dfractal[k], 2.0) + 0.964 * Dfractal[k] - 0.180; 
                double phi_r = 1 / (alpha1 * log(nNucMean[k]) + alpha2);
                rSectCol[k] = rNucMean[k] * phi_r * pow(nNucMean[k] / kfractal[k], 1/Dfractal[k]);
                // Collision cross-section (Projected Area)
                double alpha3 = 0.439 * pow(Dfractal[k], 2.0) - 2.221 * Dfractal[k] + 2.787;
                double alpha4 = -0.232 * pow(Dfractal[k], 3.0) + 1.273 * pow(Dfractal[k], 2.0) - 2.183 * Dfractal[k] + 1.906;
                double phi_p = pow(nNucMean[k], -alpha3) / alpha4;
                aSectCol[k] = Pi * pow(rNucMean[k], 2.) * phi_p * pow(nNucMean[k] / kfractal[k], 2.0/Dfractal[k]);
            }else{
                rSectCol[k] = dSectMean[k] / 2.;
                aSectCol[k] = Pi * pow(rSectCol[k], 2.);
            }
            //Collision diameter ("Smoluchowski diameter")
            dSectCol[k] = 2. * rSectCol[k];
            // Surface area (computed from smoluchowski radius, close to Rodrigues definition)
            sSectMean[k] = Pi * pow(dSectCol[k], 2.);
            // theta parameter (needed for surface reactions)
            theta_surf[k] = 3. * log(sSectMean[k]/S_C2) / log(vSectMean[k]/V_C2);
        }
    }
}   

void StFlow::sootCollisionInverse(){
    //----------------------------------------------------
    // Computes colision matrix (coagulation)
    //----------------------------------------------------

    for (size_t k = 0; k < m_nsoot; k++){
        for (size_t l = k; l < m_nsoot; l++){
            for (size_t m = l; m < m_nsoot; m++){
                double vol_min =  vSectMin[l] + vSectMin[k];
                if(vol_min <= vSectMax[m] && vol_min >= vSectMin[m]){
                    collision_mat(k,l) =  m;
                }
            }
        }
    }
}

void StFlow::loadSurface(){
    //----------------------------------------------------
    // Resizes and compute surface reactions arrays
    //----------------------------------------------------

    kpower.resize(m_nsoot,0.0);
    mpower.resize(m_nsoot,0.0);
    vc2power.resize(m_nsoot,0.0);
    vc2powerm.resize(m_nsoot,0.0);
    vmax_kpower.resize(m_nsoot,0.0);
    vmin_kpower.resize(m_nsoot,0.0);
    vmax_mpower.resize(m_nsoot,0.0);
    vmin_mpower.resize(m_nsoot,0.0);
    vmaxmc2_kpower.resize(m_nsoot,0.0);
    vminpc2_kpower.resize(m_nsoot,0.0);
    vmaxmc2_mpower.resize(m_nsoot,0.0);
    vminpc2_mpower.resize(m_nsoot,0.0);
    vc2powervect.resize(m_nsoot,0.0);
    vc2powermvect.resize(m_nsoot,0.0);
    for (size_t k = 0; k < m_nsoot; k++){
        kpower[k] = (theta_surf[k] + 3.0)/3.0;
        mpower[k] = theta_surf[k]/3.0;
        vc2power[k] = (3.0 - theta_surf[k] )/3.0;
        vc2powerm[k] = - mpower[k];
        vmax_kpower[k] = pow(vSectMax[k],kpower[k]);
        vmin_kpower[k] = pow(vSectMin[k],kpower[k]);
        vmax_mpower[k] = pow(vSectMax[k],mpower[k]);
        vmin_mpower[k] = pow(vSectMin[k],mpower[k]);
        vmaxmc2_kpower[k] = pow(vSectMax[k] - V_C2,kpower[k]);
        vminpc2_kpower[k] = pow(vSectMin[k] + V_C2,kpower[k]);
        vmaxmc2_mpower[k] = pow(vSectMax[k] - V_C2,mpower[k]);
        vminpc2_mpower[k] = pow(vSectMin[k] + V_C2,mpower[k]);
        vc2powervect[k] = pow(V_C2,vc2power[k]);
        vc2powermvect[k] = pow(V_C2,vc2powerm[k]);
    }
}

void StFlow::showSootInfo(){
    writelog("\n{:=^47}{:9}{:=^47}", "", "SOOT INFO", "");
    writelog("\n{:21}{:2}{}", "NUMBER OF SECTIONS", ":", m_nsoot);
    writelog("\n{:21}{:1}","PRECURSORS", ":");
    for (size_t k = 0; k != n_PAH; k++){
        writelog(" {}", m_thermo->speciesName(m_precursors[k]));
    }
    if (dTrashLowLim > 0){writelog("\n{:21}{:2}{} m", "MAX DIAMETER", ":", dTrashLowLim);}
    writelog("\n{:21}{:2}{}", "AGGREGATES MORPHOLOGY",  ":", m_soot_morphology);
    writelog("\n{:21}{:2}{}", "RETROACTION", ":", m_do_retroaction);
    writelog("\n{:21}{:2}{}", "CONDENSATION", ":", m_do_condensation);
    writelog("\n{:21}{:2}{}", "COAGULATION", ":", m_do_coagulation);
    writelog("\n{:21}{:2}{}", "SURFACE GROWTH", ":", m_do_sg);
    writelog("\n{:21}{:2}{}", "OXIDATION", ":", m_do_oxidation);
    if (m_do_sg || m_do_oxidation){
        writelog("\n{:21}{:2}", "HACA MODEL",  ":");
        if (m_haca_model == 1){
            writelog("Mauss");
        }else if (m_haca_model == 2){
            writelog("Blanquart");
        }else if (m_haca_model == 3){
            writelog("Kazakov"); 
            writelog("\n{:21}{:2}{}{:2}", "-> Flame temperature", ":", kazakovTad, " K");
        }
    }
    writelog("\n{:21}{:2}{}", "RADIATIVE HEAT LOSSES", ":", m_do_soot_radiation);
    writelog("\n{:21}{:2}{}", "SORET EFFECT", ":", m_do_soot_soret);
    writeline('=', 103, false, true);
    writelog("\n");
}
void StFlow::showSootSections(){
  //----------------------------------------------------
  // Displays section informations
  //----------------------------------------------------

  if (m_nsoot <= 0){
    writelog ("\nSOOT INFO : no sections to show :(");
    return;
  }
  writelog("\n{:=^36}{:30}{:=^37}", "", "SOOT SECTIONS INFO (SI UNITS)", "");
  writelog("\n Section            Vol. min.         Vol. mean         Vol. max.         Diam. mean         Diam. coll");
  writeline('-', 103, false, true);
  for (size_t k = 0; k < m_nsoot; k++) {
      writelog("\n {:10d}        {:10.4g}        {:10.4g}        {:10.4g}        {:10.4g}        {:10.4g}",
               k,vSectMin[k],vSectMean[k],vSectMax[k],dSectMean[k],dSectCol[k]);
  }
  writeline('=', 103, false, true);
  writelog("\n");
}

void StFlow::sootSource(const double* x, size_t j) {
    //----------------------------------------------------
    // Computes soot source terms
    //----------------------------------------------------
    // Dimerization
    sootDimerization(x,j);
    // Nucleation
    sootNucleation(x,j);
    // Condensation
    if (m_do_condensation){sootCondensation(x,j);}
    // Coagulation
    if (m_do_coagulation){sootCoagulation(x,j);}
    // Surface chemistry
    if (m_do_sg || m_do_oxidation){sootSurface(x,j);}
}

void StFlow::sootDimerization(const double* x, size_t j){
  //----------------------------------------------------
  // Computes dimerization source term
  //----------------------------------------------------
    double omega = 0, c_dimer = 0;
    double V_k, d2_k, omega_k, C_k, gamma_k, k_k;
    for (size_t k = 0; k < n_PAH; k++){
        // PAH volume [m3]
        V_k = n_C[k] * V_C2 / 2.0;
        // Squared PAH diameter [m2]
        d2_k = pow(6.0 * V_k / Pi, 2.0/3.0);
        // PAH concentraiton [kmol/m3]
        C_k = m_rho[j] * Y(x,m_precursors[k],j) / W_PAH[k];
        // PAH sticking coefficient [-]
        gamma_k = Cn_cst * pow(W_PAH[k], 4.0);
        // Reaction rate constant [m3/kmol/s]
        k_k = sqrt(4.0 * Pi * Boltzmann * T(x,j) / (W_C * n_C[k] / Avogadro)) * d2_k * Avogadro;
        // PAH dimerization reaction rate [kmol/m3/s]
        omega_k = gamma_k * k_k * pow(C_k, 2.0);
        // PAH consumption [kmol/m3/s] (x2 because 2 mol PAH -> 1 mol dimer)
        sootConsumption[m_precursors[k]] = 2 * omega_k;
        // Adds PAH contribution to overall local dimerization rate 
        omega += omega_k;
        // Adds PAH contribution to overall dimer carbon atoms
        c_dimer += omega_k * n_C[k];
        
    }
    // Overall dimerization rate [1/kg/s]
    r_dimer = omega * Avogadro / m_rho[j];
    // Overall dimer volume [m3]
    V_dimer = c_dimer * V_C2 / omega;

    // double K00 = 0, cdimer = 0;
    // double V_k, d_k, kdimfwd_k, K00_k, C_k;
    // for (size_t k = 0; k < n_PAH; k++){
    //     // PAH volume [m3]
    //     V_k = n_C[k] * V_C2 / 2.0;
    //     // PAH diameter [m]
    //     d_k = pow(6.0 * V_k / Pi, 1.0/3.0);
    //     // PAH concentraiton [mol/m3]
    //     C_k = m_rho[j] * Y(x,m_precursors[k],j) / W_PAH[k];
    //     // m6/s/kg/mol2
    //     kdimfwd_k = Cn_cst * pow(W_PAH[k], 4.0) / m_rho[j]
    //                * sqrt((4.0 * Pi * Boltzmann * T(x,j)) / ((W_C / Avogadro) * n_C[k]))
    //                * pow(d_k, 2.0) * pow(Avogadro, 2.0) ;
    //     // (mol.g^(-1).s^(-1))
    //     K00_k = pow(C_k, 2.0) * kdimfwd_k / Avogadro;
    //     // (mol.g^(-1).s^(-1))
    //     K00 += K00_k;
    //     cdimer += K00_k * n_C[k];
    //     // Precursor consumption (mol/m3/s)
    //     // (x2 because 2 mol PAH -> 1 mol dimer)
    //     sootConsumption[m_precursors[k]] = 2 * K00_k * m_rho[j];
    // }
    // // Dimerization rate [1/kg/s]
    // r_dimer = K00 * Avogadro;
    // // Equivalent dimer volume [m3]
    // V_dimer = cdimer * V_C2 / K00;
}

void StFlow::sootNucleation(const double* x, size_t j){
    //----------------------------------------------------
    // Compute nucleation source terms
    //----------------------------------------------------
    // [m(11/2)/s/kg]
    beta_fm = sqrt(Pi * Boltzmann * T(x,j)
                  / (rho_soot * 2.0)) / m_rho[j];
    // Dimer collision frequency [m6/s/kg]
    beta_dimer = eps_nucl * beta_fm * sqrt(2.0 / V_dimer)
                     * pow (2.0 * pow(6.0 * V_dimer / Pi, 1.0 / 3.0), 2.0);
    // // Dimer particle density [1/m3]
    // N_dimer = sqrt((r_dimer / beta_dimer));
    // // Source term for nucleation [m3/kg/s]
    // m_qdotNucleation[j] = V_dimer * beta_dimer * pow(N_dimer, 2.0);

    // Dimer diameter [m]
    double d_dimer = 2.0 * pow(3.0 * V_dimer / (4.0 * Pi), 1.0/3.0);
    // Dimer mass [kg]
    double m_dimer = V_dimer * rho_soot;
    // Reaction rate constant [m3/kmol/s]
    // k_nuc = beta_dimer * m_rho[j] * Avogadro / 2
    double k_nuc = eps_nucl * Avogadro * pow(d_dimer, 2.0) 
                 * sqrt(4.0 * Pi * Boltzmann * T(x, j) / m_dimer);
    // Dimer concentration (QSS Assumption) [kmol/m3]
    double C_dimer = sqrt((r_dimer / Avogadro * m_rho[j]) / (2.0 * k_nuc));
    // Dimer number concentration [1/m3]
    N_dimer = C_dimer * Avogadro;
    // Nucleation reaction rate [kmol/m3/s]
    double omega_nuc = k_nuc * pow(C_dimer, 2.0);
    // Nucleation rate [m3/kg/s]
    m_qdotNucleation[j] = omega_nuc * Avogadro * V_dimer / m_rho[j] * 2.0; 

}

void StFlow::sootCondensation(const double* x, size_t j){
    //----------------------------------------------------
    // Compute condensation source term and corrects
    // nucleation source term
    //----------------------------------------------------

    double kCond = 0;
    double interval, vmoy, term1, term2, term3;
    vector<double> beta_cond(m_nsoot,0.0);
    for (size_t k = 0; k < m_nsoot; k++){
        // Condensation collision kernel[m6/s/kg]
        beta_cond[k] = eps_cond * beta_fm * sqrt(1.0/vSectMean[k]+1.0/V_dimer)
                    * pow(dSectCol[k] + pow(6.0 / Pi * V_dimer, 1.0/3.0), 2.0);
        // Collisional frequency between a dimer and any soot particle [m3/s/kg]
        kCond += beta_cond[k] * (q[k] * (vSectMax[k] - vSectMin[k]) / (vSectMean[k]));
    }
    // Correction of dimer density and nucleation source term
    N_dimer = - (kCond / (2.0 * beta_dimer))
                + sqrt((r_dimer / beta_dimer)
                + pow(kCond / (2.0 * beta_dimer),2.0));
    // Nucleation source term [m3/kg/1]
    m_qdotNucleation[j] = V_dimer * beta_dimer * pow(N_dimer,2.0);

    // Condensation source terms
    for (size_t k = 0; k < m_nsoot; k++){
        interval = (vSectMax[k] - V_dimer) - vSectMin[k];
        vmoy = ((vSectMax[k] - V_dimer) + vSectMin[k]) / 2;
        term1 = N_dimer * beta_cond[k] * q[k] / vmoy * V_dimer * interval;
        if ( k==0 ){
            term2 = 0.0;
            term3 = 0.0;
        }else{
            interval = vSectMax[k-1] - (vSectMax[k-1] - V_dimer);
            vmoy = (vSectMax[k-1] + (vSectMax[k-1] - V_dimer)) / 2;
            term2 = N_dimer * beta_cond[k-1] * q[k-1] / vmoy * V_dimer * interval;
            term3 = N_dimer * beta_cond[k-1] * q[k-1] / vmoy * vmoy * interval;
        }
        m_qdotCondensation(k,j) = term1 + term2 + term3;
        if ( k>0 ){
            m_qdotCondensation(k-1,j) -= term3;
        }
    }
}


void StFlow::sootCoagulation(const double* x, size_t j){
    //----------------------------------------------------
    // Compute coagulation source terms
    //----------------------------------------------------

    for (size_t  k = 0; k < m_nsoot; k++){
        m_qdotCoagulation(k,j) = 0.0;
    }
    // Viscosity (Cantera) [Pa.s-1] [CLEAN BUT SLOW]
    // setGas(x,j);
    // double visc = m_trans->viscosity();
    // Viscosity (Sutherland, air) [Pa.s-1] [TRASH BUT QUICK]
    double visc = Csut1 * pow(T(x,j),1.5) / (T(x,j) + Csut2);
    // Mean free path [m]
    double mfp_rod = Boltzmann * T(x,j) / (sqrt(2.0) * Pi * pow(d_gaz,2.0) * m_press);
    double mfp_tha = 66.4e-7 * (101325. / m_press) * (T(x,j) / 293.) * (1. + 110./293.) / (1 + 110. / T(x,j));
    double beta_coag_lk = 0.0;
    for (size_t  k = 0; k < m_nsoot; k++){
        for (size_t  l = k; l < m_nsoot; l++){
            if (m_collision_model == "rodrigues" || m_collision_model == "sphere"){
                //Cunnningham numbers [-]
                double Cu_k = 1 + 1.257 * 2.0 * mfp_rod / dSectCol[k];
                double Cu_l = 1 + 1.257 * 2.0 * mfp_rod / dSectCol[l];
                // Coagulation rate in continuous regime [m^6.kg^(-1).s-(-1)]
                double beta_cont_lk = (2.0 * Boltzmann * T(x,j)) / (3.0 * visc)
                                        * (dSectCol[l] + dSectCol[k])
                                        * (Cu_l / dSectCol[l] + Cu_k / dSectCol[k])
                                        / m_rho[j];
                // Coagulation rate in free molecular regime [m^6.kg^(-1).s-(-1)]
                double beta_fm_lk = eps_coag * beta_fm
                                        * sqrt(1 / vSectMean[k] + 1 / vSectMean[l])
                                        * pow(dSectCol[k] + dSectCol[l], 2.0);
                // Overall coagulation rate [m^6.kg^(-1).s-(-1)]
                beta_coag_lk = (beta_fm_lk * beta_cont_lk)
                                        / (beta_fm_lk + beta_cont_lk);
            }else if (m_collision_model == "thajudeen"){
                double reduced_mass = rho_soot * vSectMean[k] * vSectMean[l] / (vSectMean[k] + vSectMean[l]);
                double smoluchowski_radius = rNucMean[k] * (1.203 - (0.4315*(nNucMean[k] + nNucMean[l]))/(nNucMean[k]*Dfractal[k] + nNucMean[l]*Dfractal[l])) * pow(rSectCol[k]/rNucMean[k] + rSectCol[l]/rNucMean[k], 0.8806 + (0.3497*(nNucMean[k] + nNucMean[l]))/(nNucMean[k]*Dfractal[k] + nNucMean[l]*Dfractal[l]));
                double kn_diff_k = mfp_tha * Pi * rSectCol[k] / aSectCol[k];
                double kn_diff_l = mfp_tha * Pi * rSectCol[l] / aSectCol[l];
                double f_k = (6.0 * Pi * rSectCol[k] * visc) / (1 + kn_diff_k * (1.257 + 0.4 * exp(-1.1/kn_diff_k)));
                double f_l = (6.0 * Pi * rSectCol[l] * visc) / (1 + kn_diff_l * (1.257 + 0.4 * exp(-1.1/kn_diff_l)));
                double friction_coefficient = f_k * f_l / (f_k + f_l);
                double projected_area = pow(smoluchowski_radius, 2.0) * Pi;
                double kn_diff = (sqrt(Boltzmann * T(x,j) * reduced_mass) * Pi * smoluchowski_radius) / (friction_coefficient * projected_area);
                double H = (4.0 * Pi * pow(kn_diff, 2.0) + 25.836 * pow(kn_diff, 3.0) + eps_coag * sqrt(8.0 * Pi) * 11.211 * pow(kn_diff, 4.0)) / (1 + 3.502 * kn_diff + 7.211 * pow(kn_diff, 2.0) + 11.211 * pow(kn_diff, 3.0));
                beta_coag_lk = (H * friction_coefficient * pow(projected_area, 2.0)) / (reduced_mass * pow(Pi, 2.0) * smoluchowski_radius) / m_rho[j];
            }
            // Section of the particle corresponding to the coagulation of (l,k)y
            size_t m = collision_mat(k,l);
            // Loss of mass in section l due to coagulation with k
            double terml = beta_coag_lk * q[l] * q[k] / vSectMean[k]
                               * (vSectMax[l] - vSectMin[l]) * (vSectMax[k] - vSectMin[k]);
            // Loss of mass in section k due to coagulation with l
            double termk = beta_coag_lk * q[l] * q[k] / vSectMean[l]
                               * (vSectMax[l] - vSectMin[l]) * (vSectMax[k] - vSectMin[k]);
            double termm, termmplus;
            // Gain in mass in m and m+1 thanks to coagulation of (l,k)
            if (vSectMax[l] + vSectMax[k] <= vSectMax[m]){
                   termm = terml + termk;
                   termmplus = 0.0;
            }else{
                   termm = (terml + termk)
                           * (vSectMax[m] - (vSectMin[k] + vSectMin[l]))
                           / (vSectMax[k]+vSectMax[l]-vSectMin[k]-vSectMin[l]);
                   termmplus = (terml + termk)
                               * (vSectMax[k]+vSectMax[l]-vSectMax[m])
                               / (vSectMax[k]+vSectMax[l]-vSectMin[k]-vSectMin[l]);
            }
            if (m < m_nsoot-1){
                m_qdotCoagulation(m+1,j) += termmplus;
            }else{
                termm = terml + termk;
            }
            m_qdotCoagulation(m,j) += termm;
            m_qdotCoagulation(l,j) -= terml;
            m_qdotCoagulation(k,j) -= termk;
        }
    }
}

void StFlow::updateSootDiffFluxes(const double* x, size_t j0, size_t j1)
{
    //----------------------------------------------------
    // Computes soot molecular and thermal diffusion fluxes
    //----------------------------------------------------
    for (size_t j = j0; j < j1; j++) {
        // Evaluate gas at midpoints : j->j+1/2
        double Tmid = 0.5 * (T(x,j)+T(x,j+1));
        // [CLEAN BUT SLOW]
        // setGasAtMidpoint(x,j);
        // double wtm = m_thermo->meanMolecularWeight();
        // double rho = m_thermo->density();
        // double visco = m_trans->viscosity();
        // [TRASH BUT QUICK]
        double wtm = 0.5*(m_wtm[j]+m_wtm[j+1]);
        double rho = 0.5*(m_rho[j]+m_rho[j+1]);
        double visco = Csut1 * pow(Tmid,1.5) / (Tmid + Csut2);
        // Evaluate centered differencies at j+1/2
        double dT = T(x,j+1)-T(x,j);
        double dz = z(j+1) - z(j);
        // Evaluate viscosity
        for (size_t k = 0; k < m_nsoot; k++){
            // Evaluate soot mass fraction at j+1/2

            double Ysmid = 0.5 * (Ys(x,k,j)+Ys(x,k,j+1));
            // Evaluate centered differencies at j+1/2
            double dYs = Ys(x,k,j+1) - Ys(x,k,j);
            //Evaluate diffusion coefficient at j+1/2
            double Ds = 3 / (2 * rho_soot * pow(dSectMean[k], 2.0)
                                 * (1 + alphaT * Pi / 8))
                            * sqrt ((wtm*Boltzmann*Tmid)/(2*Pi*Avogadro));
            // Evaluate soot mass diffusion flux at j+1/2
            m_soot_diff(k,j) = rho * Ds * dYs/dz;
            // Evaluate soot mass soret flux at j+1/2
            if (m_do_soot_soret){
                m_soot_soret(k,j) = C_soot_soret * visco * Ysmid * dT/dz / Tmid;
            }
        }
    }
}

void StFlow::sootSurface(const double* x, size_t j){
    //----------------------------------------------------
    // Compute surface chemistry source terms
    //----------------------------------------------------

    // Computes reaction rates
    sootSurfaceInitialization(x,j);
    double term1, term2, term3;
    //---------------
    // SURFACE GROWTH
    //---------------
    if (m_do_sg){
        for (size_t k = 0; k < m_nsoot; k++){
            term1 = n_sites * alpha_surf * ksg / m_rho[j]  * vc2powervect[k] / mpower[k] * q[k]
                  * (vmaxmc2_mpower[k] - vmin_mpower[k]);
            if (k==0){
                term2 = 0.0;
                term3 = 0.0;
            }else{
                term2 = n_sites * alpha_surf * ksg / m_rho[j] * vc2powervect[k-1] / mpower[k-1] * q[k-1]
                      * (vmax_mpower[k-1] - vmaxmc2_mpower[k-1]);
                term3 = n_sites * alpha_surf * ksg / m_rho[j] * vc2powermvect[k-1] / kpower[k-1] * q[k-1]
                      * (vmax_kpower[k-1] - vmaxmc2_kpower[k-1]);
            }
            // Surface growth source term [m3/kg/s]
            m_qdotSg(k,j) = term1 + term2 + term3;
            if (k > 0){m_qdotSg(k-1,j) -= term3;}
        }
    }
    //----------
    // OXYDATION
    //----------
    if (m_do_oxidation){
        for (size_t k=m_nsoot; k --> 0;){
            term1 = alpha_surf * kox / m_rho[j] * vc2powervect[k] / mpower[k] * q[k]
                  * (vmax_mpower[k] - vminpc2_mpower[k]);
            if (k==m_nsoot-1){
                term2 = 0.0;
                term3 = 0.0;
            }else if (k==m_nsoot-2 && m_trash_section) {
                term2 = 0.0;
                term3 = 0.0;
            }else{
                term2 = alpha_surf * kox / m_rho[j] * vc2powervect[k+1] / mpower[k+1] * q[k+1]
                      * (vminpc2_mpower[k+1] - vmin_mpower[k+1]);
                term3 = alpha_surf * kox / m_rho[j] * vc2powermvect[k+1] / kpower[k+1] * q[k+1]
                      * (vminpc2_kpower[k+1] - vmin_kpower[k+1]);
            }
            // Oxidation source term [m3/kg/s]
            m_qdotOxidation(k,j) = - (term1 + term2 - term3);
            if (k < m_nsoot - 1){
                m_qdotOxidation(k+1,j) -= term3;
            }
        }
    }
}

void StFlow::sootSurfaceInitialization(const double* x, size_t j){
    // Recover species concentrations [kmol/m3]
    double conc_H = 0, conc_H2 = 0, conc_OH = 0, conc_H2O = 0, conc_C2H2 = 0, conc_O2 = 0;
    double k01f = 0, k01b = 0, k02f = 0, k02b = 0, k03f = 0, k03b = 0, k04f = 0, k04b = 0;
    double k05f = 0, k05b = 0, k06f = 0, k06bisf = 0;
    double r01f = 0, r01b = 0, r02f = 0, r02b = 0, r03f = 0, r03b = 0, r04f = 0, r04b = 0;
    double r05f = 0, r05b = 0, r06f = 0, r06bisf = 0, r07f = 0;
    double ko2, koh, kd, krev, chi, chip, chic, alpha_kazakov;
    size_t k;
    // Recover molar concentrations
    k = m_thermo->speciesIndex("H");
    if (k != npos){conc_H = Y(x,k,j) * m_rho[j] / m_wt[k];}
    k = m_thermo->speciesIndex("H2");
    if (k != npos){conc_H2 = Y(x,k,j) * m_rho[j] / m_wt[k];}
    k = m_thermo->speciesIndex("H2O");
    if (k != npos){conc_H2O = Y(x,k,j) * m_rho[j] / m_wt[k];}
    k = m_thermo->speciesIndex("C2H2");
    if (k != npos){conc_C2H2 = Y(x,k,j) * m_rho[j] / m_wt[k];}
    k = m_thermo->speciesIndex("O2");
    if (k != npos){conc_O2 = Y(x,k,j) * m_rho[j] / m_wt[k];}
    k = m_thermo->speciesIndex("OH");
    if (k != npos){conc_OH = Y(x,k,j) * m_rho[j] / m_wt[k];}
    // SI units
    double slt = 1.0e3 / GasConstant / T(x,j); // 1e3 because GasConstant in J/kmol/K
    double T_log=log(T(x,j));
    if (m_haca_model == 1 || m_haca_model == 11){ //Mauss
        // Arrhenius constants
        k01f = ak01f * exp(nk01f*T_log - ek01f*slt); // [m3/kmol/s]
        k01b = ak01b * exp(nk01b*T_log - ek01b*slt); // [m3/kmol/s]
        k02f = ak02f * exp(nk02f*T_log - ek02f*slt); // [m3/kmol/s]
        k02b = ak02b * exp(nk02b*T_log - ek02b*slt); // [m3/kmol/s]
        k03f = ak03f * exp(nk03f*T_log - ek03f*slt); // [m3/kmol/s]
        k03b = ak03b * exp(nk03b*T_log - ek03b*slt); // [1/s]
        k04f = ak04f * exp(nk04f*T_log - ek04f*slt); // [m3/kmol/s]
        k04b = ak04b * exp(nk04b*T_log - ek04b*slt); // [1/s]
        k05f = ak05f * exp(nk05f*T_log - ek05f*slt); // [1/s]
        k05b = ak05b * exp(nk05b*T_log - ek05b*slt); // [m3/kmol/s]
        k06f = ak06f * exp(nk06f*T_log - ek06f*slt); // [m3/kmol/s]
        k06bisf = k06f;
        // Partial (soot or soot radical concentration not accounted) surface chemistry reaction rates [1/s]
        r01f = k01f * conc_H; 
        r01b = k01b * conc_H2;
        r02f = k02f * conc_OH;
        r02b = k02b * conc_H2O;
        r03f = k03f * conc_H;
        r03b = k03b;
        r04f = k04f * conc_C2H2;
        r04b = k04b;
        r05f = k05f;
        r05b = k05b * conc_H;
        r06f = k06f * conc_O2;
        r06bisf = k06bisf * conc_O2;
        r07f = Avogadro * S_C2 * gamma_oh * conc_OH * 
               pow(GasConstant * T(x,j)/(2* Pi * m_wt[k]), 0.5);
        // From QSS [-]
        double fk10   = r05f / (r04b+r06bisf+r05f);
        double A_QSS = (r01f + r02f + r03b + r07f + r05b * (1.0 - fk10)) / (r01b + r02b + r03f + r04f * fk10);
        double B_QSS = r04f / ( r04b + r06bisf + r05f);
        double D_QSS = r05b / (r04b + r06f + r05f);

        chi  = 1.0;
        chip = chi * A_QSS;
        chic = chi * (A_QSS*B_QSS + D_QSS);
        // Partial reaction rates[1/s]
        koh  = r07f * chi;
        kd   = r04f * chip;
        ko2  = r06f * (chip + chic);
        krev = r04b * chic;
        // Partial surface growth and oxidation rates [1/s]
        ksg = kd - krev;
        kox = ko2 + koh;
    } else if (m_haca_model == 2){ //Blanquart
        // Arrhenius constants [m3/kmol/s]
        k01f = ak01f * exp(nk01f*T_log - ek01f*slt);
        k01b = ak01b * exp(nk01b*T_log - ek01b*slt);
        k02f = ak02f * exp(nk02f*T_log - ek02f*slt);
        k02b = ak02b * exp(nk02b*T_log - ek02b*slt);
        k03f = ak03f * exp(nk03f*T_log - ek03f*slt); // [1/s]
        k03b = ak03b * exp(nk03b*T_log - ek03b*slt);
        k04f = ak04f * exp(nk04f*T_log - ek04f*slt);
        k05f = ak05f * exp(nk05f*T_log - ek05f*slt);
        // Partial (soot or soot radical concentration not accounted) surface chemistry reaction rates [1/s]
        r01f = k01f * conc_H;
        r01b = k01b * conc_H2;
        r02f = k02f * conc_OH;
        r02b = k02b * conc_H2O;
        r03f = k03f;
        r03b = k03b * conc_H;
        r04f = k04f * conc_C2H2;
        r05f = k05f * conc_O2; 
        r06f = Avogadro * S_C2 * gamma_oh * conc_OH * 
               pow(GasConstant * T(x,j)/(2* Pi * m_wt[k]), 0.5);
        // From QSS [-]
        chi = (r01f + r02f + r03f) / (r01b + r02b + r03b + r04f + r05f); // * 1.7e15; 
        // Oxidation O2 [must be 1/s]
        ko2 = chi * r05f;
        // Oxidation OH [1/s]
        koh = r06f * chi;
        // Oxidation 
        kox = ko2 + koh;
        // Surface growth [1/s]
        ksg = chi * r04f;
    } else if (m_haca_model == 3){ //Kazakov
        // Fraction of sites available at the surface of the particle for reaction [-]
        alpha_kazakov = (1.0 + tanh(8.168e3/kazakovTad - 4.57)) / 2.0;
        // Arrhenius constants [m3/kmol/s]
        k01f = ak01f * exp(nk01f*T_log - ek01f*slt);
        k01b = ak01b * exp(nk01b*T_log - ek01b*slt);
        k02f = ak02f * exp(nk02f*T_log - ek02f*slt);
        k03f = ak03f * exp(nk03f*T_log - ek03f*slt);
        k04f = ak04f * exp(nk04f*T_log - ek04f*slt);
        // Partial (soot or soot radical concentration not accounted) surface chemistry reaction rates [1/s]
        r01f = k01f * conc_H;
        r01b = k01b * conc_H2;
        r02f = k02f * conc_H;
        r03f = k03f * conc_C2H2;
        r04f = k04f  * conc_O2;
        r05f = Avogadro * S_C2 * gamma_oh * conc_OH *
               pow(GasConstant * T(x,j)/(2* Pi * m_wt[k]), 0.5);
        // From QSS [-]
        chi = r01f / (r01b + r02f + r03f + r04f); // * 2.3e14
        // Oxidation O2 [1/s]
        ko2 = alpha_kazakov * chi * r04f;
        // Oxidation OH [1/s]
        koh = r05f * chi;
        //Oxidation [1/s]
        kox = ko2 + koh;
        // Surface growth [1/s]
        ksg = alpha_kazakov * chi * r03f;
    }
    //-------------------------
    // RETROACTION ON GAS PHASE
    //-------------------------
    if (m_do_retroaction){
        double deltaQ = 0.0, soot_surface = 0.0;
        for (size_t k = 0; k < m_nsoot; k++){
            deltaQ += vc2powervect[k] / mpower[k] * q[k] * (vmax_mpower[k]-vmin_mpower[k]);
        }
        //Overall active sites concentration [kmol/m3]
        soot_surface = alpha_surf * n_sites * deltaQ / (V_C2 * Avogadro); // n_sites added here
        // Surface chemistry reaction rates [kmol/s/m3]
        double w01f = 0, w01b = 0, w02f = 0, w02b = 0, w03f = 0, w03b = 0;
        double w04f = 0, w04b = 0, w05f = 0, w05b = 0, w06f = 0, w06bisf = 0;
        double w07f = 0;
        if (m_do_sg || m_do_oxidation){
            // Reaction rates for reactions occuring on both oxidation & surface growth
            if (m_haca_model == 1){
                w01f = r01f * soot_surface * chi;
                w01b = r01b * soot_surface * chip;
                w02f = r02f * soot_surface * chi;
                w02b = r02b * soot_surface * chip;
                w03f = r03f * soot_surface * chip;
                w03b = r03b * soot_surface * chi;
                w05f = r05f * soot_surface * chic;
                w05b = r05b * soot_surface * chi;
                w06bisf = r06bisf * soot_surface * chic;
            }else if (m_haca_model == 2){
                w01f = r01f * soot_surface * chi;
                w01b = r01b * soot_surface * chi;
                w02f = r02f * soot_surface * chi;
                w02b = r02b * soot_surface * chi;
                w03f = r03f * soot_surface * chi;
                w03b = r03b * soot_surface * chi;
            }else if (m_haca_model == 3){
                w01f = r01f * soot_surface * chi;
                w01b = r01b * soot_surface * chi;
                w02f = r02f * soot_surface * chi;
            }
        }
        // Reaction rates for reactions occuring on surfacce growth only
        if (m_do_sg){
            if (m_haca_model == 1){
                w04f = r04f * soot_surface * chip;
                w04b = r04b * soot_surface * chic;
            }else if (m_haca_model == 2){
                w04f = r04f * soot_surface * chi;
            }else if (m_haca_model == 3){
                w03f = r03f * soot_surface * chi;
            }
        }
        // Reaction rates for reactions occuring on oxidation only
        if (m_do_oxidation){
            if (m_haca_model == 1){
                // Reaction rates [kmol/m3/s]
                w06f = r06f * soot_surface * chip;
                w07f = r07f * soot_surface * chi;
            }else if (m_haca_model == 2){
                w05f = r05f * soot_surface * chi;
                w06f = r06f * soot_surface * chi;
            }else if (m_haca_model == 3){
                w04f = r04f * soot_surface * chi;
                w05f = r05f * soot_surface * chi;
            }
        }
        
        // Gaseous species consumption [kmol/s/m3]
        if (m_haca_model == 1){
            k = m_thermo->speciesIndex("H");
            if (k != npos){sootConsumption[k] = - (w01f - w01b) - (w03f - w03b) + (w05f - w05b);}
            k = m_thermo->speciesIndex("C2H2");
            if (k != npos){sootConsumption[k] = (w04f - w04b);}
            k = m_thermo->speciesIndex("H2");
            if (k != npos){sootConsumption[k] =  (w01f - w01b);}
            k = m_thermo->speciesIndex("O2");
            if (k != npos){sootConsumption[k] = - w06f - w06bisf;}
            k = m_thermo->speciesIndex("CO");
            if (k != npos){sootConsumption[k] = 2.0 * w06f;}
            k = m_thermo->speciesIndex("OH");
            if (k != npos){sootConsumption[k] = - (w02f - w02b) - w07f;}
            k = m_thermo->speciesIndex("H2O");
            if (k != npos){sootConsumption[k] = + (w02f - w02b);}
            k = m_thermo->speciesIndex("HCO");
            if (k != npos){sootConsumption[k] = 2.0 * w06bisf + w07f;}
            k = m_thermo->speciesIndex("CH");
            if (k != npos){sootConsumption[k] = w07f;}
        }else if (m_haca_model == 2){
            k = m_thermo->speciesIndex("H");
            if (k != npos){sootConsumption[k] = -(w01f-w01b) + (w03f-w03b) + w04f;}
            k = m_thermo->speciesIndex("C2H2");
            if (k != npos){sootConsumption[k] = -w04f;}
            k = m_thermo->speciesIndex("H2");
            if (k != npos){sootConsumption[k] =  (w01f-w01b) + w05f + w06f;}
            k = m_thermo->speciesIndex("O2");
            if (k != npos){sootConsumption[k] = -w05f;}
            k = m_thermo->speciesIndex("CO");
            if (k != npos){sootConsumption[k] = 2.0 * w05f + w06f;}
            k = m_thermo->speciesIndex("OH");
            if (k != npos){sootConsumption[k] = -(w02f-w02b) - w06f;}
            k = m_thermo->speciesIndex("H2O");
            if (k != npos){sootConsumption[k] = (w02f-w02b);}        
        }else if (m_haca_model == 3){
            k = m_thermo->speciesIndex("H");
            if (k != npos){sootConsumption[k] = -(w01f-w01b) - w02f + w03f;}
            k = m_thermo->speciesIndex("C2H2");
            if (k != npos){sootConsumption[k] = -w03f;}
            k = m_thermo->speciesIndex("H2");
            if (k != npos){sootConsumption[k] =  (w01f-w01b) + w05f;}
            k = m_thermo->speciesIndex("O2");
            if (k != npos){sootConsumption[k] = -w04f;}
            k = m_thermo->speciesIndex("CO");
            if (k != npos){sootConsumption[k] = 2.0 * w04f + w05f;}
            k = m_thermo->speciesIndex("OH");
            if (k != npos){sootConsumption[k] = - w05f;}
        }
    }
}

void StFlow::getDistributionOrdre0(const double* x, size_t j){
    for (size_t k = 0; k < m_nsoot; k++){
        // Soot volume fraction density [1/m3]
        q[k] = m_rho[j] * Ys(x,k,j) / (rho_soot * (vSectMax[k] - vSectMin[k]));
    }
}
//-----------------------------//
// END OF SOOT RELATED METHODS //
//-----------------------------//

void StFlow::evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
{
    size_t j = m_points - 1;
    AVBPcompute_local_thick(x,j);
    double thick_prev = avbp_thick[j-1];

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    rsd[index(c_offset_V,j)] = V(x,j);
    diag[index(c_offset_V,j)] = 0;
    double sum = 0.0;
    // set residual of poisson's equ to zero
    rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    if (m_nsoot > 0){
        for (size_t k = 0; k < m_nsoot; k++){
            if (m_do_retroaction){sum += Ys(x,k,j);}
            rsd[index(c_offset_S + k, j)] =
            m_soot_diff(k,j-1) + rho_u(x,j)*Ys(x,k,j);
        }
    }
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1)*thick_prev + rho_u(x,j)*Y(x,k,j);
// no thick         rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
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
