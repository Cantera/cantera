/**
 *  @file MultiTransport.cpp
 *  Implementation file for class MultiTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/MultiTransport.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

///////////////////// helper functions /////////////////////////

/**
 * The Parker temperature correction to the rotational collision number.
 *
 * @param tr Reduced temperature \f$ \epsilon/kT \f$
 * @param sqtr square root of tr.
 */
doublereal Frot(doublereal tr, doublereal sqtr)
{
    const doublereal c1 = 0.5*sqrt(Pi)*Pi;
    const doublereal c2 = 0.25*Pi*Pi + 2.0;
    const doublereal c3 = sqrt(Pi)*Pi;
    return 1.0 + c1*sqtr + c2*tr + c3*sqtr*tr;
}

//////////////////// class MultiTransport methods //////////////

MultiTransport::MultiTransport(thermo_t* thermo)
    : GasTransport(thermo)
{
}

void MultiTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    GasTransport::init(thermo, mode, log_level);

    // the L matrix
    m_Lmatrix.resize(3*m_nsp, 3*m_nsp);
    m_a.resize(3*m_nsp, 1.0);
    m_b.resize(3*m_nsp, 0.0);
    m_aa.resize(m_nsp, m_nsp, 0.0);
    m_molefracs_last.resize(m_nsp, -1.0);
    m_frot_298.resize(m_nsp);
    m_rotrelax.resize(m_nsp);
    m_cinternal.resize(m_nsp);
    m_om22.resize(m_nsp, m_nsp);
    m_astar.resize(m_nsp, m_nsp);
    m_bstar.resize(m_nsp, m_nsp);
    m_cstar.resize(m_nsp, m_nsp);

    // set flags all false
    m_abc_ok = false;
    m_l0000_ok = false;
    m_lmatrix_soln_ok = false;
    m_thermal_tlast = 0.0;

    // some work space
    m_spwork1.resize(m_nsp);
    m_spwork2.resize(m_nsp);
    m_spwork3.resize(m_nsp);

    // precompute and store log(epsilon_ij/k_B)
    m_log_eps_k.resize(m_nsp, m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            m_log_eps_k(i,j) = log(m_epsilon(i,j)/Boltzmann);
            m_log_eps_k(j,i) = m_log_eps_k(i,j);
        }
    }

    // precompute and store constant parts of the Parker rotational
    // collision number temperature correction
    const doublereal sq298 = sqrt(298.0);
    const doublereal kb298 = Boltzmann * 298.0;
    m_sqrt_eps_k.resize(m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        m_sqrt_eps_k[k] = sqrt(m_eps[k]/Boltzmann);
        m_frot_298[k] = Frot(m_eps[k]/kb298, m_sqrt_eps_k[k]/sq298);
    }
}

doublereal MultiTransport::thermalConductivity()
{
    solveLMatrixEquation();
    doublereal sum = 0.0;
    for (size_t k = 0; k  < 2*m_nsp; k++) {
        sum += m_b[k + m_nsp] * m_a[k + m_nsp];
    }
    return -4.0*sum;
}

void MultiTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    solveLMatrixEquation();
    const doublereal c = 1.6/GasConstant;
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = c * m_mw[k] * m_molefracs[k] * m_a[k];
    }
}

void MultiTransport::solveLMatrixEquation()
{
    // if T has changed, update the temperature-dependent properties.
    updateThermal_T();
    update_C();
    if (m_lmatrix_soln_ok) {
        return;
    }

    // Copy the mole fractions twice into the last two blocks of the right-hand-
    // side vector m_b. The first block of m_b was set to zero when it was
    // created, and is not modified so doesn't need to be reset to zero.
    for (size_t k = 0; k < m_nsp; k++) {
        m_b[k] = 0.0;
        m_b[k + m_nsp] = m_molefracs[k];
        m_b[k + 2*m_nsp] = m_molefracs[k];
    }

    // Set the right-hand side vector to zero in the 3rd block for all species
    // with no internal energy modes.  The corresponding third-block rows and
    // columns will be set to zero, except on the diagonal of L01,01, where they
    // are set to 1.0. This has the effect of eliminating these equations from
    // the system, since the equation becomes: m_a[2*m_nsp + k] = 0.0.

    // Note that this differs from the Chemkin procedure, where all *monatomic*
    // species are excluded. Since monatomic radicals can have non-zero internal
    // heat capacities due to electronic excitation, they should be retained.
    for (size_t k = 0; k < m_nsp; k++) {
        if (!hasInternalModes(k)) {
            m_b[2*m_nsp + k] = 0.0;
        }
    }

    // evaluate the submatrices of the L matrix
    m_Lmatrix.resize(3*m_nsp, 3*m_nsp, 0.0);

    //! Evaluate the upper-left block of the L matrix.
    eval_L0000(m_molefracs.data());
    eval_L0010(m_molefracs.data());
    eval_L0001();
    eval_L1000();
    eval_L1010(m_molefracs.data());
    eval_L1001(m_molefracs.data());
    eval_L0100();
    eval_L0110();
    eval_L0101(m_molefracs.data());

    // Solve it using GMRES or LU decomposition. The last solution in m_a should
    // provide a good starting guess, so convergence should be fast.
    m_a = m_b;
    solve(m_Lmatrix, m_a.data());
    m_lmatrix_soln_ok = true;
    m_molefracs_last = m_molefracs;
    // L matrix is overwritten with LU decomposition
    m_l0000_ok = false;
}

void MultiTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                      size_t ldx, const doublereal* const grad_X,
                                      size_t ldf, doublereal* const fluxes)
{
    // update the binary diffusion coefficients if necessary
    update_T();
    updateDiff_T();

    // If any component of grad_T is non-zero, then get the
    // thermal diffusion coefficients
    bool addThermalDiffusion = false;
    for (size_t i = 0; i < ndim; i++) {
        if (grad_T[i] != 0.0) {
            addThermalDiffusion = true;
        }
    }
    if (addThermalDiffusion) {
        getThermalDiffCoeffs(m_spwork.data());
    }

    const doublereal* y = m_thermo->massFractions();
    doublereal rho = m_thermo->density();

    for (size_t i = 0; i < m_nsp; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            m_aa(i,j) = m_molefracs[j]*m_molefracs[i]/m_bdiff(i,j);
            sum += m_aa(i,j);
        }
        m_aa(i,i) -= sum;
    }

    // enforce the condition \sum Y_k V_k = 0. This is done by replacing
    // the flux equation with the largest gradx component in the first
    // coordinate direction with the flux balance condition.
    size_t jmax = 0;
    doublereal gradmax = -1.0;
    for (size_t j = 0; j < m_nsp; j++) {
        if (fabs(grad_X[j]) > gradmax) {
            gradmax = fabs(grad_X[j]);
            jmax = j;
        }
    }

    // set the matrix elements in this row to the mass fractions,
    // and set the entry in gradx to zero
    for (size_t j = 0; j < m_nsp; j++) {
        m_aa(jmax,j) = y[j];
    }
    vector_fp gsave(ndim), grx(ldx*m_nsp);
    for (size_t n = 0; n < ldx*ndim; n++) {
        grx[n] = grad_X[n];
    }

    // copy grad_X to fluxes
    for (size_t n = 0; n < ndim; n++) {
        const double* gx = grad_X + ldx*n;
        copy(gx, gx + m_nsp, fluxes + ldf*n);
        fluxes[jmax + n*ldf] = 0.0;
    }

    // solve the equations
    solve(m_aa, fluxes, ndim, ldf);
    doublereal pp = pressure_ig();

    // multiply diffusion velocities by rho * V to create mass fluxes, and
    // restore the gradx elements that were modified
    for (size_t n = 0; n < ndim; n++) {
        size_t offset = n*ldf;
        for (size_t i = 0; i < m_nsp; i++) {
            fluxes[i + offset] *= rho * y[i] / pp;
        }
    }

    // thermal diffusion
    if (addThermalDiffusion) {
        for (size_t n = 0; n < ndim; n++) {
            size_t offset = n*ldf;
            doublereal grad_logt = grad_T[n]/m_temp;
            for (size_t i = 0; i < m_nsp; i++) {
                fluxes[i + offset] -= m_spwork[i]*grad_logt;
            }
        }
    }
}

void MultiTransport::getMassFluxes(const doublereal* state1, const doublereal* state2, doublereal delta,
                                   doublereal* fluxes)
{
    double* x1 = m_spwork1.data();
    double* x2 = m_spwork2.data();
    double* x3 = m_spwork3.data();
    size_t nsp = m_thermo->nSpecies();
    m_thermo->restoreState(nsp+2, state1);
    double p1 = m_thermo->pressure();
    double t1 = state1[0];
    m_thermo->getMoleFractions(x1);

    m_thermo->restoreState(nsp+2, state2);
    double p2 = m_thermo->pressure();
    double t2 = state2[0];
    m_thermo->getMoleFractions(x2);

    double p = 0.5*(p1 + p2);
    double t = 0.5*(state1[0] + state2[0]);

    for (size_t n = 0; n < nsp; n++) {
        x3[n] = 0.5*(x1[n] + x2[n]);
    }
    m_thermo->setState_TPX(t, p, x3);
    m_thermo->getMoleFractions(m_molefracs.data());

    // update the binary diffusion coefficients if necessary
    update_T();
    updateDiff_T();

    // If there is a temperature gradient, then get the
    // thermal diffusion coefficients
    bool addThermalDiffusion = false;
    if (state1[0] != state2[0]) {
        addThermalDiffusion = true;
        getThermalDiffCoeffs(m_spwork.data());
    }

    const doublereal* y = m_thermo->massFractions();
    doublereal rho = m_thermo->density();
    for (size_t i = 0; i < m_nsp; i++) {
        doublereal sum = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            m_aa(i,j) = m_molefracs[j]*m_molefracs[i]/m_bdiff(i,j);
            sum += m_aa(i,j);
        }
        m_aa(i,i) -= sum;
    }

    // enforce the condition \sum Y_k V_k = 0. This is done by replacing the
    // flux equation with the largest gradx component with the flux balance
    // condition.
    size_t jmax = 0;
    doublereal gradmax = -1.0;
    for (size_t j = 0; j < m_nsp; j++) {
        if (fabs(x2[j] - x1[j]) > gradmax) {
            gradmax = fabs(x1[j] - x2[j]);
            jmax = j;
        }
    }

    // set the matrix elements in this row to the mass fractions,
    // and set the entry in gradx to zero
    for (size_t j = 0; j < m_nsp; j++) {
        m_aa(jmax,j) = y[j];
        fluxes[j] = x2[j] - x1[j];
    }
    fluxes[jmax] = 0.0;

    // Solve the equations
    solve(m_aa, fluxes);

    doublereal pp = pressure_ig();
    // multiply diffusion velocities by rho * Y_k to create
    // mass fluxes, and divide by pressure
    for (size_t i = 0; i < m_nsp; i++) {
        fluxes[i] *= rho * y[i] / pp;
    }

    // thermal diffusion
    if (addThermalDiffusion) {
        doublereal grad_logt = (t2 - t1)/m_temp;
        for (size_t i = 0; i < m_nsp; i++) {
            fluxes[i] -= m_spwork[i]*grad_logt;
        }
    }
}

void MultiTransport::getMolarFluxes(const doublereal* const state1,
                                    const doublereal* const state2,
                                    const doublereal delta,
                                    doublereal* const fluxes)
{
    getMassFluxes(state1, state2, delta, fluxes);
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        fluxes[k] /= m_mw[k];
    }
}

void MultiTransport::getMultiDiffCoeffs(const size_t ld, doublereal* const d)
{
    doublereal p = pressure_ig();

    // update the mole fractions
    update_C();

    // update the binary diffusion coefficients
    update_T();
    updateThermal_T();

    // evaluate L0000 if the temperature or concentrations have
    // changed since it was last evaluated.
    if (!m_l0000_ok) {
        eval_L0000(m_molefracs.data());
    }

    // invert L00,00
    int ierr = invert(m_Lmatrix, m_nsp);
    if (ierr != 0) {
        throw CanteraError("MultiTransport::getMultiDiffCoeffs",
                           "invert returned ierr = {}", ierr);
    }
    m_l0000_ok = false; // matrix is overwritten by inverse
    m_lmatrix_soln_ok = false;

    doublereal prefactor = 16.0 * m_temp
                           * m_thermo->meanMolecularWeight()/(25.0 * p);
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            double c = prefactor/m_mw[j];
            d[ld*j + i] = c*m_molefracs[i]*
                          (m_Lmatrix(i,j) - m_Lmatrix(i,i));
        }
    }
}

void MultiTransport::update_T()
{
    if (m_temp == m_thermo->temperature() && m_nsp == m_thermo->nSpecies()) {
        return;
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be
    // redone, and the L matrix reevaluated.
    m_abc_ok = false;
    m_lmatrix_soln_ok = false;
    m_l0000_ok = false;
}

void MultiTransport::update_C()
{
    // Update the local mole fraction array
    m_thermo->getMoleFractions(m_molefracs.data());

    for (size_t k = 0; k < m_nsp; k++) {
        // add an offset to avoid a pure species condition
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
        if (m_molefracs[k] != m_molefracs_last[k]) {
            // If any mole fractions have changed, signal that concentration-
            // dependent quantities will need to be recomputed before use.
            m_l0000_ok = false;
            m_lmatrix_soln_ok = false;
        }
    }
}

void MultiTransport::updateThermal_T()
{
    if (m_thermal_tlast == m_thermo->temperature()) {
        return;
    }
    // we need species viscosities and binary diffusion coefficients
    updateSpeciesViscosities();
    updateDiff_T();

    // evaluate polynomial fits for A*, B*, C*
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            double z = m_logt - m_log_eps_k(i,j);
            int ipoly = m_poly[i][j];
            if (m_mode == CK_Mode) {
                m_om22(i,j) = poly6(z, m_omega22_poly[ipoly].data());
                m_astar(i,j) = poly6(z, m_astar_poly[ipoly].data());
                m_bstar(i,j) = poly6(z, m_bstar_poly[ipoly].data());
                m_cstar(i,j) = poly6(z, m_cstar_poly[ipoly].data());
            } else {
                m_om22(i,j) = poly8(z, m_omega22_poly[ipoly].data());
                m_astar(i,j) = poly8(z, m_astar_poly[ipoly].data());
                m_bstar(i,j) = poly8(z, m_bstar_poly[ipoly].data());
                m_cstar(i,j) = poly8(z, m_cstar_poly[ipoly].data());
            }
            m_om22(j,i) = m_om22(i,j);
            m_astar(j,i) = m_astar(i,j);
            m_bstar(j,i) = m_bstar(i,j);
            m_cstar(j,i) = m_cstar(i,j);
        }
    }
    m_abc_ok = true;

    // evaluate the temperature-dependent rotational relaxation rate
    for (size_t k = 0; k < m_nsp; k++) {
        double tr = m_eps[k]/ m_kbt;
        double sqtr = m_sqrt_eps_k[k] / m_sqrt_t;
        m_rotrelax[k] = std::max(1.0,m_zrot[k]) * m_frot_298[k]/Frot(tr, sqtr);
    }

    doublereal c = 1.2*GasConstant*m_temp;
    for (size_t k = 0; k < m_nsp; k++) {
        m_bdiff(k,k) = c * m_visc[k] * m_astar(k,k)/m_mw[k];
    }

    // Calculate the internal heat capacities by subtracting off the translational contributions
    /*
     *  HKM Exploratory comment:
     *       The translational component is 1.5
     *       The rotational component is 1.0 for a linear molecule and 1.5 for a nonlinear molecule
     *           and zero for a monatomic.
     *       Chemkin has traditionally subtracted 1.5 here (SAND86-8246).
     *       The original Dixon-Lewis paper subtracted 1.5 here.
     */
    vector_fp cp(m_thermo->nSpecies());
    m_thermo->getCp_R_ref(&cp[0]);
    for (size_t k = 0; k < m_nsp; k++) {
        m_cinternal[k] = cp[k] - 2.5;
    }
    m_thermal_tlast = m_thermo->temperature();
}

//! Constant to compare dimensionless heat capacities against zero
static const doublereal Min_C_Internal = 0.001;

bool MultiTransport::hasInternalModes(size_t j)
{
    return (m_cinternal[j] > Min_C_Internal);
}

void MultiTransport::eval_L0000(const doublereal* const x)
{
    doublereal prefactor = 16.0*m_temp/25.0;
    doublereal sum;
    for (size_t i = 0; i < m_nsp; i++) {
        // subtract-off the k=i term to account for the first delta
        // function in Eq. (12.121)
        sum = -x[i]/m_bdiff(i,i);
        for (size_t k = 0; k < m_nsp; k++) {
            sum += x[k]/m_bdiff(i,k);
        }

        sum /= m_mw[i];
        for (size_t j = 0; j != m_nsp; ++j) {
            m_Lmatrix(i,j) = prefactor * x[j]
                             * (m_mw[j] * sum + x[i]/m_bdiff(i,j));
        }
        // diagonal term is zero
        m_Lmatrix(i,i) = 0.0;
    }
}

void MultiTransport::eval_L0010(const doublereal* const x)
{
    doublereal prefactor = 1.6*m_temp;
    for (size_t j = 0; j < m_nsp; j++) {
        double xj = x[j];
        double wj = m_mw[j];
        double sum = 0.0;
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i,j + m_nsp) = - prefactor * x[i] * xj * m_mw[i] *
                                     (1.2 * m_cstar(j,i) - 1.0) /
                                     ((wj + m_mw[i]) * m_bdiff(j,i));

            // the next term is independent of "j";
            // need to do it for the "j,j" term
            sum -= m_Lmatrix(i,j+m_nsp);
        }
        m_Lmatrix(j,j+m_nsp) += sum;
    }
}

void MultiTransport::eval_L1000()
{
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+m_nsp,j) = m_Lmatrix(j,i+m_nsp);
        }
    }
}

void MultiTransport::eval_L1010(const doublereal* x)
{
    const doublereal fiveover3pi = 5.0/(3.0*Pi);
    doublereal prefactor = (16.0*m_temp)/25.0;

    for (size_t j = 0; j < m_nsp; j++) {
        // get constant terms that depend on just species "j"
        double constant1 = prefactor*x[j];
        double wjsq = m_mw[j]*m_mw[j];
        double constant2 = 13.75*wjsq;
        double constant3 = m_crot[j]/m_rotrelax[j];
        double constant4 = 7.5*wjsq;
        double fourmj = 4.0*m_mw[j];
        double threemjsq = 3.0*m_mw[j]*m_mw[j];
        double sum = 0.0;
        for (size_t i = 0; i < m_nsp; i++) {
            double sumwij = m_mw[i] + m_mw[j];
            double term1 = m_bdiff(i,j) * sumwij*sumwij;
            double term2 = fourmj*m_astar(i,j)*(1.0 + fiveover3pi*
                (constant3 + (m_crot[i]/m_rotrelax[i]))); //  see Eq. (12.125)

            m_Lmatrix(i+m_nsp,j+m_nsp) = constant1*x[i]*m_mw[i] /(m_mw[j]*term1) *
                                         (constant2 - threemjsq*m_bstar(i,j)
                                          - term2*m_mw[j]);

            sum += x[i] /(term1) *
                   (constant4 + m_mw[i]*m_mw[i]*
                    (6.25 - 3.0*m_bstar(i,j)) + term2*m_mw[i]);
        }

        m_Lmatrix(j+m_nsp,j+m_nsp) -= sum*constant1;
    }
}

void MultiTransport::eval_L1001(const doublereal* x)
{
    doublereal prefactor = 32.00*m_temp/(5.00*Pi);
    for (size_t j = 0; j < m_nsp; j++) {
        // collect terms that depend only on "j"
        if (hasInternalModes(j)) {
            double constant = prefactor*m_mw[j]*x[j]*m_crot[j]/(m_cinternal[j]*m_rotrelax[j]);
            double sum = 0.0;
            for (size_t i = 0; i < m_nsp; i++) {
                // see Eq. (12.127)
                m_Lmatrix(i+m_nsp,j+2*m_nsp) = constant * m_astar(j,i) * x[i] /
                                          ((m_mw[j] + m_mw[i]) * m_bdiff(j,i));
                sum += m_Lmatrix(i+m_nsp,j+2*m_nsp);
            }
            m_Lmatrix(j+m_nsp,j+2*m_nsp) += sum;
        } else {
            for (size_t i = 0; i < m_nsp; i++) {
                m_Lmatrix(i+m_nsp,j+2*m_nsp) = 0.0;
            }
        }
    }
}

void MultiTransport::eval_L0001()
{
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i,j+2*m_nsp) = 0.0;
        }
    }
}

void MultiTransport::eval_L0100()
{
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+2*m_nsp,j) = 0.0; // see Eq. (12.123)
        }
    }
}

void MultiTransport::eval_L0110()
{
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+2*m_nsp,j+m_nsp) = m_Lmatrix(j+m_nsp,i+2*m_nsp); // see Eq. (12.123)
        }
    }
}

void MultiTransport::eval_L0101(const doublereal* x)
{
    for (size_t i = 0; i < m_nsp; i++) {
        if (hasInternalModes(i)) {
            // collect terms that depend only on "i"
            double constant1 = 4*m_temp*x[i]/m_cinternal[i];
            double constant2 = 12*m_mw[i]*m_crot[i] /
                               (5*Pi*m_cinternal[i]*m_rotrelax[i]);
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                // see Eq. (12.131)
                double diff_int = m_bdiff(i,k);
                m_Lmatrix(k+2*m_nsp,i+2*m_nsp) = 0.0;
                sum += x[k]/diff_int;
                if (k != i) {
                    sum += x[k]*m_astar(i,k)*constant2 / (m_mw[k]*diff_int);
                }
            }
            // see Eq. (12.130)
            m_Lmatrix(i+2*m_nsp,i+2*m_nsp) =
                - 8/Pi*m_mw[i]*x[i]*x[i]*m_crot[i] /
                (m_cinternal[i]*m_cinternal[i]*GasConstant*m_visc[i]*m_rotrelax[i])
                - constant1*sum;
        } else {
            for (size_t k = 0; k < m_nsp; k++) {
                m_Lmatrix(i+2*m_nsp,i+2*m_nsp) = 1.0;
            }
        }
    }
}

}
