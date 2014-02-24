/**
 *  @file MultiTransport.cpp
 *  Implementation file for class MultiTransport
 */
/*
 *  Copyright 2001 California Institute of Technology
 *  See file License.txt for licensing information
 */

#include "cantera/transport/MultiTransport.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "L_matrix.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

///////////////////// helper functions /////////////////////////

/**
 *  The Parker temperature correction to the rotational collision number.
 *
 *  @param tr Reduced temperature \f$ \epsilon/kT \f$
 *  @param sqtr square root of tr.
 */
inline doublereal Frot(doublereal tr, doublereal sqtr)
{
    const doublereal c1 = 0.5*SqrtPi*Pi;
    const doublereal c2 = 0.25*Pi*Pi + 2.0;
    const doublereal c3 = SqrtPi*Pi;
    return 1.0 + c1*sqtr + c2*tr + c3*sqtr*tr;
}

//////////////////// class MultiTransport methods //////////////

MultiTransport::MultiTransport(thermo_t* thermo)
    : GasTransport(thermo)
{
}

bool MultiTransport::initGas(GasTransportParams& tr)
{
    GasTransport::initGas(tr);

    // copy polynomials and parameters into local storage
    m_poly       = tr.poly;
    m_astar_poly = tr.astar_poly;
    m_bstar_poly = tr.bstar_poly;
    m_cstar_poly = tr.cstar_poly;
    m_om22_poly  = tr.omega22_poly;
    m_zrot       = tr.zrot;
    m_crot       = tr.crot;
    m_eps        = tr.eps;
    m_sigma      = tr.sigma;
    m_alpha      = tr.alpha;
    m_dipole     = tr.dipole;
    m_zrot       = tr.zrot;

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
    //        int j;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            m_log_eps_k(i,j) = log(tr.epsilon(i,j)/Boltzmann);
            m_log_eps_k(j,i) = m_log_eps_k(i,j);
        }
    }

    // precompute and store constant parts of the Parker rotational
    // collision number temperature correction
    const doublereal sq298 = sqrt(298.0);
    const doublereal kb298 = Boltzmann * 298.0;
    m_sqrt_eps_k.resize(m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        m_sqrt_eps_k[k] = sqrt(tr.eps[k]/Boltzmann);
        m_frot_298[k] = Frot(tr.eps[k]/kb298,
                             m_sqrt_eps_k[k]/sq298);
    }

    return true;
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

    // Copy the mole fractions twice into the last two blocks of
    // the right-hand-side vector m_b. The first block of m_b was
    // set to zero when it was created, and is not modified so
    // doesn't need to be reset to zero.
    for (size_t k = 0; k < m_nsp; k++) {
        m_b[k] = 0.0;
        m_b[k + m_nsp] = m_molefracs[k];
        m_b[k + 2*m_nsp] = m_molefracs[k];
    }

    // Set the right-hand side vector to zero in the 3rd block for
    // all species with no internal energy modes.  The
    // corresponding third-block rows and columns will be set to
    // zero, except on the diagonal of L01,01, where they are set
    // to 1.0. This has the effect of eliminating these equations
    // from the system, since the equation becomes: m_a[2*m_nsp +
    // k] = 0.0.

    // Note that this differs from the Chemkin procedure, where
    // all *monatomic* species are excluded. Since monatomic
    // radicals can have non-zero internal heat capacities due to
    // electronic excitation, they should be retained.

    for (size_t k = 0; k < m_nsp; k++) {
        if (!hasInternalModes(k)) {
            m_b[2*m_nsp + k] = 0.0;
        }
    }

    // evaluate the submatrices of the L matrix
    m_Lmatrix.resize(3*m_nsp, 3*m_nsp, 0.0);

    //! Evaluate the upper-left block of the L matrix.
    eval_L0000(DATA_PTR(m_molefracs));
    eval_L0010(DATA_PTR(m_molefracs));
    eval_L0001();
    eval_L1000();
    eval_L1010(DATA_PTR(m_molefracs));
    eval_L1001(DATA_PTR(m_molefracs));
    eval_L0100();
    eval_L0110();
    eval_L0101(DATA_PTR(m_molefracs));

    // Solve it using GMRES or LU decomposition. The last solution
    // in m_a should provide a good starting guess, so convergence
    // should be fast.

    //if (m_gmres) {
    //    gmres(m_mgmres, 3*m_nsp, m_Lmatrix, m_b.begin(),
    //        m_a.begin(), m_eps_gmres);
    //    m_lmatrix_soln_ok = true;
    //    m_l0000_ok = true;            // L matrix not modified by GMRES
    //}
    //else {
    copy(m_b.begin(), m_b.end(), m_a.begin());
    try {
        solve(m_Lmatrix, DATA_PTR(m_a));
    } catch (CanteraError& err) {
        err.save();
        //if (info != 0) {
        throw CanteraError("MultiTransport::solveLMatrixEquation",
                           "error in solving L matrix.");
    }
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
        getThermalDiffCoeffs(DATA_PTR(m_spwork));
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
    //for (n = 0; n < ndim; n++) {
    //    gsave[n] = grad_X[jmax + n*ldx];   // save the input mole frac gradient
    //grad_X[jmax + n*ldx] = 0.0;
    //    grx[jmax + n*ldx] = 0.0;
    // }

    // copy grad_X to fluxes
    const doublereal* gx;
    for (size_t n = 0; n < ndim; n++) {
        gx = grad_X + ldx*n;
        copy(gx, gx + m_nsp, fluxes + ldf*n);
        fluxes[jmax + n*ldf] = 0.0;
    }

    // use LAPACK to solve the equations
    int info=0;
    ct_dgetrf(static_cast<int>(m_aa.nRows()),
              static_cast<int>(m_aa.nColumns()), m_aa.ptrColumn(0),
              static_cast<int>(m_aa.nRows()),
              &m_aa.ipiv()[0], info);
    if (info == 0) {
        ct_dgetrs(ctlapack::NoTranspose,
                  static_cast<int>(m_aa.nRows()), ndim,
                  m_aa.ptrColumn(0), static_cast<int>(m_aa.nRows()),
                  &m_aa.ipiv()[0], fluxes, ldf, info);
        if (info != 0) {
            info += 100;
        }
    } else
        throw CanteraError("MultiTransport::getSpeciesFluxes",
                           "Error in DGETRF");
    if (info > 50)
        throw CanteraError("MultiTransport::getSpeciesFluxes",
                           "Error in DGETRS");

    size_t offset;
    doublereal pp = pressure_ig();

    // multiply diffusion velocities by rho * V to create
    // mass fluxes, and restore the gradx elements that were
    // modified
    for (size_t n = 0; n < ndim; n++) {
        offset = n*ldf;
        for (size_t i = 0; i < m_nsp; i++) {
            fluxes[i + offset] *= rho * y[i] / pp;
        }
        //grad_X[jmax + n*ldx] = gsave[n];
    }

    // thermal diffusion
    if (addThermalDiffusion) {
        for (size_t n = 0; n < ndim; n++) {
            offset = n*ldf;
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
    double* x1 = DATA_PTR(m_spwork1);
    double* x2 = DATA_PTR(m_spwork2);
    double* x3 = DATA_PTR(m_spwork3);
    size_t n, nsp = m_thermo->nSpecies();
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

    for (n = 0; n < nsp; n++) {
        x3[n] = 0.5*(x1[n] + x2[n]);
    }
    m_thermo->setState_TPX(t, p, x3);
    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // update the binary diffusion coefficients if necessary
    update_T();
    updateDiff_T();

    // If there is a temperature gradient, then get the
    // thermal diffusion coefficients

    bool addThermalDiffusion = false;
    if (state1[0] != state2[0]) {
        addThermalDiffusion = true;
        getThermalDiffCoeffs(DATA_PTR(m_spwork));
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

    // enforce the condition \sum Y_k V_k = 0. This is done by
    // replacing the flux equation with the largest gradx
    // component with the flux balance condition.
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

    // use LAPACK to solve the equations
    int info=0;
    size_t nr = m_aa.nRows();
    size_t nc = m_aa.nColumns();

    ct_dgetrf(nr, nc, m_aa.ptrColumn(0), nr, &m_aa.ipiv()[0], info);
    if (info == 0) {
        int ndim = 1;
        ct_dgetrs(ctlapack::NoTranspose, nr, ndim,
                  m_aa.ptrColumn(0), nr, &m_aa.ipiv()[0], fluxes, nr, info);
        if (info != 0)
            throw CanteraError("MultiTransport::getMassFluxes",
                               "Error in DGETRS. Info = "+int2str(info));
    } else
        throw CanteraError("MultiTransport::getMassFluxes",
                           "Error in DGETRF.  Info = "+int2str(info));

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
        eval_L0000(DATA_PTR(m_molefracs));
    }

    // invert L00,00
    int ierr = invert(m_Lmatrix, m_nsp);
    if (ierr != 0) {
        throw CanteraError("MultiTransport::getMultiDiffCoeffs",
                           string(" invert returned ierr = ")+int2str(ierr));
    }
    m_l0000_ok = false;           // matrix is overwritten by inverse
    m_lmatrix_soln_ok = false;

    //doublereal pres = m_thermo->pressure();
    doublereal prefactor = 16.0 * m_temp
                           * m_thermo->meanMolecularWeight()/(25.0 * p);
    doublereal c;

    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            c = prefactor/m_mw[j];
            d[ld*j + i] = c*m_molefracs[i]*
                          (m_Lmatrix(i,j) - m_Lmatrix(i,i));
        }
    }
}

void MultiTransport::update_T()
{
    if (m_temp == m_thermo->temperature()) {
        return;
    }

    GasTransport::update_T();

    // temperature has changed, so polynomial fits will need to be
    // redone, and the L matrix reevaluated.
    m_abc_ok  = false;
    m_lmatrix_soln_ok = false;
    m_l0000_ok = false;
}

void MultiTransport::update_C()
{
    // Update the local mole fraction array
    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

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
    doublereal z;
    int ipoly;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            z = m_logt - m_log_eps_k(i,j);
            ipoly = m_poly[i][j];
            if (m_mode == CK_Mode) {
                m_om22(i,j) = poly6(z, DATA_PTR(m_om22_poly[ipoly]));
                m_astar(i,j) = poly6(z, DATA_PTR(m_astar_poly[ipoly]));
                m_bstar(i,j) = poly6(z, DATA_PTR(m_bstar_poly[ipoly]));
                m_cstar(i,j) = poly6(z, DATA_PTR(m_cstar_poly[ipoly]));
            } else {
                m_om22(i,j) = poly8(z, DATA_PTR(m_om22_poly[ipoly]));
                m_astar(i,j) = poly8(z, DATA_PTR(m_astar_poly[ipoly]));
                m_bstar(i,j) = poly8(z, DATA_PTR(m_bstar_poly[ipoly]));
                m_cstar(i,j) = poly8(z, DATA_PTR(m_cstar_poly[ipoly]));
            }
            m_om22(j,i)  = m_om22(i,j);
            m_astar(j,i) = m_astar(i,j);
            m_bstar(j,i) = m_bstar(i,j);
            m_cstar(j,i) = m_cstar(i,j);
        }
    }
    m_abc_ok = true;

    // evaluate the temperature-dependent rotational relaxation rate
    doublereal tr, sqtr;
    for (size_t k = 0; k < m_nsp; k++) {
        tr = m_eps[k]/ m_kbt;
        sqtr = m_sqrt_eps_k[k] / m_sqrt_t;
        m_rotrelax[k] = std::max(1.0,m_zrot[k]) * m_frot_298[k]/Frot(tr, sqtr);
    }

    doublereal d;
    doublereal c = 1.2*GasConstant*m_temp;
    for (size_t k = 0; k < m_nsp; k++) {
        d = c * m_visc[k] * m_astar(k,k)/m_mw[k];
        m_bdiff(k,k) = d;
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
    const vector_fp& cp = ((IdealGasPhase*)m_thermo)->cp_R_ref();
    for (size_t k = 0; k < m_nsp; k++) {
        m_cinternal[k] = cp[k] - 2.5;
    }

    // m_thermo->update_T(m_update_thermal_T);
    m_thermal_tlast = m_thermo->temperature();
}

}
