/**
 *  @file DustyGasTransport.cpp
 *  Implementation file for class DustyGasTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/DustyGasTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
DustyGasTransport::DustyGasTransport(thermo_t* thermo) :
    Transport(thermo),
    m_temp(-1.0),
    m_gradP(0.0),
    m_knudsen_ok(false),
    m_bulk_ok(false),
    m_porosity(0.0),
    m_tortuosity(1.0),
    m_pore_radius(0.0),
    m_diam(0.0),
    m_perm(-1.0)
{
}

void DustyGasTransport::setThermo(thermo_t& thermo)
{
    Transport::setThermo(thermo);
    m_gastran->setThermo(thermo);
}

void DustyGasTransport::initialize(ThermoPhase* phase, Transport* gastr)
{
    // constant mixture attributes
    m_thermo = phase;
    m_nsp = m_thermo->nSpecies();
    if (m_gastran.get() != gastr) {
        m_gastran.reset(gastr);
    }

    // make a local copy of the molecular weights
    m_mw = m_thermo->molecularWeights();

    m_multidiff.resize(m_nsp, m_nsp);
    m_d.resize(m_nsp, m_nsp);
    m_dk.resize(m_nsp, 0.0);

    m_x.resize(m_nsp, 0.0);
    m_thermo->getMoleFractions(m_x.data());

    // set flags all false
    m_knudsen_ok = false;
    m_bulk_ok = false;

    m_spwork.resize(m_nsp);
    m_spwork2.resize(m_nsp);
}

void DustyGasTransport::updateBinaryDiffCoeffs()
{
    if (m_bulk_ok) {
        return;
    }

    // get the gaseous binary diffusion coefficients
    m_gastran->getBinaryDiffCoeffs(m_nsp, m_d.ptrColumn(0));
    doublereal por2tort = m_porosity / m_tortuosity;
    for (size_t n = 0; n < m_nsp; n++) {
        for (size_t m = 0; m < m_nsp; m++) {
            m_d(n,m) *= por2tort;
        }
    }
    m_bulk_ok = true;
}

void DustyGasTransport::updateKnudsenDiffCoeffs()
{
    if (m_knudsen_ok) {
        return;
    }
    doublereal K_g = m_pore_radius * m_porosity / m_tortuosity;
    for (size_t k = 0; k < m_nsp; k++) {
        m_dk[k] = 2.0/3.0 * K_g * sqrt((8.0 * GasConstant * m_temp)/
                                         (Pi * m_mw[k]));
    }
    m_knudsen_ok = true;
}

void DustyGasTransport::eval_H_matrix()
{
    updateBinaryDiffCoeffs();
    updateKnudsenDiffCoeffs();
    for (size_t k = 0; k < m_nsp; k++) {
        // evaluate off-diagonal terms
        for (size_t j = 0; j < m_nsp; j++) {
            m_multidiff(k,j) = -m_x[k]/m_d(k,j);
        }

        // evaluate diagonal term
        double sum = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            if (j != k) {
                sum += m_x[j]/m_d(k,j);
            }
        }
        m_multidiff(k,k) = 1.0/m_dk[k] + sum;
    }
}

void DustyGasTransport::getMolarFluxes(const doublereal* const state1,
                                       const doublereal* const state2,
                                       const doublereal delta,
                                       doublereal* const fluxes)
{
    // cbar will be the average concentration between the two points
    doublereal* const cbar = m_spwork.data();
    doublereal* const gradc = m_spwork2.data();
    const doublereal t1 = state1[0];
    const doublereal t2 = state2[0];
    const doublereal rho1 = state1[1];
    const doublereal rho2 = state2[1];
    const doublereal* const y1 = state1 + 2;
    const doublereal* const y2 = state2 + 2;
    doublereal c1sum = 0.0, c2sum = 0.0;

    for (size_t k = 0; k < m_nsp; k++) {
        double conc1 = rho1 * y1[k] / m_mw[k];
        double conc2 = rho2 * y2[k] / m_mw[k];
        cbar[k] = 0.5*(conc1 + conc2);
        gradc[k] = (conc2 - conc1) / delta;
        c1sum += conc1;
        c2sum += conc2;
    }

    // Calculate the pressures at p1 p2 and pbar
    doublereal p1 = c1sum * GasConstant * t1;
    doublereal p2 = c2sum * GasConstant * t2;
    doublereal pbar = 0.5*(p1 + p2);
    doublereal gradp = (p2 - p1)/delta;
    doublereal tbar = 0.5*(t1 + t2);
    m_thermo->setState_TPX(tbar, pbar, cbar);
    updateMultiDiffCoeffs();

    // Multiply m_multidiff and gradc together and store the result in fluxes[]
    multiply(m_multidiff, gradc, fluxes);
    for (size_t k = 0; k < m_nsp; k++) {
        cbar[k] /= m_dk[k];
    }

    // if no permeability has been specified, use result for
    // close-packed spheres
    double b = 0.0;
    if (m_perm < 0.0) {
        double p = m_porosity;
        double d = m_diam;
        double t = m_tortuosity;
        b = p*p*p*d*d/(72.0*t*(1.0-p)*(1.0-p));
    } else {
        b = m_perm;
    }
    b *= gradp / m_gastran->viscosity();
    scale(cbar, cbar + m_nsp, cbar, b);

    // Multiply m_multidiff with cbar and add it to fluxes
    increment(m_multidiff, cbar, fluxes);
    scale(fluxes, fluxes + m_nsp, fluxes, -1.0);
}

void DustyGasTransport::updateMultiDiffCoeffs()
{
    // see if temperature has changed
    updateTransport_T();

    // update the mole fractions
    updateTransport_C();
    eval_H_matrix();

    // invert H
    int ierr = invert(m_multidiff);
    if (ierr != 0) {
        throw CanteraError("DustyGasTransport::updateMultiDiffCoeffs",
                           "invert returned ierr = {}", ierr);
    }
}

void DustyGasTransport::getMultiDiffCoeffs(const size_t ld, doublereal* const d)
{
    updateMultiDiffCoeffs();
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = m_multidiff(i,j);
        }
    }
}

void DustyGasTransport::updateTransport_T()
{
    if (m_temp == m_thermo->temperature()) {
        return;
    }
    m_temp = m_thermo->temperature();
    m_knudsen_ok = false;
    m_bulk_ok = false;
}

void DustyGasTransport::updateTransport_C()
{
    m_thermo->getMoleFractions(m_x.data());

    // add an offset to avoid a pure species condition
    // (check - this may be unnecessary)
    for (size_t k = 0; k < m_nsp; k++) {
        m_x[k] = std::max(Tiny, m_x[k]);
    }
    // diffusion coeffs depend on Pressure
    m_bulk_ok = false;
}

void DustyGasTransport::setPorosity(doublereal porosity)
{
    m_porosity = porosity;
    m_knudsen_ok = false;
    m_bulk_ok = false;
}

void DustyGasTransport::setTortuosity(doublereal tort)
{
    m_tortuosity = tort;
    m_knudsen_ok = false;
    m_bulk_ok = false;
}

void DustyGasTransport::setMeanPoreRadius(doublereal rbar)
{
    m_pore_radius = rbar;
    m_knudsen_ok = false;
}

void DustyGasTransport::setMeanParticleDiameter(doublereal dbar)
{
    m_diam = dbar;
}

void DustyGasTransport::setPermeability(doublereal B)
{
    m_perm = B;
}

Transport& DustyGasTransport::gasTransport()
{
    return *m_gastran;
}

}
