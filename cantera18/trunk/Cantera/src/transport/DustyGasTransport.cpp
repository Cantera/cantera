/**
 *
 *  @file DustyGasTransport.cpp
 *  Implementation file for class DustyGasTransport
 *
 *  @ingroup transportProps
 *
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2009/03/27 18:24:39 $
 *  $Revision: 1.12 $
 *
 *  Copyright 2003 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "DustyGasTransport.h"

using namespace std;

/** 
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-20


namespace Cantera {


    //////////////////// class DustyGasTransport methods //////////////


    DustyGasTransport::DustyGasTransport(thermo_t* thermo) 
        : Transport(thermo),
	  m_temp(-1.0), 
          m_porosity(0.0),
          m_tortuosity(1.0),
          m_pore_radius(0.0),
          m_diam(0.0),
          m_perm(-1.0),
          m_gastran(0)
    {}

        void DustyGasTransport::setParameters(const int type, const int k, const doublereal* const p) {
            switch(type) {
            case 0:
                setPorosity(p[0]); break;
            case 1:
                setTortuosity(p[0]); break;
            case 2:
                setMeanPoreRadius(p[0]); break;
            case 3:
                setMeanParticleDiameter(p[0]); break;
            case 4:
                setPermeability(p[0]); break;
            default:
                throw CanteraError("DustyGasTransport::init",
                    "unknown parameter");
            }
        }

    void DustyGasTransport::initialize(ThermoPhase* phase, Transport* gastr) {

        // constant mixture attributes
        m_thermo = phase;
        m_nsp   = m_thermo->nSpecies();
        m_tmin  = m_thermo->minTemp();
        m_tmax  = m_thermo->maxTemp();

        m_gastran = gastr;

        // make a local copy of the molecular weights
        m_mw.resize(m_nsp);
        copy(m_thermo->molecularWeights().begin(), 
            m_thermo->molecularWeights().end(), m_mw.begin());

        m_multidiff.resize(m_nsp, m_nsp);
        m_d.resize(m_nsp, m_nsp);
        m_dk.resize(m_nsp, 0.0);
        m_x.resize(m_nsp);
        //m_gradConc.resize(m_nsp);
        //m_conc.resize(m_nsp);

        // set flags all false
        m_knudsen_ok = false;
        m_bulk_ok = false;

        // some work space
        m_spwork.resize(m_nsp);
        m_spwork2.resize(m_nsp);
    }


    /******************* binary diffusion coefficients **************/


    void DustyGasTransport::updateBinaryDiffCoeffs() {
        if (m_bulk_ok) return;
        int n,m;

        // get the gaseous binary diffusion coefficients
        m_gastran->getBinaryDiffCoeffs(m_nsp, m_d.ptrColumn(0));
        doublereal por2tort = m_porosity / m_tortuosity;
        for (n = 0; n < m_nsp; n++) 
            for (m = 0; m < m_nsp; m++) 
                m_d(n,m) *= por2tort;
        m_bulk_ok = true;
    }

    void DustyGasTransport::updateKnudsenDiffCoeffs() {
        if (m_knudsen_ok) return;
        doublereal K_g = m_pore_radius * m_porosity / m_tortuosity;
        const doublereal TwoThirds = 2.0/3.0;
        for (int k = 0; k < m_nsp; k++) {
            m_dk[k] = TwoThirds * K_g * sqrt((8.0 * GasConstant * m_temp)/
                (Pi * m_mw[k]));
        }
        m_knudsen_ok = true;
    }

        
    void DustyGasTransport::eval_H_matrix() {
        updateBinaryDiffCoeffs();
        updateKnudsenDiffCoeffs();
        int k,l,j;
        doublereal sum;
        for (k = 0; k < m_nsp; k++) {

            // evaluate off-diagonal terms
            for (l = 0; l < m_nsp; l++) m_multidiff(k,l) = -m_x[k]/m_d(k,l);

            // evaluate diagonal term
            sum = 0.0;
            for (j = 0; j < m_nsp; j++) if (j != k) sum += m_x[j]/m_d(k,j);
            m_multidiff(k,k) = 1.0/m_dk[k] + sum;
        }
    }

//     void DustyGasTransport::getMolarFluxes(const double* grad_conc,
//         double grad_P, double* fluxes) {
//         updateMultiDiffCoeffs();
//         copy(grad_conc, grad_conc + m_nsp, m_spwork.begin());
//         multiply(m_multidiff, m_spwork.begin(), fluxes);
//         m_thermo->getConcentrations(m_spwork.begin());
//         divide_each(m_spwork.begin(), m_spwork.end(), m_dk.begin());

//         // if no permeability has been specified, use result for 
//         // close-packed spheres
//         double b = 0.0;
//         if (m_perm < 0.0) {
//             double p = m_porosity;
//             double d = m_diam;
//             double t = m_tortuosity;
//             b = p*p*p*d*d/(72.0*t*(1.0-p)*(1.0-p));
//         }
//         else {
//             b = m_perm;
//         }
//         b *= grad_P / m_gastran->viscosity();
//         scale(m_spwork.begin(), m_spwork.end(), m_spwork.begin(), b);
//         increment(m_multidiff, m_spwork.begin(), fluxes);
//         scale(fluxes, fluxes + m_nsp, fluxes, -1.0);
//     }


    void DustyGasTransport::getMolarFluxes(const doublereal* state1,
        const doublereal* state2, double delta, double* fluxes) {
        int k;
        doublereal conc1, conc2;
        doublereal* cbar = DATA_PTR(m_spwork);
        doublereal* gradc = DATA_PTR(m_spwork2);
        doublereal t1 = state1[0];
        doublereal t2 = state2[0];
        doublereal rho1 = state1[1];
        doublereal rho2 = state2[1];
        const doublereal* y1 = state1 + 2;
        const doublereal* y2 = state2 + 2;
        doublereal c1sum = 0.0, c2sum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            conc1 = rho1*y1[k]/m_mw[k];
            conc2 = rho2*y2[k]/m_mw[k];
            cbar[k] = 0.5*(conc1 + conc2);
            gradc[k] = (conc2 - conc1)/delta;
            c1sum += conc1;
            c2sum += conc2;
        }
        doublereal p1 = c1sum * GasConstant * state1[0];
        doublereal p2 = c2sum * GasConstant * state2[0];
        doublereal pbar = 0.5*(p1 + p2);
        doublereal gradp = (p2 - p1)/delta;
        doublereal tbar = 0.5*(t1 + t2);

        m_thermo->setState_TPX(tbar, pbar, cbar);

        updateMultiDiffCoeffs();

        multiply(m_multidiff, gradc, fluxes);
        divide_each(cbar, cbar + m_nsp, m_dk.begin());

        // if no permeability has been specified, use result for 
        // close-packed spheres
        double b = 0.0;
        if (m_perm < 0.0) {
            double p = m_porosity;
            double d = m_diam;
            double t = m_tortuosity;
            b = p*p*p*d*d/(72.0*t*(1.0-p)*(1.0-p));
        }
        else {
            b = m_perm;
        }
        b *= gradp / m_gastran->viscosity();
        scale(cbar, cbar + m_nsp, cbar, b);
        increment(m_multidiff, cbar, fluxes);
        scale(fluxes, fluxes + m_nsp, fluxes, -1.0);
    }


    void DustyGasTransport::updateMultiDiffCoeffs() {
        // see if temperature has changed
        updateTransport_T();

        // update the mole fractions
        updateTransport_C();

        eval_H_matrix();

        // invert H
        int ierr = invert(m_multidiff);

        if (ierr != 0) {
            throw CanteraError("DustyGasTransport::updateMultiDiffCoeffs",
                "invert returned ierr = "+int2str(ierr));
        }
    }

    void DustyGasTransport::getMultiDiffCoeffs(const int ld, doublereal* const d) {
        int i,j;
        updateMultiDiffCoeffs();
        for (i = 0; i < m_nsp; i++) {
            for (j = 0; j < m_nsp; j++) {            
                d[ld*j + i] = m_multidiff(i,j);
            }
        }
    }

              
    /**
     *  Update temperature-dependent quantities.
     */ 
    void DustyGasTransport::updateTransport_T() 
    {
        if (m_temp == m_thermo->temperature()) return;
        m_temp = m_thermo->temperature();
        m_knudsen_ok = false;
        m_bulk_ok = false;
    }                 

    void DustyGasTransport::updateTransport_C()  
    {
        m_thermo->getMoleFractions(DATA_PTR(m_x));

        // add an offset to avoid a pure species condition
        // (check - this may be unnecessary)
        int k;
        for (k = 0; k < m_nsp; k++) {
            m_x[k] = fmaxx(MIN_X, m_x[k]);
        }
    }
}
