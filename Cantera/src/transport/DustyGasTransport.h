/**
 *
 *  @file DustyGasTransport.h
 *  Interface for class DustyGasTransport
 *
 */

// Copyright 2003  California Institute of Technology


#ifndef CT_DUSTYGASTRAN_H
#define CT_DUSTYGASTRAN_H


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

using namespace std;

// Cantera includes
#include "TransportBase.h"
#include "../DenseMatrix.h"


namespace Cantera {


    /**
     * Class DustyGasTransport implements the Dusty Gas model for
     * transport in porous media. As implemented here, only species
     * transport is handled. The viscosity, thermal conductivity, and
     * thermal diffusion coefficients are not implemented.
    */
    class DustyGasTransport : public Transport {

    public:


        /// default constructor
      DustyGasTransport(thermo_t* thermo=0);
        virtual ~DustyGasTransport() {}

        // overloaded base class methods

        virtual int model() { return cDustyGasTransport; }

        virtual void setParameters(int type, int k, doublereal* p);


        //virtual void getBinaryDiffCoeffs(int ld, doublereal* d);

        /**
         * Get the multicomponent effective diffusion coefficients.
         */
        virtual void getMultiDiffCoeffs(int ld, doublereal* d);


        // new methods added in this class

        /**
         * Get the molar gas species fluxes. These fluxes include both the ordinary mass diffusion component
         * and the Darcy (pressure-driven) commponent. 
         */
        void getMolarFluxes(const double* grad_conc,
            double grad_P, double* fluxes);

        /// Set the porosity (dimensionless)
        void setPorosity(doublereal porosity) {
            m_porosity = porosity;
            m_knudsen_ok = false;
            m_bulk_ok = false;
        }

        /// Set the tortuosity (dimensionless)
        void setTortuosity(doublereal tort) {
            m_tortuosity = tort;
            m_knudsen_ok = false;
            m_bulk_ok = false;
        }

        /// Set the mean pore radius (m)
        void setMeanPoreRadius(doublereal rbar) {
            m_pore_radius = rbar;
            m_knudsen_ok = false;
        }

        /// Set the mean particle diameter
        void setMeanParticleDiameter(doublereal dbar) {
            m_diam = dbar;
        }

        /** 
         * Set the permeability. If not set, the value for
         * close-packed spheres will be used by default. 
         */ 
        void setPermeability(doublereal B) {
            m_perm = B;
        }

        /**
         * @internal
         */

        friend class TransportFactory;

    protected:

        // called by TransportFactory
        void initialize(ThermoPhase* phase, Transport* gastr);

    private:

        void updateTransport_T();
        void updateTransport_C();

        void updateBinaryDiffCoeffs();
        void updateMultiDiffCoeffs();
        void updateKnudsenDiffCoeffs();
        void eval_H_matrix();


        // gas attributes
        int m_nsp;
        doublereal m_tmin, m_tmax;
        vector_fp  m_mw;

        // property values

        /// binary diffusion coefficients
        DenseMatrix                  m_d;

        /// mole fractions
        vector_fp                    m_x;

        /// Knudsen diffusion coefficients
        vector_fp                    m_dk;

        /// temperature
        doublereal m_temp;

        /// multicomponent diffusion coefficients
        DenseMatrix  m_multidiff;

        // work space
        vector_fp  m_spwork;

        bool m_knudsen_ok;
        bool m_bulk_ok;

        doublereal m_porosity;
        doublereal m_tortuosity;
        doublereal m_pore_radius;
        doublereal m_diam;
        doublereal m_perm;

        Transport* m_gastran;

    };
}
#endif






