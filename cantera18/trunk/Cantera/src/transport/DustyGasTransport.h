///
///
///  @file DustyGasTransport.h
///  Interface for class DustyGasTransport
///
///


// Copyright 2003  California Institute of Technology


#ifndef CT_DUSTYGASTRAN_H
#define CT_DUSTYGASTRAN_H

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"


namespace Cantera {

    ///
    /// Class DustyGasTransport implements the Dusty Gas model for
    /// transport in porous media. As implemented here, only species
    /// transport is handled. The viscosity, thermal conductivity, and
    /// thermal diffusion coefficients are not implemented.
    ///
    class DustyGasTransport : public Transport {

    public:

        /// default constructor
        DustyGasTransport(thermo_t* thermo=0);
        
        /// Destructor. Does nothing, since class allocates no memory
        /// on the heap.
        virtual ~DustyGasTransport() {}
        
        
        //---------------------------------------------------------
        // overloaded base class methods

        virtual int model() { return cDustyGasTransport; }

        virtual void setParameters(const int type, const int k, const doublereal* const p);
        
        virtual void getMultiDiffCoeffs(const int ld, doublereal* const d);
        
        virtual void getMolarFluxes(const doublereal* state1,
            const doublereal* state2, doublereal delta, 
            doublereal* fluxes);
        
        //-----------------------------------------------------------
        // new methods added in this class
        
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
        
        /// Set the permeability. If not set, the value for
        /// close-packed spheres will be used by default. 
        void setPermeability(doublereal B) {
            m_perm = B;
        }
        
        /// Return a reference to the transport manager used to compute the gas
        /// binary diffusion coefficients and the visdcosity.
        Transport& gasTransport() { return *m_gastran; }
        
        
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
        doublereal                   m_temp;

        /// multicomponent diffusion coefficients
        DenseMatrix                  m_multidiff;

        // work space
        vector_fp  m_spwork;
        vector_fp  m_spwork2;

        // concentration gradients
        //vector_fp  m_gradConc; 
        //vector_fp  m_conc;

        doublereal m_gradP;   /// pressure gradient

        bool m_knudsen_ok;
        bool m_bulk_ok;
        bool m_conc_set;
        bool m_gradConc_set;
        bool m_gradP_set;

        doublereal m_porosity;      /// porosity
        doublereal m_tortuosity;    /// tortuosity
        doublereal m_pore_radius;   /// pore radius (m)
        doublereal m_diam;          /// particle diameter (m)
        doublereal m_perm;          /// permeability

        Transport* m_gastran;       /// pointer to gas transport manager

    };
}
#endif






