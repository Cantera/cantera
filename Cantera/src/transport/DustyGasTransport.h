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


    class TransportParams;

    class DustyGasTransport : public Transport {

    public:


        virtual ~DustyGasTransport() {}

        // overloaded base class methods
        virtual int model() { return cDustyGasTransport; }

        virtual void setParameters(int type, int k, doublereal* p) {
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

        virtual void getBinaryDiffCoeffs(int ld, doublereal* d);
        virtual void getMultiDiffCoeffs(int ld, doublereal* d);


        // new methods

        void getMolarFluxes(const double* grad_conc,
            double grad_P, double* fluxes);

        void setPorosity(doublereal porosity) {
            m_porosity = porosity;
            m_knudsen_ok = false;
            m_bulk_ok = false;
        }

        void setTortuosity(doublereal tort) {
            m_tortuosity = tort;
            m_knudsen_ok = false;
            m_bulk_ok = false;
        }

        void setMeanPoreRadius(doublereal rbar) {
            m_pore_radius = rbar;
            m_knudsen_ok = false;
        }

        void setMeanParticleDiameter(doublereal dbar) {
            m_diam = dbar;
        }

        void setPermeability(doublereal B) {
            m_perm = B;
        }

        /**
         * @internal
         */
        virtual bool init(TransportParams& tr);

        void updateTransport_T();
        void updateTransport_C();

        friend class TransportFactory;

    protected:

        void updateBinaryDiffCoeffs();
        void updateMultiDiffCoeffs();
        void updateKnudsenDiffCoeffs();
        void eval_H_matrix();
        /// default constructor
        DustyGasTransport(thermo_t* thermo=0);
        void initialize(ThermoPhase* phase, Transport* gastr);
    private:


        // mixture attributes
        int m_nsp;
        doublereal m_tmin, m_tmax;
        vector_fp  m_mw;

        // property values
        DenseMatrix                  m_d;
        vector_fp                    m_visc;
        vector_fp                    m_x;
        vector_fp                    m_dk;
        doublereal m_temp;

        // H matrix quantities
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

        doublereal pressure_ig() {
            return m_thermo->molarDensity() * GasConstant * m_thermo->temperature();
        }
    };
}
#endif






