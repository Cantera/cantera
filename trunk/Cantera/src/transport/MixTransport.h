/**
 *
 *  @file MixTransport.h
 *   Header file defining class MixTransport
 */

/* $Author: dggoodwin $
 * $Revision: 1.4 $
 * $Date: 2005/08/18 14:44:35 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_MIXTRAN_H
#define CT_MIXTRAN_H


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

    /**
     * Class MixTransport implements mixture-averaged transport
     * properties for ideal gas mixtures. The model is based on that
     * described by Kee, Coltrin, and Glarborg, "Theoretical and
     * Practical Aspects of Chemically Reacting Flow Modeling."
     */
    class MixTransport : public Transport {

    public:

        virtual ~MixTransport() {}

        virtual int model() { return cMixtureAveraged; }

        // overloaded base class methods
        virtual doublereal viscosity();


        virtual void getSpeciesViscosities(doublereal* visc)
            { updateViscosity_T(); copy(m_visc.begin(), m_visc.end(), visc); }

        virtual void getThermalDiffCoeffs(doublereal* dt);
        virtual doublereal thermalConductivity();

        virtual void getBinaryDiffCoeffs(int ld, doublereal* d);
        virtual void getMixDiffCoeffs(doublereal* d);
        virtual void getMobilities(doublereal* mobil);
        virtual void update_T();
        virtual void update_C();

        virtual void getSpeciesFluxes(int ndim, 
        doublereal* grad_T, int ldx, const doublereal* grad_X, 
            int ldf, doublereal* fluxes);

        virtual bool init(TransportParams& tr);

        friend class TransportFactory;

	/**
	 * Return a structure containing all of the pertinent parameters
	 * about a species that was used to construct the Transport
	 * properties in this object.
	 *
	 * @param k Species number to obtain the properties from.
	 */
	struct GasTransportData getGasTransportData(int);

    protected:

        /// default constructor
        MixTransport();

    private:


        doublereal pressure_ig() {
            return (m_thermo->molarDensity() * GasConstant *
		    m_thermo->temperature());
        }

        // mixture attributes
        int m_nsp;
        doublereal m_tmin, m_tmax;
        vector_fp  m_mw;

        // polynomial fits
        vector<vector_fp>            m_visccoeffs;
        vector<vector_fp>            m_condcoeffs;
        vector<vector_fp>            m_diffcoeffs;
        vector_fp                    m_polytempvec;

        // property values
        DenseMatrix                  m_bdiff;
        vector_fp                    m_visc;
        vector_fp                    m_sqvisc;
        vector_fp                    m_cond;

        array_fp                    m_molefracs;

        vector<vector<int> > m_poly;
        vector<vector_fp >   m_astar_poly;
        vector<vector_fp >   m_bstar_poly;
        vector<vector_fp >   m_cstar_poly;
        vector<vector_fp >   m_om22_poly;
        DenseMatrix          m_astar;
        DenseMatrix          m_bstar;
        DenseMatrix          m_cstar;
        DenseMatrix          m_om22;

        DenseMatrix m_phi;            // viscosity weighting functions
        DenseMatrix m_wratjk, m_wratkj1;

        vector_fp   m_zrot;
        vector_fp   m_crot;
        vector_fp   m_cinternal;
        vector_fp   m_eps;
	vector_fp   m_alpha;
	vector_fp   m_dipoleDiag;

        doublereal m_temp, m_logt, m_kbt, m_t14, m_t32;
        doublereal m_sqrt_kbt, m_sqrt_t;

        vector_fp  m_sqrt_eps_k;
        DenseMatrix m_log_eps_k;
        vector_fp  m_frot_298;
        vector_fp  m_rotrelax;

        doublereal m_lambda;
        doublereal m_viscmix;

        // work space
        vector_fp  m_spwork;

        void updateThermal_T();
        void updateViscosity_T();
        void updateCond_T();
        void updateSpeciesViscosities();
        void updateDiff_T();
        void correctBinDiffCoeffs();
        bool m_viscmix_ok;
        bool m_viscwt_ok;
        bool m_spvisc_ok;
        bool m_diffmix_ok;
        bool m_bindiff_ok;
        bool m_abc_ok;
        bool m_spcond_ok;
        bool m_condmix_ok;

        int m_mode;

        DenseMatrix m_epsilon;
        DenseMatrix m_diam;
        DenseMatrix incl;
        bool m_debug;
    };
}
#endif






