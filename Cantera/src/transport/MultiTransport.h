/**
 *
 *  @file MultiTransport.h
 *  Interface for class MultiTransport
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_MULTITRAN_H
#define CT_MULTITRAN_H


// Define this for better agreement with Chemkin TRANLIB results, even
// if the results are less correct.
//#undef CHEMKIN_COMPATIBILITY_MODE


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"


namespace Cantera {


    class TransportParams;

    /////////////////////////////////////////////////////////////

    /**
     * Class L_Matrix is used to represent the "L" matrix.  This class
     * is used instead of DenseMatrix so that a version of mult can be
     * used that knows about the structure of the L matrix,
     * specifically that the upper-right and lower-left blocks are
     * zero.
     * @ingroup transportProps
     */
    class L_Matrix : public DenseMatrix {
    public:
        L_Matrix() {}
        virtual ~L_Matrix(){}

        /**
         * This method is used by GMRES to multiply the L matrix by a
         * vector b.  The L matrix has a 3x3 block structure, where each
         * block is a K x K matrix.  The elements of the upper-right and
         * lower-left blocks are all zero.  This method is defined so
         * that the multiplication only involves the seven non-zero
         * blocks.
         */
        virtual void mult(const doublereal* b, doublereal* prod) const;
    };


    const int GMRES = 1, LU = 2;

    /**
     * Class MultiTransport implements multicomponent transport
     * properties for ideal gas mixtures. The implementation generally
     * follows the procedure outlined in Kee, Coltrin, and Glarborg,
     * "Theoretical and Practical Aspects of Chemically Reacting Flow
     * Modeling," Wiley Interscience.
     * @ingroup transportProps
     */
    class MultiTransport : public Transport {

    public:


      virtual ~MultiTransport();

        // overloaded base class methods
        virtual int model() {
            if (m_mode == CK_Mode)
                return CK_Multicomponent;
            else
                return cMulticomponent;
        }

        virtual doublereal viscosity();

        virtual void getSpeciesViscosities(doublereal* const visc)
            { updateViscosity_T(); std::copy(m_visc.begin(), m_visc.end(), visc); }

        virtual void getThermalDiffCoeffs(doublereal* const dt);
        virtual doublereal thermalConductivity();

        virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d);
        virtual void getMultiDiffCoeffs(const int ld, doublereal* const d);

        //! Although this class implements a multicomponent diffusion
        //! model, it is convenient to be able to compute
        //! mixture-averaged diffusion coefficients too.
        /*!
         * @param d Mixture averaged diffusion coefficients
	 *          Length = m_msp, units = m2/sec
         */
        virtual void getMixDiffCoeffs(doublereal* const d);

        //! Get the species diffusive mass fluxes wrt to 
        //! the mass averaged velocity, 
        //! given the gradients in mole fraction and temperature
        /*!
         *  Units for the returned fluxes are kg m-2 s-1.
         * 
         *  @param ndim Number of dimensions in the flux expressions
         *  @param grad_T Gradient of the temperature
	 *                 (length = ndim)
	 * @param ldx  Leading dimension of the grad_X array 
	 *              (usually equal to m_nsp but not always)
	 * @param grad_X Gradients of the mole fraction
	 *             Flat vector with the m_nsp in the inner loop.
	 *             length = ldx * ndim
	 * @param ldf  Leading dimension of the fluxes array 
	 *              (usually equal to m_nsp but not always)
	 * @param fluxes  Output of the diffusive mass fluxes
	 *             Flat vector with the m_nsp in the inner loop.
	 *             length = ldx * ndim
	 */
        virtual void getSpeciesFluxes(int ndim,
				      const doublereal* grad_T, 
				      int ldx, 
				      const doublereal* grad_X,
				      int ldf,
				      doublereal* fluxes);

        virtual void getMolarFluxes(const doublereal* state1,
            const doublereal* state2, doublereal delta,
            doublereal* fluxes);

        virtual void getMassFluxes(const doublereal* state1,
            const doublereal* state2, doublereal delta,
            doublereal* fluxes);

        virtual void setSolutionMethod(int method) {
            if (method == GMRES) m_gmres = true;
            else m_gmres = false;
        }

        virtual void setOptions_GMRES(int m, doublereal eps) {
            if (m > 0) m_mgmres = m;
            if (eps > 0.0) m_eps_gmres = eps;
        }

        void save(std::string outfile);

        /**
         * @internal
         */
        virtual bool init(TransportParams& tr);


        /**
         * @name Property Updating This methods are used to update
         * temperature- or concentration-dependent quantities. The
         * methods of the first group (with names that do not begin
         * with an underscore) invoke the 'update' method of the
         * relevant property updater. These methods are the ones that
         * are called by other methods of the class to update
         * properties.  The methods that actually perform the updates
         * are the ones with names beginning with an underscore. These
         * are only called by the property updaters.
         */
        void updateTransport_T();
        void updateTransport_C();

        void updateThermal_T();
        void updateViscosity_T();
        void updateSpeciesViscosities_T();
        void updateDiff_T();


        void _update_transport_T();
        void _update_transport_C();
        void _update_species_visc_T();
        void _update_visc_T();
        void _update_diff_T();
        void _update_thermal_T();

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
        MultiTransport(thermo_t* thermo=0);

    private:

//         int m_update_transport_T;
//         int m_update_transport_C;
//         int m_update_spvisc_T;
//         int m_update_visc_T;
//         int m_update_diff_T;
//         int m_update_thermal_T;

        doublereal m_diff_tlast, m_spvisc_tlast, m_visc_tlast,
            m_thermal_tlast;

        // mixture attributes
        int m_nsp;
        doublereal m_tmin, m_tmax;
        vector_fp  m_mw;

        // polynomial fits
        std::vector<vector_fp>            m_visccoeffs;
        std::vector<vector_fp>            m_diffcoeffs;
        vector_fp                    m_polytempvec;

        // property values
        DenseMatrix                  m_bdiff;
        vector_fp                    m_visc;
        vector_fp                    m_sqvisc;

        array_fp                    m_molefracs;


        std::vector<std::vector<int> > m_poly;
        std::vector<vector_fp >   m_astar_poly;
        std::vector<vector_fp >   m_bstar_poly;
        std::vector<vector_fp >   m_cstar_poly;
        std::vector<vector_fp >   m_om22_poly;
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

        // L matrix quantities
        L_Matrix  m_Lmatrix;
        DenseMatrix m_aa;
        //DenseMatrix m_Lmatrix;
        vector_fp m_a;
        vector_fp m_b;

        bool m_gmres;
        int  m_mgmres;
        doublereal m_eps_gmres;

        // work space
        vector_fp  m_spwork, m_spwork1, m_spwork2, m_spwork3;

        void correctBinDiffCoeffs();
        bool m_visc_ok;
        bool m_spvisc_ok;
        bool m_diff_ok;
        bool m_abc_ok;
        bool m_l0000_ok;
        bool m_lmatrix_soln_ok;
        int m_mode;

         void eval_L0000(const doublereal* x);
         void eval_L0010(const doublereal* x);
         void eval_L1000();
         void eval_L0100();
         void eval_L0001();
         void eval_L1010(const doublereal* x);
         void eval_L1001(const doublereal* x);
         void eval_L0110();
         void eval_L0101(const doublereal* x);
         bool hasInternalModes(int j);

        doublereal pressure_ig() {
            return m_thermo->molarDensity() * GasConstant * m_thermo->temperature();
        }

        void solveLMatrixEquation();
        DenseMatrix m_epsilon;
        DenseMatrix m_diam;
        DenseMatrix incl;
        bool m_debug;
    };
}
#endif
