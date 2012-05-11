/**
 *  @file MultiTransport.h
 *  Interface for class MultiTransport
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_MULTITRAN_H
#define CT_MULTITRAN_H

// Define this for better agreement with Chemkin TRANLIB results, even
// if the results are less correct.
//#undef CHEMKIN_COMPATIBILITY_MODE

// Cantera includes
#include "TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{

//====================================================================================================================
//! Transport solve options
//! @deprecated GMRES option is unimplemented.
enum TRANSOLVE_TYPE {
    //!  Solve the dense matrix via a gmres iteration
    TRANSOLVE_GMRES = 1,
    //!  Solve the dense matrix via an LU gauss elimination
    TRANSOLVE_LU
};
//====================================================================================================================
class GasTransportParams;
//====================================================================================================================
//!   Class L_Matrix is used to represent the "L" matrix.
/*!
 *  This class is used instead of DenseMatrix so that a version of mult can be
 * used that knows about the structure of the L matrix,
 * specifically that the upper-right and lower-left blocks are
 * zero.
 * @ingroup transportProps
 */
class L_Matrix : public DenseMatrix
{
public:

    //! default constructor
    L_Matrix() {}

    //! destructor
    virtual ~L_Matrix() {}

    //! Conduct a multiply with the Dense matrix
    /*!
     * This method is used by GMRES to multiply the L matrix by a
     * vector b.  The L matrix has a 3x3 block structure, where each
     * block is a K x K matrix.  The elements of the upper-right and
     * lower-left blocks are all zero.  This method is defined so
     * that the multiplication only involves the seven non-zero
     * blocks.
     *
     *   @param b
     *   @param prod
     *   @deprecated GMRES method is not implemented
     */
    DEPRECATED(virtual void mult(const doublereal* b, doublereal* prod) const);
};


//====================================================================================================================
//! Class MultiTransport implements multicomponent transport
//! properties for ideal gas mixtures.
/*!
 *
 *  The implementation generally
 * follows the procedure outlined in Kee, Coltrin, and Glarborg,
 * "Theoretical and Practical Aspects of Chemically Reacting Flow
 * Modeling," Wiley Interscience.
 *
 * @ingroup transportProps
 */
class MultiTransport : public Transport
{

protected:

    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    MultiTransport(thermo_t* thermo=0);

public:

    //! Destructor
    virtual ~MultiTransport();

    // overloaded base class methods
    virtual int model() const {
        if (m_mode == CK_Mode) {
            return CK_Multicomponent;
        } else {
            return cMulticomponent;
        }
    }

    virtual doublereal viscosity();

    virtual void getSpeciesViscosities(doublereal* const visc) {
        updateViscosity_T();
        std::copy(m_visc.begin(), m_visc.end(), visc);
    }


    //! Return the thermal diffusion coefficients (kg/m/s)
    /*!
     *  Eqn. (12.126) displays how they are calculated. The reference work is from
     *  Dixon-Lewis.
     *
     *  Eqns. (12.168) shows how they are used in an expression for the species flux.
     *
     * @param dt  Vector of thermal diffusion coefficients. Units = kg/m/s
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    virtual doublereal thermalConductivity();

    virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);
    virtual void getMultiDiffCoeffs(const size_t ld, doublereal* const d);

    //! Although this class implements a multicomponent diffusion
    //! model, it is convenient to be able to compute
    //! mixture-averaged diffusion coefficients too.
    /*!
     * @param d Mixture averaged diffusion coefficients
     *          Length = m_msp, units = m2/sec
     */
    virtual void getMixDiffCoeffs(doublereal* const d);

    //! Get the species diffusive mass fluxes wrt to  the mass averaged velocity,
    //! given the gradients in mole fraction and temperature
    /*!
     *  Units for the returned fluxes are kg m-2 s-1.
     *
     *  @param ndim     Number of dimensions in the flux expressions
     *  @param grad_T   Gradient of the temperature
     *                   (length = ndim)
     * @param ldx       Leading dimension of the grad_X array
     *                   (usually equal to m_nsp but not always)
     * @param grad_X    Gradients of the mole fraction
     *                  Flat vector with the m_nsp in the inner loop.
     *                   length = ldx * ndim
     * @param ldf       Leading dimension of the fluxes array
     *                   (usually equal to m_nsp but not always)
     * @param fluxes    Output of the diffusive mass fluxes
     *                  Flat vector with the m_nsp in the inner loop.
     *                   length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                  int ldx,  const doublereal* const grad_X,
                                  int ldf, doublereal* const fluxes);

    //! Get the molar diffusional fluxes [kmol/m^2/s] of the species, given the thermodynamic
    //! state at two nearby points.
    /*!
     * The molar diffusional fluxes are calculated with reference to the mass averaged
     * velocity. This is a one-dimensional vector
     *
     * @param state1 Array of temperature, density, and mass
     *               fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     *               fractions for state 2.
     * @param delta  Distance from state 1 to state 2 (m).
     * @param fluxes Output molar fluxes of the species.
     *               (length = m_nsp)
     */
    virtual void getMolarFluxes(const doublereal* const state1,
                                const doublereal* const state2,
                                const doublereal delta,
                                doublereal* const fluxes);

    //! Get the mass diffusional fluxes [kg/m^2/s] of the species, given the thermodynamic
    //! state at two nearby points.
    /*!
     * The specific diffusional fluxes are calculated with reference to the mass averaged
     * velocity. This is a one-dimensional vector
     *
     * @param state1 Array of temperature, density, and mass
     *               fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     *               fractions for state 2.
     * @param delta  Distance from state 1 to state 2 (m).
     * @param fluxes Output mass fluxes of the species.
     *               (length = m_nsp)
     */
    virtual void getMassFluxes(const doublereal* state1,
                               const doublereal* state2, doublereal delta,
                               doublereal* fluxes);

    //! Set the solution method for inverting the L matrix
    /*!
     *      @param method enum TRANSOLVE_TYPE Either use direct or TRANSOLVE_GMRES
     *      @deprecated GMRES option is unimplemented.
     */
    DEPRECATED(virtual void setSolutionMethod(TRANSOLVE_TYPE method));

    //! Set the options for the GMRES solution
    /*!
     *      @param m    set the mgmres param
     *      @param eps  Set the eps parameter
     *      @deprecated GMRES option is unimplemented.
     */
    DEPRECATED(virtual void setOptions_GMRES(int m, doublereal eps));

    /**
     * @internal
     */

    //! Initialize the transport operator with parameters from GasTransportParams object
    /*!
     *  @param tr  input GasTransportParams object
     */
    virtual bool initGas(GasTransportParams& tr);

    friend class TransportFactory;

    //! Return a structure containing all of the pertinent parameters
    //! about a species that was used to construct the Transport properties in this object
    /*!
     * @param k        Species index
     */
    struct GasTransportData getGasTransportData(int k);

protected:

    //! Update basic temperature-dependent quantities if the temperature has changed.
    void updateTransport_T();

    //! Update basic concentration-dependent quantities if the concentrations have changed.
    void updateTransport_C();

    //! Update the temperature-dependent terms needed to compute the thermal
    //! conductivity and thermal diffusion coefficients.
    void updateThermal_T();

    //! Update the temperature-dependent viscosity terms
    void updateViscosity_T();

    //! Update the temperature-dependent viscosity terms.
    //! Updates the array of pure species viscosities and the weighting
    //! functions in the viscosity mixture rule.
    //! The flag m_visc_ok is set to true.
    void updateSpeciesViscosities_T();

    //! Update the binary diffusion coefficients.
    //! These are evaluated from the polynomial fits at unit pressure (1 Pa).
    void updateDiff_T();

private:

    doublereal m_diff_tlast;
    doublereal m_spvisc_tlast;
    doublereal m_visc_tlast;
    doublereal m_thermal_tlast;

    doublereal m_tmin;
    doublereal m_tmax;
    vector_fp  m_mw;

    // polynomial fits
    std::vector<vector_fp>            m_visccoeffs;
    std::vector<vector_fp>            m_diffcoeffs;
    vector_fp                    m_polytempvec;

    // property values
    DenseMatrix                  m_bdiff;
    vector_fp                    m_visc;
    vector_fp                    m_sqvisc;

    vector_fp                    m_molefracs;


    std::vector<std::vector<int> > m_poly;
    std::vector<vector_fp>   m_astar_poly;
    std::vector<vector_fp>   m_bstar_poly;
    std::vector<vector_fp>   m_cstar_poly;
    std::vector<vector_fp>   m_om22_poly;

    //! Dense matrix for astar
    DenseMatrix          m_astar;

    //! Dense matrix for bstar
    DenseMatrix          m_bstar;

    //! Dense matrix for cstar
    DenseMatrix          m_cstar;

    //! Dense matrix for omega22
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

    bool m_gmres; //!< @deprecated
    int  m_mgmres; //!< @deprecated
    doublereal m_eps_gmres; //!< @deprecated

    // work space
    vector_fp  m_spwork, m_spwork1, m_spwork2, m_spwork3;

    void correctBinDiffCoeffs();

    //! Boolean indicating viscosity is up to date
    bool m_visc_ok;
    bool m_spvisc_ok;
    bool m_diff_ok;
    bool m_abc_ok;
    bool m_l0000_ok;
    bool m_lmatrix_soln_ok;
    int m_mode;

    //! Evaluate the L0000 matrices
    /*!
     *  Evaluate the upper-left block of the L matrix.
     *  @param x vector of species mole fractions
     */
    void eval_L0000(const doublereal* const x);

    //! Evaluate the L0010 matrices
    /*!
     *  @param x vector of species mole fractions
     */
    void eval_L0010(const doublereal* const x);

    //! Evaluate the L1000 matrices
    /*!
     *
     */
    void eval_L1000();

    void eval_L0100();
    void eval_L0001();
    void eval_L1010(const doublereal* x);
    void eval_L1001(const doublereal* x);
    void eval_L0110();
    void eval_L0101(const doublereal* x);
    bool hasInternalModes(size_t j);

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
