/**
 *  @file MultiTransport.h
 *  Interface for class MultiTransport
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_MULTITRAN_H
#define CT_MULTITRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{

class GasTransportParams;

//! Class MultiTransport implements multicomponent transport properties for
//! ideal gas mixtures.
/*!
 * The implementation generally follows the procedure outlined in Kee,
 * Coltrin, and Glarborg, "Theoretical and Practical Aspects of Chemically
 * Reacting Flow Modeling," Wiley Interscience.
 *
 * @ingroup tranprops
 */
class MultiTransport : public GasTransport
{
protected:

    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    MultiTransport(thermo_t* thermo=0);

public:
    virtual int model() const {
        if (m_mode == CK_Mode) {
            return CK_Multicomponent;
        } else {
            return cMulticomponent;
        }
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

    virtual void getMultiDiffCoeffs(const size_t ld, doublereal* const d);

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
                                  size_t ldx,  const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

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

    //! Get the mass diffusional fluxes [kg/m^2/s] of the species, given the
    //! thermodynamic state at two nearby points.
    /*!
     * The specific diffusional fluxes are calculated with reference to the
     * mass averaged velocity. This is a one-dimensional vector
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

    //! Initialize the transport operator with parameters from GasTransportParams object
    /*!
     *  @param tr  input GasTransportParams object
     */
    virtual bool initGas(GasTransportParams& tr);

    friend class TransportFactory;

protected:
    //! Update basic temperature-dependent quantities if the temperature has changed.
    void update_T();

    //! Update basic concentration-dependent quantities if the concentrations have changed.
    void update_C();

    //! Update the temperature-dependent terms needed to compute the thermal
    //! conductivity and thermal diffusion coefficients.
    void updateThermal_T();

    doublereal m_thermal_tlast;

    // property values
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

public:
    vector_fp   m_crot;
    vector_fp   m_cinternal;
    vector_fp   m_zrot;
    vector_fp   m_eps;
    vector_fp   m_sigma;
    vector_fp   m_alpha;
    DenseMatrix   m_dipole;

protected:
    vector_fp  m_sqrt_eps_k;
    DenseMatrix m_log_eps_k;
    vector_fp  m_frot_298;
    vector_fp  m_rotrelax;

    doublereal m_lambda;

    // L matrix quantities
    DenseMatrix  m_Lmatrix;
    DenseMatrix m_aa;
    //DenseMatrix m_Lmatrix;
    vector_fp m_a;
    vector_fp m_b;

    // work space
    vector_fp m_spwork1, m_spwork2, m_spwork3;

    //! Mole fraction vector from last L-matrix evaluation
    vector_fp m_molefracs_last;

    void correctBinDiffCoeffs();

    //! Boolean indicating viscosity is up to date
    bool m_abc_ok;
    bool m_l0000_ok;
    bool m_lmatrix_soln_ok;

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

    virtual void solveLMatrixEquation();
    DenseMatrix incl;
    bool m_debug;
};
}
#endif
