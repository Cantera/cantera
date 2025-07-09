/**
 *  @file MultiTransport.h
 *  Interface for class MultiTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTITRAN_H
#define CT_MULTITRAN_H

// Cantera includes
#include "GasTransport.h"

namespace Cantera
{
//! Class MultiTransport implements multicomponent transport properties for
//! ideal gas mixtures.
/*!
 * The implementation generally follows the procedure outlined in Kee, et al.
 * @cite kee2003.
 *
 * @ingroup tranprops
 */
class MultiTransport : public GasTransport
{
public:
    //! default constructor
    MultiTransport() = default;

    string transportModel() const override {
        return (m_mode == CK_Mode) ? "multicomponent-CK" : "multicomponent";
    }

    //! Return the thermal diffusion coefficients (kg/m/s)
    /*!
     * Eqn. (12.126) of Kee et al. @cite kee2003 displays how they are calculated. The
     * reference work is from Dixon-Lewis @cite dixon-lewis1968.
     *
     * Eqns. (12.168) of Kee et al. @cite kee2003 shows how they are used in an
     * expression for the species flux.
     *
     * @param dt  Vector of thermal diffusion coefficients. Units = kg/m/s
     */
    void getThermalDiffCoeffs(double* const dt) override;

    double thermalConductivity() override;

    void getMultiDiffCoeffs(const size_t ld, double* const d) override;

    //! Get the species diffusive mass fluxes wrt to the mass averaged velocity,
    //! given the gradients in mole fraction and temperature
    /*!
     * Units for the returned fluxes are kg m-2 s-1.
     *
     * @param ndim     Number of dimensions in the flux expressions
     * @param grad_T   Gradient of the temperature (length = ndim)
     * @param ldx      Leading dimension of the grad_X array. (usually equal to
     *                 m_nsp but not always)
     * @param grad_X   Gradients of the mole fraction. Flat vector with the
     *                 m_nsp in the inner loop. length = ldx * ndim
     * @param ldf      Leading dimension of the fluxes array. (usually equal to
     *                 m_nsp but not always)
     * @param fluxes   Output of the diffusive mass fluxes. Flat vector with the
     *                 m_nsp in the inner loop. length = ldx * ndim
     */
    void getSpeciesFluxes(size_t ndim, const double* const grad_T,
                          size_t ldx, const double* const grad_X,
                          size_t ldf, double* const fluxes) override;

    //! Get the molar diffusional fluxes [kmol/m^2/s] of the species, given the
    //! thermodynamic state at two nearby points.
    /*!
     * The molar diffusional fluxes are calculated with reference to the mass
     * averaged velocity. This is a one-dimensional vector
     *
     * @param state1 Array of temperature, density, and mass
     *               fractions for state 1.
     * @param state2 Array of temperature, density, and mass
     *               fractions for state 2.
     * @param delta  Distance from state 1 to state 2 (m).
     * @param fluxes Output molar fluxes of the species. (length = m_nsp)
     */
    void getMolarFluxes(const double* const state1, const double* const state2,
                        const double delta, double* const fluxes) override;

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
     * @param fluxes Output mass fluxes of the species. (length = m_nsp)
     */
    void getMassFluxes(const double* state1, const double* state2, double delta,
                       double* fluxes) override;

    void init(ThermoPhase* thermo, int mode=0, int log_level=-7) override;

    void invalidateCache() override;

protected:
    //! Update basic temperature-dependent quantities if the temperature has
    //! changed.
    void update_T() override;

    //! Update basic concentration-dependent quantities if the concentrations
    //! have changed.
    void update_C() override;

    //! Update the temperature-dependent terms needed to compute the thermal
    //! conductivity and thermal diffusion coefficients.
    void updateThermal_T();

    double m_thermal_tlast;

    //! Dense matrix for astar
    DenseMatrix m_astar;

    //! Dense matrix for bstar
    DenseMatrix m_bstar;

    //! Dense matrix for cstar
    DenseMatrix m_cstar;

    vector<double> m_cinternal;

    vector<double> m_sqrt_eps_k;
    DenseMatrix m_log_eps_k;
    vector<double> m_frot_298;
    vector<double> m_rotrelax;

    double m_lambda;

    // L matrix quantities
    DenseMatrix m_Lmatrix;
    DenseMatrix m_aa;
    vector<double> m_a;
    vector<double> m_b;

    // work space
    vector<double> m_spwork1, m_spwork2, m_spwork3;

    //! Mole fraction vector from last L-matrix evaluation
    vector<double> m_molefracs_last;

    //! Boolean indicating viscosity is up to date
    bool m_l0000_ok;
    bool m_lmatrix_soln_ok;

    //! Evaluate the L0000 matrices
    /*!
     *  Evaluate the upper-left block of the L matrix.
     *  @param x vector of species mole fractions
     */
    void eval_L0000(const double* const x);

    //! Evaluate the L0010 matrices
    /*!
     *  @param x vector of species mole fractions
     */
    void eval_L0010(const double* const x);

    //! Evaluate the L1000 matrices
    void eval_L1000();

    void eval_L0100();
    void eval_L0001();
    void eval_L1010(const double* x);
    void eval_L1001(const double* x);
    void eval_L0110();
    void eval_L0101(const double* x);
    bool hasInternalModes(size_t j);

    double pressure_ig();

    virtual void solveLMatrixEquation();
    bool m_debug;
};
}
#endif
