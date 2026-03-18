/**
 *  @file MixTransport.h
 *    Headers for the MixTransport object, which models transport properties
 *    in ideal gas solutions using a mixture averaged approximation
 *    (see @ref tranprops and @link Cantera::MixTransport MixTransport @endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MIXTRAN_H
#define CT_MIXTRAN_H

#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{
//! Class MixTransport implements mixture-averaged transport properties for
//! ideal gas mixtures.
/*!
 * The model is based on that described in Kee, et al. @cite kee2003.
 *
 * Specific mixture-averaged formulas are implemented by:
 * - viscosity()
 * - thermalConductivity()
 * - getMixDiffCoeffs()
 * - getMixDiffCoeffsMole()
 * - getMixDiffCoeffsMass()
 * - getThermalDiffCoeffs()
 * - getMobilities()
 *
 * @ingroup tranprops
 */
class MixTransport : public GasTransport
{
public:
    //! Default constructor.
    MixTransport() = default;

    string transportModel() const override {
        return (m_mode == CK_Mode) ? "mixture-averaged-CK" : "mixture-averaged";
    }

    //! Return the thermal diffusion coefficients [kg/m/s]
    /*!
     * Model by S. Chapman and T.G. Cowling @cite chapman1970.
     * For more information about this implementation and its validation,
     * see T. Zirwes and A. Kronenburg @cite zirwes2025.
     *
     * The thermal diffusion coefficient of species @f$ k @f$ is computed from
     * @f[
     *      D_k^{T}= \frac{1}{2}\rho\frac{M_k}{\bar{M}}D_{mk}'\Theta_k
     * @f]
     * with
     * @f[
     *      \Theta_k=\frac{15}{2}\frac{\bar{M}^2}{\rho}\sum_i\left(\frac{1.2C_{ki}^*-1}{D_{ki}}\right)\left(\frac{Y_k\frac{\eta_i}{M_i}a_i-Y_i\frac{\eta_k}{M_k}a_k}{M_k+M_i}\right)
     * @f]
     * where @f$ C_{k,i}^* @f$ is a reduced collision integral and
     * @f[
     *      a_k=\left(1+\frac{1.065}{2\sqrt{2}X_k}\sum_{i\ne k}X_i\Phi_{k,i}\right)^{-1},
     * @f]
     * with @f$ \Phi_{k,i} @f$ the Wilke mixing operator. The thermodiffusion
     * coefficients are then normalized with
     * @f[
     *      \hat{D}^T_k=D^T_k-Y_k\sum_i D^T_i.
     * @f]
     * This ensures that the sum of all thermodiffusion coefficients
     * and thus the sum of all thermodiffusion fluxes are zero.
     *
     * @param[out] dt  Vector of thermal diffusion coefficients
     */
    void getThermalDiffCoeffs(span<double> dt) override;

    //! Returns the mixture thermal conductivity [W/m/K]
    /*!
     * The thermal conductivity is computed from the following mixture rule:
     * @f[
     *     \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
     * @f]
     *
     * It's used to compute the flux of energy due to a thermal gradient
     *
     * @f[
     *     \mathbf{q} =  - \lambda \nabla T
     * @f]
     */
    double thermalConductivity() override;

    //! Get the electrical mobilities [m²/V/s]
    /*!
     * This function returns the mobilities. Here, the mobility is calculated from the
     * diffusion coefficient using the Einstein relation:
     *
     * @f[
     *     \mu^e_k = \frac{F D_{km}'}{R T}
     * @f]
     *
     * @param mobil  Returns the mobilities of the species in array @c mobil.
     *               The array must be dimensioned at least as large as the
     *               number of species.
     */
    void getMobilities(span<double> mobil) override;

    //! Update the internal parameters whenever the temperature has changed
    /*!
     * This is called whenever a transport property is requested if the
     * temperature has changed since the last call to update_T().
     */
    void update_T() override;

    //! Update the internal parameters whenever the concentrations have changed
    /*!
     * This is called whenever a transport property is requested if the
     * concentrations have changed since the last call to update_C().
     */
    void update_C() override;

    //! Get the species diffusive mass fluxes [kg/m²/s] with respect to the mass
    //! averaged velocity, given the gradients in mole fraction and temperature.
    /*!
     * The diffusive mass flux of species @e k is computed from
     * @f[
     *     \mathbf{j}_k = -\rho \frac{M_k}{\overline{M}} D_{km}' \nabla X_k.
     * @f]
     *
     * @param ndim  Number of dimensions in the flux expressions
     * @param[in] grad_T  Gradient of the temperature (length `ndim`)
     * @param ldx  Leading dimension of the `grad_X` array (usually equal to the number
     *     of species)
     * @param[in] grad_X  Gradients of the mole fractions; flattened matrix such that
     *     @f$ dX_k/dx_n = \tt{ grad\_X[n*ldx+k]} @f$ is the gradient of species *k*
     *     in dimension *n*. Length is `ldx` * `ndim`.
     * @param ldf  Leading dimension of the `fluxes` array (usually equal to the number
     *     of species)
     * @param[out] fluxes  The diffusive mass fluxes; flattened matrix such that
     *     @f$ j_{kn} = \tt{ fluxes[n*ldf+k]} @f$ is the flux of species *k*
     *     in dimension *n*. Length is `ldf` * `ndim`.
     */
    void getSpeciesFluxes(size_t ndim, span<const double> grad_T,
                          size_t ldx, span<const double> grad_X,
                          size_t ldf, span<double> fluxes) override;

    void init(shared_ptr<ThermoPhase> thermo, int mode=0) override;

protected:
    //! Update the temperature dependent parts of the species thermal
    //! conductivities
    /*!
     * These are evaluated from the polynomial fits of the temperature and are
     * assumed to be independent of pressure
     */
    void updateCond_T();

    //! vector of species thermal conductivities [W/m/K]
    /*!
     * These are used in Wilke's rule to calculate the viscosity of the
     * solution. length = #m_nsp.
     */
    vector<double> m_cond;

    //! Internal storage for the calculated mixture thermal conductivity [W/m/K]
    double m_lambda = 0.0;

    //! Update boolean for the species thermal conductivities
    bool m_spcond_ok = false;

    //! Update boolean for the mixture rule for the mixture thermal conductivity
    bool m_condmix_ok = false;
};
}
#endif
