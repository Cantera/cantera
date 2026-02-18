/**
 *  @file UnityLewisTransport.h
 *    Headers for the UnityLewisTransport object, which models transport
 *    properties in ideal gas solutions using the unity Lewis number
 *    approximation
 *    (see @ref tranprops and @link Cantera::UnityLewisTransport UnityLewisTransport @endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITYLEWISTRAN_H
#define CT_UNITYLEWISTRAN_H

#include "MixTransport.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{
//! Class UnityLewisTransport implements the unity Lewis number approximation
//! for the mixture-averaged species diffusion coefficients. Mixture-averaged
//! transport properties for viscosity and thermal conductivity are inherited
//! from the MixTransport class.
//! @ingroup tranprops
class UnityLewisTransport : public MixTransport
{
public:
    UnityLewisTransport() = default;

    string transportModel() const override {
        return "unity-Lewis-number";
    }

    //! Returns the unity Lewis number approximation based diffusion
    //! coefficients [m²/s].
    /*!
     * Returns the unity Lewis number approximation based diffusion coefficients
     * for a gas, appropriate for calculating the mass averaged diffusive flux
     * with respect to the mass averaged velocity using gradients of the mole
     * fraction.
     *
     * @f[
     *     D^\prime_{km} = \frac{\lambda}{\rho c_p}
     * @f]
     *
     * In order to obtain the expected behavior from a unity Lewis number model,
     * this formulation requires that the correction velocity be computed as
     *
     * @f[
     *     V_c = \sum \frac{W_k}{\overline{W}} D^\prime_{km} \nabla X_k
     * @f]
     *
     * @param[out] d  Vector of diffusion coefficients for each species. length #m_nsp.
     */
    void getMixDiffCoeffs(span<double> d) override {
        double Dm = thermalConductivity() / (m_thermo->density() * m_thermo->cp_mass());
        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm;
        }
    }

    //! Thermal diffusion is not enabled in the unity Lewis number model.
    /*!
     * @param[out] dt  Thermal diffusion coefficients all set to zero.
     */
    void getThermalDiffCoeffs(span<double> dt) override {
        for (size_t k = 0; k < m_nsp; k++) {
            dt[k] = 0.0;
        }
    }

    //! Not implemented for unity Lewis number approximation
    void getMixDiffCoeffsMole(span<double> d) override {
        throw NotImplementedError("UnityLewisTransport::getMixDiffCoeffsMole");
    }

    //! Returns the unity Lewis number approximation based diffusion
    //! coefficients [m²/s].
    /*!
     * These are the coefficients for calculating the diffusive mass fluxes
     * from the species mass fraction gradients, computed as
     *
     * @f[
     *     D_{km} = \frac{\lambda}{\rho c_p}
     * @f]
     *
     * @param[out] d  Vector of diffusion coefficients for each species; length #m_nsp.
     */
    void getMixDiffCoeffsMass(span<double> d) override {
        double Dm = thermalConductivity() / (m_thermo->density() * m_thermo->cp_mass());
        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm;
        }
    }
};
}
#endif
