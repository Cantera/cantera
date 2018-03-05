/**
 *  @file UnityLewisTransport.h
 *    Headers for the UnityLewisTransport object, which models transport
 *    properties in ideal gas solutions using the unity Lewis number
 *    approximation
 *    (see \ref tranprops and \link Cantera::UnityLewisTransport UnityLewisTransport \endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITYLEWISTRAN_H
#define CT_UNITYLEWISTRAN_H

#include "MixTransport.h"

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
//    UnityLewisTransport() {}

    virtual std::string transportType() const {
        return "UnityLewis";
    }

    //! Returns the unity Lewis number approximation based diffusion
    //! coefficients [m^2/s].
    /*!
     * Returns the unity Lewis number approximation based diffusion coefficients
     * for a gas, appropriate for calculating the mass averaged diffusive flux
     * with respect to the mass averaged velocity using gradients of the mole
     * fraction.
     *
     * \f[
     *     D^\prime_{km} = \frac{\lambda}{\rho c_p}
     * \f]
     *
     * @param[out] d  Vector of diffusion coefficients for each species (m^2/s).
     * length m_nsp.
     */
    virtual void getMixDiffCoeffs(double* const d) {
        double Dm = thermalConductivity() / (m_thermo->density() * m_thermo->cp_mass());
        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm;
        }
    }

    //! Not implemented for unity Lewis number approximation
    virtual void getMixDiffCoeffsMole(double* const d){
        throw NotImplementedError("UnityLewisTransport::getMixDiffCoeffsMole");
    }

    //! Not implemented for unity Lewis number approximation
    virtual void getMixDiffCoeffsMass(double* const d){
        throw NotImplementedError("UnityLewisTransport::getMixDiffCoeffsMass");
    }
};
}
#endif
