/**
 * @file DustyGasTransport.h Headers for the DustyGasTransport object, which
 *   models transport properties in porous media using the dusty gas
 *   approximation (see @ref tranprops and @link Cantera::DustyGasTransport
 *   DustyGasTransport @endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DUSTYGASTRAN_H
#define CT_DUSTYGASTRAN_H

// Cantera includes
#include "Transport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{
//! Class DustyGasTransport implements the Dusty Gas model for transport in porous media.
/*!
 * As implemented here, only species transport is handled. The viscosity,
 * thermal conductivity, and thermal diffusion coefficients are not implemented.
 *
 * The dusty gas model includes the effects of Darcy's law. There is a net flux
 * of species due to a pressure gradient that is part of Darcy's law.
 *
 * The dusty gas model expresses the value of the molar flux of species
 * @f$ k @f$, @f$ J_k @f$ by the following formula.
 *
 * @f[
 *     \sum_{j \ne k}{\frac{X_j J_k - X_k J_j}{D^e_{kj}}} + \frac{J_k}{\mathcal{D}^{e}_{k,knud}} =
 *             - \nabla C_k  - \frac{C_k}{\mathcal{D}^{e}_{k,knud}} \frac{\kappa}{\mu} \nabla p
 * @f]
 *
 * @f$ j @f$ is a sum over all species in the gas.
 *
 * The effective Knudsen diffusion coefficients are given by the following form
 *
 * @f[
 *    \mathcal{D}^e_{k,knud} =  \frac{2}{3} \frac{r_{pore} \phi}{\tau} \left( \frac{8 R T}{\pi W_k}  \right)^{1/2}
 * @f]
 *
 * The effective knudsen diffusion coefficients take into account the effects of
 * collisions of gas-phase molecules with the wall.
 *
 * For references on the Dusty Gas Model, see Zhu and Kee @cite zhu2006; Zhu, et al.
 * @cite zhu2005; Mason and Malinauskas @cite mason1983; and Veldsink, et al.
 * @cite veldsink1995.
 * @ingroup tranprops
 */
class DustyGasTransport : public Transport
{
public:
    //! default constructor
    DustyGasTransport() = default;

    //  overloaded base class methods

    string transportModel() const override {
        return "DustyGas";
    }

    void getMultiDiffCoeffs(const size_t ld, double* const d) override;

    //! Get the molar fluxes [kmol/m^2/s], given the thermodynamic state at two nearby points.
    /*!
     *   @f[
     *       J_k = - \sum_{j = 1, N} \left[D^{multi}_{kj}\right]^{-1} \left( \nabla C_j  + \frac{C_j}{\mathcal{D}^{knud}_j} \frac{\kappa}{\mu} \nabla p \right)
     *   @f]
     *
     * @param  state1  Array of temperature, density, and mass fractions for state 1.
     * @param  state2  Array of temperature, density, and mass fractions for state 2.
     * @param  delta   Distance from state 1 to state 2 (m).
     *
     * @param fluxes   Vector of species molar fluxes due to diffusional driving force
     */
    void getMolarFluxes(const double* const state1, const double* const state2,
                        const double delta, double* const fluxes) override;

    // new methods added in this class

    //! Set the porosity (dimensionless)
    /*!
     * @param porosity  Set the value of the porosity
     */
    void setPorosity(double porosity);

    //! Set the tortuosity (dimensionless)
    /*!
     * Tortuosity is considered to be constant within the object
     *
     * @param tort  Value of the tortuosity
     */
    void setTortuosity(double tort);

    //! Set the mean pore radius (m)
    /*!
     * @param rbar  Value of the pore radius ( m)
     */
    void setMeanPoreRadius(double rbar);

    //! Set the mean particle diameter
    /*!
     * @param dbar  Set the mean particle diameter (m)
     */
    void setMeanParticleDiameter(double dbar);

    //! Set the permeability of the media
    /*!
     * If not set, the value for close-packed spheres will be used by default.
     *
     * The value for close-packed spheres is given below, where p is the
     * porosity, t is the tortuosity, and d is the diameter of the sphere
     *
     * @f[
     *     \kappa = \frac{p^3 d^2}{72 t (1 - p)^2}
     * @f]
     *
     * @param B  set the permeability of the media (units = m^2)
     */
    void setPermeability(double B);

    //! Return a reference to the transport manager used to compute the gas
    //! binary diffusion coefficients and the viscosity.
    /*!
     * @returns a reference to the gas transport object
     */
    Transport& gasTransport();

    //! Make the TransportFactory object a friend, because this object has
    //! restricted its instantiation to classes which are friends.
    friend class TransportFactory;

protected:
    //! Initialization routine called by TransportFactory
    /*!
     * The DustyGas model is a subordinate model to the gas phase transport
     * model. Here we set the gas phase models.
     *
     * This is a protected routine, so that initialization of the Model must
     * occur within Cantera's setup
     *
     * @param  phase  Pointer to the underlying ThermoPhase model for the gas phase
     * @param  gastr  Pointer to the underlying Transport model for transport in
     *     the gas phase.
     */
    void initialize(ThermoPhase* phase, Transport* gastr);

private:
    //! Update temperature-dependent quantities within the object
    /*!
     * The object keeps a value m_temp, which is the temperature at which
     * quantities were last evaluated at. If the temperature is changed, update
     * Booleans are set false, triggering recomputation.
     */
    void updateTransport_T();

    //! Update concentration-dependent quantities within the object
    /*!
     * The object keeps a value m_temp, which is the temperature at which
     * quantities were last evaluated at. If the temperature is changed, update
     * Booleans are set false, triggering recomputation.
     */
    void updateTransport_C();

    //! Private routine to update the dusty gas binary diffusion coefficients
    /*!
     * The dusty gas binary diffusion coefficients @f$  D^{dg}_{i,j} @f$ are
     * evaluated from the binary gas-phase diffusion coefficients @f$
     * D^{bin}_{i,j} @f$  using the following formula
     *
     * @f[
     *     D^{dg}_{i,j} =  \frac{\phi}{\tau} D^{bin}_{i,j}
     * @f]
     *
     * where @f$ \phi @f$ is the porosity of the media and @f$ \tau @f$ is the
     * tortuosity of the media.
     */
    void updateBinaryDiffCoeffs();

    //! Update the Multicomponent diffusion coefficients that are used in the
    //! approximation
    /*!
     * This routine updates the H matrix and then inverts it.
     */
    void updateMultiDiffCoeffs();

    //! Update the Knudsen diffusion coefficients
    /*!
     * The Knudsen diffusion coefficients are given by the following form
     *
     * @f[
     *     \mathcal{D}^{knud}_k =  \frac{2}{3} \frac{r_{pore} \phi}{\tau} \left( \frac{8 R T}{\pi W_k}  \right)^{1/2}
     * @f]
     */
    void updateKnudsenDiffCoeffs();

    //! Calculate the H matrix
    /*!
     * The multicomponent diffusion H matrix @f$  H_{k,l} @f$ is given by the following form
     *
     * @f[
     *    H_{k,l} = - \frac{X_k}{D_{k,l}}
     * @f]
     * @f[
     *    H_{k,k} = \frac{1}{\mathcal(D)^{knud}_{k}} + \sum_{j \ne k}^N{ \frac{X_j}{D_{k,j}} }
     * @f]
     */
    void eval_H_matrix();

    //! Local copy of the species molecular weights
    /*!
     *  units kg /kmol
     *  length = m_nsp;
     */
    vector<double> m_mw;

    //! binary diffusion coefficients
    DenseMatrix m_d;

    //! mole fractions
    vector<double> m_x;

    //! Knudsen diffusion coefficients. @see updateKnudsenDiffCoeffs()
    vector<double> m_dk;

    //! temperature
    double m_temp = -1.0;

    //! Multicomponent diffusion coefficients. @see eval_H_matrix()
    DenseMatrix m_multidiff;

    //! work space of size m_nsp;
    vector<double> m_spwork;

    //! work space of size m_nsp;
    vector<double> m_spwork2;

    //! Pressure Gradient
    double m_gradP = 0.0;

    //! Update-to-date variable for Knudsen diffusion coefficients
    bool m_knudsen_ok = false;

    //! Update-to-date variable for Binary diffusion coefficients
    bool m_bulk_ok = false;

    //! Porosity
    double m_porosity = 0.0;

    //! Tortuosity
    double m_tortuosity = 1.0;

    //! Pore radius (meter)
    double m_pore_radius = 0.0;

    //! Particle diameter
    /*!
     * The medium is assumed to consist of particles of size m_diam. units =  m
     */
    double m_diam = 0.0;

    //! Permeability of the media
    /*!
     * The permeability is the proportionality constant for Darcy's law which
     * relates discharge rate and viscosity to the applied pressure gradient.
     *
     * Below is Darcy's law, where @f$ \kappa @f$ is the permeability
     *
     * @f[
     *     v = \frac{\kappa}{\mu} \frac{\delta P}{\delta x}
     * @f]
     *
     * units are m2
     */
    double m_perm = -1.0;

    //! Pointer to the transport object for the gas phase
    unique_ptr<Transport> m_gastran;
};
}
#endif
