/**
 * @file Transport.h Headers for the Transport object, which is the virtual
 *     base class for all transport property evaluators and also includes the
 *     tranprops group definition (see @ref tranprops and @link
 *     Cantera::Transport Transport @endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

/**
 * @defgroup tranprops Transport Properties
 *
 * These classes provide transport properties, including diffusion coefficients,
 * thermal conductivity, and viscosity.
 */

#ifndef CT_TRANSPORT_H
#define CT_TRANSPORT_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class ThermoPhase;

/**
 * @addtogroup tranprops
 */
//!  @cond

const int CK_Mode = 10;

//!   @endcond

//! Base class for transport property managers.
/*!
 * All classes that compute transport properties for a single phase derive from
 * this class.  Class Transport is meant to be used as a base class only. It is
 * possible to instantiate it, but its methods throw exceptions if called.
 *
 * ## Relationship of the Transport class to the ThermoPhase Class
 *
 * This section describes how calculations are carried out within the Transport
 * class. The Transport class and derived classes of the the Transport class
 * necessarily use the ThermoPhase class to obtain the list of species and the
 * thermodynamic state of the phase.
 *
 * No state information is stored within Transport classes. Queries to the
 * underlying ThermoPhase object must be made to obtain the state of the system.
 *
 * An exception to this however is the state information concerning the the
 * gradients of variables. This information is not stored within the ThermoPhase
 * objects. It may be collected within the Transport objects. In fact, the
 * meaning of const operations within the Transport class refers to calculations
 * which do not change the state of the system nor the state of the first order
 * gradients of the system.
 *
 * When a const operation is evoked within the Transport class, it is also
 * implicitly assumed that the underlying state within the ThermoPhase object
 * has not changed its values.
 *
 * @todo Provide a general mechanism to store the gradients of state variables
 *        within the system.
 *
 * @ingroup tranprops
 */
class Transport
{
public:
    //!  Constructor.
    /*!
     * New transport managers should be created using TransportFactory, not by
     * calling the constructor directly.
     *
     * @see TransportFactory
     */
    Transport() = default;

    virtual ~Transport() {}

    // Transport objects are not copyable or assignable
    Transport(const Transport&) = delete;
    Transport& operator=(const Transport&) = delete;

    //! Identifies the model represented by this Transport object. Each derived class
    //! should override this method to return a meaningful identifier.
    //! @since New in %Cantera 3.0. The name returned by this method corresponds
    //!     to the canonical name used in the YAML input format.
    virtual string transportModel() const {
        return "none";
    }

    /**
     * Phase object. Every transport manager is designed to compute properties
     * for a specific phase of a mixture, which might be a liquid solution, a
     * gas mixture, a surface, etc. This method returns a reference to the
     * object representing the phase itself.
     */
    ThermoPhase& thermo() {
        return *m_thermo;
    }

    //! Check that the specified species index is in range. Throws an exception
    //! if k is greater than nSpecies()
    void checkSpeciesIndex(size_t k) const;

    //! Check that an array size is at least nSpecies(). Throws an exception if
    //! kk is less than nSpecies(). Used before calls which take an array
    //! pointer.
    void checkSpeciesArraySize(size_t kk) const;

    //! @name Transport Properties
    //! @{

    /**
     * The viscosity in Pa-s.
     */
    virtual double viscosity() {
        throw NotImplementedError("Transport::viscosity",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns the pure species viscosities
    /*!
     * The units are Pa-s and the length is the number of species
     *
     * @param visc   Vector of viscosities
     */
    virtual void getSpeciesViscosities(double* const visc) {
        throw NotImplementedError("Transport::getSpeciesViscosities",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! The bulk viscosity in Pa-s.
    /*!
     * The bulk viscosity is only non-zero in rare cases. Most transport
     * managers either overload this method to return zero, or do not implement
     * it, in which case an exception is thrown if called.
     */
    virtual double bulkViscosity() {
        throw NotImplementedError("Transport::bulkViscosity",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns the mixture thermal conductivity in W/m/K.
    /*!
     * Units are in W / m K  or equivalently kg m / s3 K
     *
     * @returns thermal conductivity in W/m/K.
     */
    virtual double thermalConductivity() {
        throw NotImplementedError("Transport::thermalConductivity",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! The electrical conductivity (Siemens/m).
    virtual double electricalConductivity() {
        throw NotImplementedError("Transport::electricalConductivity",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Get the Electrical mobilities (m^2/V/s).
    /*!
     * This function returns the mobilities. In some formulations this is equal
     * to the normal mobility multiplied by Faraday's constant.
     *
     * Frequently, but not always, the mobility is calculated from the diffusion
     * coefficient using the Einstein relation
     *
     * @f[
     *      \mu^e_k = \frac{F D_k}{R T}
     * @f]
     *
     * @param mobil_e  Returns the mobilities of the species in array @c
     *               mobil_e. The array must be dimensioned at least as large as
     *               the number of species.
     */
    virtual void getMobilities(double* const mobil_e) {
        throw NotImplementedError("Transport::getMobilities",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! @}

    //! Get the species diffusive mass fluxes wrt to the specified solution
    //! averaged velocity, given the gradients in mole fraction and temperature
    /*!
     * Units for the returned fluxes are kg m-2 s-1.
     *
     * Usually the specified solution average velocity is the mass averaged
     * velocity. This is changed in some subclasses, however.
     *
     * @param ndim       Number of dimensions in the flux expressions
     * @param grad_T     Gradient of the temperature (length = ndim)
     * @param ldx        Leading dimension of the grad_X array (usually equal to
     *                   m_nsp but not always)
     * @param grad_X     Gradients of the mole fraction Flat vector with the
     *                   m_nsp in the inner loop. length = ldx * ndim
     * @param ldf        Leading dimension of the fluxes array (usually equal to
     *                   m_nsp but not always)
     * @param fluxes     Output of the diffusive mass fluxes. Flat vector with
     *                   the m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim, const double* const grad_T,
                                  size_t ldx, const double* const grad_X,
                                  size_t ldf, double* const fluxes) {
        throw NotImplementedError("Transport::getSpeciesFluxes",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Get the molar fluxes [kmol/m^2/s], given the thermodynamic state at two
    //! nearby points.
    /*!
     * @param[in] state1 Array of temperature, density, and mass fractions for
     *               state 1.
     * @param[in] state2 Array of temperature, density, and mass fractions for
     *               state 2.
     * @param[in] delta  Distance from state 1 to state 2 (m).
     * @param[out] cfluxes Output array containing the diffusive molar fluxes of
     *               species from state1 to state2. This is a flat vector with
     *               m_nsp in the inner loop. length = ldx * ndim. Units are
     *               [kmol/m^2/s].
     */
    virtual void getMolarFluxes(const double* const state1,
                                const double* const state2, const double delta,
                                double* const cfluxes) {
        throw NotImplementedError("Transport::getMolarFluxes",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Get the mass fluxes [kg/m^2/s], given the thermodynamic state at two
    //! nearby points.
    /*!
     * @param[in] state1 Array of temperature, density, and mass
     *               fractions for state 1.
     * @param[in] state2 Array of temperature, density, and mass fractions for
     *               state 2.
     * @param[in] delta Distance from state 1 to state 2 (m).
     * @param[out] mfluxes Output array containing the diffusive mass fluxes of
     *               species from state1 to state2. This is a flat vector with
     *               m_nsp in the inner loop. length = ldx * ndim. Units are
     *               [kg/m^2/s].
     */
    virtual void getMassFluxes(const double* state1,
                               const double* state2, double delta,
                               double* mfluxes) {
        throw NotImplementedError("Transport::getMassFluxes",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return a vector of Thermal diffusion coefficients [kg/m/sec].
    /*!
     * The thermal diffusion coefficient @f$ D^T_k @f$ is defined so that the
     * diffusive mass flux of species *k* induced by the local temperature
     * gradient is given by the following formula:
     *
     * @f[
     *     M_k J_k = -D^T_k \nabla \ln T.
     * @f]
     *
     * The thermal diffusion coefficient can be either positive or negative.
     *
     * @param dt On return, dt will contain the species thermal diffusion
     *           coefficients.  Dimension dt at least as large as the number of
     *           species. Units are kg/m/s.
     */
    virtual void getThermalDiffCoeffs(double* const dt) {
        throw NotImplementedError("Transport::getThermalDiffCoeffs",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns the matrix of binary diffusion coefficients [m^2/s].
    /*!
     * @param[in] ld   Inner stride for writing the two dimension diffusion
     *                 coefficients into a one dimensional vector
     * @param[out] d   Diffusion coefficient matrix (must be at least m_k * m_k
     *                 in length.
     */
    virtual void getBinaryDiffCoeffs(const size_t ld, double* const d) {
        throw NotImplementedError("Transport::getBinaryDiffCoeffs",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the Multicomponent diffusion coefficients. Units: [m^2/s].
    /*!
     * If the transport manager implements a multicomponent diffusion
     * model, then this method returns the array of multicomponent
     * diffusion coefficients. Otherwise it throws an exception.
     *
     * @param[in] ld  The dimension of the inner loop of d (usually equal to m_nsp)
     * @param[out] d  flat vector of diffusion coefficients, fortran ordering.
     *            d[ld*j+i] is the D_ij diffusion coefficient (the diffusion
     *            coefficient for species i due to concentration gradients in
     *            species j). Units: m^2/s
     */
    virtual void getMultiDiffCoeffs(const size_t ld, double* const d) {
        throw NotImplementedError("Transport::getMultiDiffCoeffs",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns a vector of mixture averaged diffusion coefficients
    /**
     * Mixture-averaged diffusion coefficients [m^2/s].  If the transport
     * manager implements a mixture-averaged diffusion model, then this method
     * returns the array of mixture-averaged diffusion coefficients. Otherwise
     * it throws an exception.
     *
     * @param d  Return vector of mixture averaged diffusion coefficients
     *           Units = m2/s. Length = n_sp
     */
    virtual void getMixDiffCoeffs(double* const d) {
        throw NotImplementedError("Transport::getMixDiffCoeffs",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns a vector of mixture averaged diffusion coefficients
    virtual void getMixDiffCoeffsMole(double* const d) {
        throw NotImplementedError("Transport::getMixDiffCoeffsMole",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Returns a vector of mixture averaged diffusion coefficients
    virtual void getMixDiffCoeffsMass(double* const d) {
        throw NotImplementedError("Transport::getMixDiffCoeffsMass",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the polynomial fits to the viscosity of species i
    virtual void getViscosityPolynomial(size_t i, double* coeffs) const{
        throw NotImplementedError("Transport::getViscosityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the temperature fits of the heat conductivity of species i
    virtual void getConductivityPolynomial(size_t i, double* coeffs) const{
        throw NotImplementedError("Transport::getConductivityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the polynomial fits to the binary diffusivity of species pair (i, j)
    virtual void getBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs) const{
        throw NotImplementedError("Transport::getBinDiffusivityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the polynomial fits to the collision integral of species pair (i, j)
    virtual void getCollisionIntegralPolynomial(size_t i, size_t j,
                                                double* astar_coeffs,
                                                double* bstar_coeffs,
                                                double* cstar_coeffs) const{
        throw NotImplementedError("Transport::getCollisionIntegralPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Modify the polynomial fits to the viscosity of species i
    virtual void setViscosityPolynomial(size_t i, double* coeffs){
        throw NotImplementedError("Transport::setViscosityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Modify the temperature fits of the heat conductivity of species i
    virtual void setConductivityPolynomial(size_t i, double* coeffs){
        throw NotImplementedError("Transport::setConductivityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Modify the polynomial fits to the binary diffusivity of species pair (i, j)
    virtual void setBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs){
        throw NotImplementedError("Transport::setBinDiffusivityPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Modify the polynomial fits to the collision integral of species pair (i, j)
    virtual void setCollisionIntegralPolynomial(size_t i, size_t j,
                                                double* astar_coeffs,
                                                double* bstar_coeffs,
                                                double* cstar_coeffs, bool flag){
        throw NotImplementedError("Transport::setCollisionIntegralPolynomial",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! Return the parameters for a phase definition which are needed to
    //! reconstruct an identical object using the newTransport function. This
    //! excludes the individual species transport properties, which are handled
    //! separately.
    AnyMap parameters() const;

    //! Get error metrics about any functional fits calculated for pure species
    //! transport properties.
    //!
    //! See GasTransport::fitDiffCoeffs and GasTransport::fitProperties.
    //!
    //! @warning  This method is an experimental part of the %Cantera API and may be
    //!      changed or removed without notice.
    //! @since New in %Cantera 3.1.
    AnyMap fittingErrors() const { return m_fittingErrors; };

    //! @name Transport manager construction
    //!
    //! These methods are used during construction.
    //! @{

    //! Initialize a transport manager
    /*!
     * This routine sets up a transport manager. It calculates the collision
     * integrals and populates species-dependent data structures.
     *
     * @param thermo  Pointer to the ThermoPhase object
     * @param mode    Chemkin compatible mode or not. This alters the
     *                 specification of the collision integrals. defaults to no.
     * @param log_level Defaults to zero, no logging
     * @deprecated The `log_level` parameter is deprecated and will be removed after
     *     %Cantera 3.1.
     */
    virtual void init(ThermoPhase* thermo, int mode=0, int log_level=-7) {}

    //! Boolean indicating the form of the transport properties polynomial fits.
    //! Returns true if the Chemkin form is used.
    virtual bool CKMode() const {
        throw NotImplementedError("Transport::CK_Mode",
            "Not implemented for transport model '{}'.", transportModel());
    }

    //! @}

    //! Invalidate any cached values which are normally updated only when a
    //! change in state is detected
    //! @since New in %Cantera 3.1.
    virtual void invalidateCache() {}

protected:
    //! pointer to the object representing the phase
    ThermoPhase* m_thermo;

    //! Number of species
    size_t m_nsp = 0;

    //! Maximum errors associated with fitting pure species transport properties.
    AnyMap m_fittingErrors;
};

}

#endif
