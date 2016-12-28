/**
 *  @file LiquidTransport.h
 *   Header file defining class LiquidTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_LIQUIDTRAN_H
#define CT_LIQUIDTRAN_H

#include "TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "LiquidTransportParams.h"

namespace Cantera
{
//! Class LiquidTransport implements models for transport
//! properties for liquid phases.
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * Liquid Transport is set up with some flexibility in this class.  Transport
 * properties like viscosity and thermal conductivity are allowed flexibility
 * within the constraints of the LiquidTransportProperty and
 * LiquidTransportInteractions classes. For species diffusion, the
 * LiquidTransport class focuses on the Stefan-Maxwell equation to determine the
 * diffusion velocities.  Other options for liquid diffusion include solvent-
 * dominated diffusion, and a class SolventTransport should be forthcoming.
 *
 * The class LiquidTransport has several roles.
 * -# It brings together the individual species transport properties, expressed
 *    as subclasses of LTPspecies (Liquid Transport Properties of Species)
 *    through LiquidTransportData, with models for the composition dependence of
 *    liquid transport properties expressed as subclasses of
 *    LiquidTranInteraction (mixing rules) through LiquidTransportParams.
 *    Calculating mixture properties generally consists of calling the
 *    getMixTansProp member of LiquidTranInteraction by passing a vector of
 *    LTPSpecies
 * -# It calculates the bulk velocity \f$ \vec{v} \f$ and individual species
 *    diffusion velocities, \f$ \vec{V_i} \f$ using the Stefan-Maxwell
 *    equations.  It is possible to set a flag to calculate relative to a mass-
 *    averaged bulk velocity, relative to a mole-averaged bulk velocity or
 *    relative to a single species velocity using the \<velocityBasis
 *    basis="mass"\>, \<velocityBasis basis="mass"\>, or \<velocityBasis
 *    basis="Cl-"\> keyword. Mass-averaged velocities are the default for which
 *    the diffusion velocities satisfy
 *     \f[
 *        \sum_{i} Y_i \vec{V_i} = 0
 *     \f]
 *     for mass fraction \f$ Y_i \f$.  For mole-averaged velocities
 *     \f[
 *        \sum_{i} X_i \vec{V_i} = 0
 *     \f]
 *     for mole fraction \f$ X_i \f$. or
 *     \f[
 *        \vec{V_i} = 0
 *     \f]
 *     for reference species \f$ i \f$.
 * -# It provides access to a number of derived quantities related to transport
 *    properties as described in the various methods below.
 *
 * Within LiquidTransport, the state is presumed to be defined in terms of the
 * species mole fraction, temperature and pressure.  Charged species are
 * expected and quantities like the electric current are computed based on a
 * combined electrochemical potential.
 *
 * @ingroup tranprops
 */
class LiquidTransport : public Transport
{
public:
    //! Default constructor.
    /*!
     * This requires call to initLiquid(LiquidTransportParams& tr)
     * after filling LiquidTransportParams to complete instantiation.
     * The filling of LiquidTransportParams is currently carried out
     * in the TransportFactory class, but might be moved at some point.
     *
     * @param thermo  ThermoPhase object holding species information.
     * @param ndim    Number of spatial dimensions.
     */
    LiquidTransport(thermo_t* thermo = 0, int ndim = 1);

    LiquidTransport(const LiquidTransport& right);
    LiquidTransport& operator=(const LiquidTransport& right);
    virtual Transport* duplMyselfAsTransport() const;
    virtual ~LiquidTransport();

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient. We get
     * the object ready to do property evaluations. A lot of the input required
     * to do property evaluations is contained in the LiquidTransportParams
     * class that is filled in TransportFactory.
     *
     * @param tr  Transport parameters for all of the species in the phase.
     */
    virtual bool initLiquid(LiquidTransportParams& tr);

    friend class TransportFactory;

    virtual int model() const {
        warn_deprecated("LiquidTransport::model",
                        "To be removed after Cantera 2.3.");
        return cLiquidTransport;
    }

    virtual std::string transportType() const {
        return "Liquid";
    }

    //! Returns the viscosity of the solution
    /*!
     * The viscosity calculation is handled by subclasses of
     * LiquidTranInteraction as specified in the input file. These in turn
     * employ subclasses of LTPspecies to determine the individual species
     * viscosities.
     */
    virtual doublereal viscosity();

    //! Returns the pure species viscosities for all species
    /*!
     * The pure species viscosities are evaluated using the appropriate
     * subclasses of LTPspecies as specified in the input file.
     *
     * @param visc  array of length "number of species"
     *              to hold returned viscosities.
     */
    virtual void getSpeciesViscosities(doublereal* const visc);

    //! Returns the ionic conductivity of the solution
    /*!
     * The ionic conductivity calculation is handled by subclasses of
     * LiquidTranInteraction as specified in the input file. These in turn
     * employ subclasses of LTPspecies to determine the individual species ionic
     * conductivities.
     */
    virtual doublereal ionConductivity();

    //! Returns the pure species ionic conductivities for all species
    /*!
     * The pure species ionic conductivities are evaluated using the appropriate
     * subclasses of LTPspecies as specified in the input file.
     *
     * @param ionCond  Array of length "number of species" to hold returned
     *     ionic conductivities.
     */
    virtual void getSpeciesIonConductivity(doublereal* const ionCond);

    //! Returns the pointer to the mobility ratios of the binary
    //! combinations of the transported species for the solution
    //! Has size of the number of binary interactions = nsp*nsp
    /*!
     * The mobility ratio calculation is handled by subclasses of
     * LiquidTranInteraction as specified in the input file. These in turn
     * employ subclasses of LTPspecies to determine the mobility ratios in the
     * pure species.
     *
     * @param mobRat  Vector of mobility ratios
     */
    virtual void mobilityRatio(doublereal* mobRat);

    //! Returns a double pointer to the mobility ratios of the transported
    //! species in each pure species phase.
    /*!
     * Has size of the number of binary interactions by the number of species
     * (nsp*nsp X nsp). The pure species mobility ratios are evaluated using the
     * appropriate subclasses of LTPspecies as specified in the input file.
     *
     * @param mobRat  array of length "number of species" to hold returned
     *                mobility ratios.
     */
    virtual void getSpeciesMobilityRatio(doublereal** mobRat);

    //! Returns the self diffusion coefficients of the species in the phase.
    //! Has size of nsp(coeffs)
    /*!
     * The self diffusion coefficient is the diffusion coefficient of a tracer
     * species at the current temperature and composition of the species.
     * Therefore, the dilute limit of transport is assumed for the tracer
     * species. The effective formula may be calculated from the Stefan-Maxwell
     * formulation by adding another row for the tracer species, assigning all
     * D's to be equal to the respective species D's, and then taking the limit
     * as the tracer species mole fraction goes to zero. The corresponding flux
     * equation for the tracer species k in units of kmol m-2 s-1 is.
     *
     * \f[
     *      J_k = - D^{sd}_k \frac{C_k}{R T}  \nabla \mu_k
     * \f]
     *
     * The derivative is taken at constant T and P.
     *
     * The self diffusion calculation is handled by subclasses of
     * LiquidTranInteraction as specified in the input file. These in turn
     * employ subclasses of LTPspecies to determine the individual species self
     * diffusion coeffs.
     *
     * @param selfDiff Vector of self-diffusion coefficients. Length = number of
     *                  species in phase. units = m**2 s-1
     */
    virtual void selfDiffusion(doublereal* const selfDiff);

    //! Returns the self diffusion coefficients in the pure species phases.
    //! Has size of nsp(coeffs) x nsp(phases)
    /*!
     * The pure species molar volumes are evaluated using the appropriate
     * subclasses of LTPspecies as specified in the input file.
     *
     * @param selfDiff  array of length "number of species" to hold returned
     *              self diffusion coeffs.
     */
    virtual void getSpeciesSelfDiffusion(doublereal** selfDiff);

    //! Returns the hydrodynamic radius for all species
    /*!
     * The species hydrodynamic radii are evaluated using the appropriate
     * subclasses of LTPspecies as specified in the input file.
     *
     * @param radius  array of length "number of species" to hold returned
     *                radii.
     */
    virtual void getSpeciesHydrodynamicRadius(doublereal* const radius);

    //! Returns the binary diffusion coefficients
    /*!
     * The binary diffusion coefficients are specified in the input file through
     * the LiquidTransportInteractions class.  These are the binary interaction
     * coefficients employed in the Stefan-Maxwell equation.
     *
     * @param ld  number of species in system
     * @param d   vector of binary diffusion coefficients. units = m2 s-1.
     *          length = ld*ld = (number of species)^2
     */
    virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);

    //! Get the Mixture diffusion coefficients
    /*!
     * The mixture diffusion coefficients are not well defined in the context of
     * LiquidTransport because the Stefan Maxwell equation is solved.  Here the
     * mixture diffusion coefficients are defined according to Ficks law:
     * \f[
     *     X_i \vec{V_i} = -D_i \nabla X_i.
     * \f]
     * Solving Ficks Law for \f$ D_i \f$ gives a mixture diffusion coefficient
     * \f[
     *     D_i = - X_i \vec{V_i} / ( \nabla X_i ).
     * \f]
     * If \f$ \nabla X_i = 0 \f$ this is undefined and the nonsensical value -1
     * is returned.
     *
     * Note that this evaluation of \f$ \vec{V_i} \f$ requires a solve of the
     * Stefan Maxwell equation making this determination of the mixture averaged
     * diffusion coefficients a \e slow method for obtaining diffusion
     * coefficients.
     *
     * Also note that the Stefan Maxwell solve will be based upon the
     * thermodynamic state (including gradients) most recently set.  Gradients
     * can be set specifically using set_Grad_V, set_Grad_X and set_Grad_T or
     * through calls to getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff,
     * getSpeciesVdiffES, etc.
     *
     * @param d vector of mixture diffusion coefficients. units = m2 s-1.
     *          length = number of species
     */
    virtual void getMixDiffCoeffs(doublereal* const d);

    //! Return the thermal diffusion coefficients
    /*!
     * These are all zero for this simple implementation
     *
     * @param dt thermal diffusion coefficients
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Return the thermal conductivity of the solution
    /*!
     * The thermal conductivity calculation is handled by subclasses of
     * LiquidTranInteraction as specified in the input file. These in turn
     * employ subclasses of LTPspecies to determine the individual species
     * thermal conductivities.
     */
    virtual doublereal thermalConductivity();

    //! Get the Electrical mobilities (m^2/V/s).
    /*!
     * The electrical mobilities are not well defined in the context of
     * LiquidTransport because the Stefan Maxwell equation is solved.  Here the
     * electrical mobilities are calculated from the mixture-averaged diffusion
     * coefficients through a call to getMixDiffCoeffs() using the Einstein
     * relation
     *
     *     \f[
     *          \mu^e_k = \frac{F D_k}{R T}
     *     \f]
     *
     * Note that this call to getMixDiffCoeffs() requires a solve of the Stefan
     * Maxwell equation making this determination of the mixture averaged
     * diffusion coefficients a \e slow method for obtaining diffusion
     * coefficients.
     *
     * Also note that the Stefan Maxwell solve will be based upon the
     * thermodynamic state (including gradients) most recently set. Gradients
     * can be set specifically using set_Grad_V, set_Grad_X and set_Grad_T or
     * through calls to getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff,
     * getSpeciesVdiffES, etc.
     *
     * @param mobil_e  Returns the electrical mobilities of the species in array
     *                 \c mobil_e. The array must be dimensioned at least as
     *                 large as the number of species.
     */
    virtual void getMobilities(doublereal* const mobil_e);

    //! Get the fluid mobilities (s kmol/kg).
    /*!
     * The fluid mobilities are not well defined in the context of
     * LiquidTransport because the Stefan Maxwell equation is solved.  Here the
     * fluid mobilities are calculated from the mixture-averaged diffusion
     * coefficients through a call to getMixDiffCoeffs() using the Einstein
     * relation
     *
     * \f[
     *     \mu^f_k = \frac{D_k}{R T}
     * \f]
     *
     * Note that this call to getMixDiffCoeffs() requires a solve of the Stefan
     * Maxwell equation making this determination of the mixture averaged
     * diffusion coefficients a \e slow method for obtaining diffusion
     * coefficients.
     *
     * Also note that the Stefan Maxwell solve will be based upon the
     * thermodynamic state (including gradients) most recently set. Gradients
     * can be set specifically using set_Grad_V, set_Grad_X and set_Grad_T or
     * through calls to getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff,
     * getSpeciesVdiffES, etc.
     *
     * @param mobil_f  Returns the fluid mobilities of the species in array \c
     *                 mobil_f. The array must be dimensioned at least as large
     *                 as the number of species.
     */
    virtual void getFluidMobilities(doublereal* const mobil_f);

    //! Specify the value of the gradient of the voltage
    /*!
     * @param grad_V Gradient of the voltage (length num dimensions);
     */
    virtual void set_Grad_V(const doublereal* const grad_V);

    //! Specify the value of the gradient of the temperature
    /*!
     * @param grad_T Gradient of the temperature (length num dimensions);
     */
    virtual void set_Grad_T(const doublereal* const grad_T);

    //! Specify the value of the gradient of the MoleFractions
    /*!
     * @param grad_X Gradient of the mole fractions(length nsp * num dimensions);
     */
    virtual void set_Grad_X(const doublereal* const grad_X);

    //! Compute the mixture electrical conductivity from
    //! the Stefan-Maxwell equation.
    /*!
     * To compute the mixture electrical conductance, the Stefan Maxwell
     * equation is solved for zero species gradients and for unit potential
     * gradient, \f$ \nabla V \f$. The species fluxes are converted to current
     * by summing over the charge-weighted fluxes according to
     * \f[
     *     \vec{i} = \sum_{i} z_i F \rho \vec{V_i} / W_i
     * \f]
     * where \f$ z_i \f$ is the charge on species i, \f$ F \f$ is Faraday's
     * constant, \f$ \rho \f$  is the density, \f$ W_i \f$ is the molecular mass
     * of species i. The conductance, \f$ \kappa \f$ is obtained from
     * \f[
     *     \kappa = \vec{i} / \nabla V.
     * \f]
     */
    virtual doublereal getElectricConduct();

    //! Compute the electric current density in A/m^2
    /*!
     * The electric current is computed first by computing the species diffusive
     * fluxes using the Stefan Maxwell solution and then the current, \f$
     * \vec{i} \f$ by summing over the charge-weighted fluxes according to
     * \f[
     *     \vec{i} = \sum_{i} z_i F \rho \vec{V_i} / W_i
     * \f]
     * where \f$ z_i \f$ is the charge on species i, \f$ F \f$ is Faraday's
     * constant, \f$ \rho \f$  is the density, \f$ W_i \f$ is the molecular mass
     * of species \c i.
     *
     * @param ndim       The number of spatial dimensions (1, 2, or 3).
     * @param grad_T     The temperature gradient (ignored in this model).
     * @param ldx        Leading dimension of the grad_X array.
     * @param grad_X     Gradients of the mole fraction. Flat vector with the
     *                   m_nsp in the inner loop. length = ldx * ndim
     * @param ldf        Leading dimension of the grad_V and current vectors.
     * @param grad_V     The electrostatic potential gradient.
     * @param current    The electric current in A/m^2.
     */
    virtual void getElectricCurrent(int ndim,
                                    const doublereal* grad_T,
                                    int ldx,
                                    const doublereal* grad_X,
                                    int ldf,
                                    const doublereal* grad_V,
                                    doublereal* current);

    //! Get the species diffusive velocities wrt to the averaged velocity,
    //! given the gradients in mole fraction and temperature
    /*!
     * The average velocity can be computed on a mole-weighted
     * or mass-weighted basis, or the diffusion velocities may
     * be specified as relative to a specific species (i.e. a
     * solvent) all according to the velocityBasis input parameter.
     *
     * Units for the returned velocities are m s-1.
     *
     * @param ndim Number of dimensions in the flux expressions
     * @param grad_T Gradient of the temperature (length = ndim)
     * @param ldx  Leading dimension of the grad_X array (usually equal to m_nsp
     *              but not always)
     * @param grad_X Gradients of the mole fraction. Flat vector with the m_nsp
     *             in the inner loop. length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array (usually equal to m_nsp
     *              but not always)
     * @param Vdiff  Output of the diffusive velocities. Flat vector with the
     *             m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesVdiff(size_t ndim,
                                 const doublereal* grad_T,
                                 int ldx,
                                 const doublereal* grad_X,
                                 int ldf,
                                 doublereal* Vdiff);

    //! Get the species diffusive velocities wrt to the averaged velocity, given
    //! the gradients in mole fraction, temperature and electrostatic potential.
    /*!
     * The average velocity can be computed on a mole-weighted or mass-weighted
     * basis, or the diffusion velocities may be specified as relative to a
     * specific species (i.e. a solvent) all according to the velocityBasis
     * input parameter.
     *
     * Units for the returned velocities are m s-1.
     *
     * @param ndim       Number of dimensions in the flux expressions
     * @param grad_T     Gradient of the temperature (length = ndim)
     * @param ldx         Leading dimension of the grad_X array (usually equal
     *                       to m_nsp but not always)
     * @param grad_X      Gradients of the mole fraction. Flat vector with the
     *                    m_nsp in the inner loop. length = ldx * ndim
     * @param ldf         Leading dimension of the fluxes array (usually equal
     *                        to m_nsp but not always)
     * @param grad_Phi   Gradients of the electrostatic potential (length =
     *                        ndim)
     * @param Vdiff      Output of the species diffusion velocities. Flat vector
     *                   with the m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesVdiffES(size_t ndim, const doublereal* grad_T,
                                   int ldx, const doublereal* grad_X,
                                   int ldf, const doublereal* grad_Phi,
                                   doublereal* Vdiff);

    //! Return the species diffusive mass fluxes wrt to the averaged velocity in
    //! [kmol/m^2/s].
    /*!
     * The diffusive mass flux of species \e k [kmol/m^2/s] is computed
     * using the Stefan-Maxwell equation
     *
     * \f[
     *     X_i \nabla \mu_i  = RT \sum_i \frac{X_i X_j}{D_{ij}}
     *                           ( \vec{V}_j - \vec{V}_i )
     * \f]
     *
     * to determine the diffusion velocity and
     *
     * \f[
     *      \vec{N}_i = C_T X_i \vec{V}_i
     * \f]
     *
     * to determine the diffusion flux.  Here \f$ C_T \f$ is the total
     * concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$ are the Stefa-
     * Maxwell interaction parameters in [m^2/s], \f$ \vec{V}_{i} \f$ is the
     * diffusion velocity of species \e i, \f$ \mu_i \f$ is the electrochemical
     * potential of species \e i.
     *
     * Note that for this method, there is no argument for the gradient of the
     * electric potential (voltage).  Electric potential gradients can be set
     * with set_Grad_V() or method getSpeciesFluxesES() can be called.x
     *
     * The diffusion velocity is relative to an average velocity that can be
     * computed on a mole-weighted or mass-weighted basis, or the diffusion
     * velocities may be specified as relative to a specific species (i.e. a
     * solvent) all according to the `velocityBasis` input parameter.
     *
     * @param ndim       The number of spatial dimensions (1, 2, or 3).
     * @param grad_T     The temperature gradient (ignored in this model).
     *                      (length = ndim)
     * @param ldx        Leading dimension of the grad_X array.
     *                       (usually equal to m_nsp but not always)
     * @param grad_X     Gradients of the mole fraction. Flat vector with the
     *                   m_nsp in the inner loop. length = ldx * ndim
     * @param ldf        Leading dimension of the fluxes array (usually equal to
     *                        m_nsp but not always)
     * @param fluxes     Output of the diffusive mass fluxes. Flat vector with
     *                   the m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                  size_t ldx, const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

    //! Return the species diffusive mass fluxes wrt to the averaged velocity in
    //! [kmol/m^2/s].
    /*!
     * The diffusive mass flux of species \e k is computed using the Stefan-
     * Maxwell equation
     * \f[
     *     X_i \nabla \mu_i
     *                     = RT \sum_i \frac{X_i X_j}{D_{ij}}
     *                           ( \vec{V}_j - \vec{V}_i )
     * \f]
     * to determine the diffusion velocity and
     * \f[
     *      \vec{N}_i = C_T X_i \vec{V}_i
     * \f]
     * to determine the diffusion flux.  Here \f$ C_T \f$ is the total
     * concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$ are the Stefa-
     * Maxwell interaction parameters in [m^2/s], \f$ \vec{V}_{i} \f$ is the
     * diffusion velocity of species \e i, \f$ \mu_i \f$ is the electrochemical
     * potential of species \e i.
     *
     * The diffusion velocity is relative to an average velocity that can be
     * computed on a mole-weighted or mass-weighted basis, or the diffusion
     * velocities may be specified as relative to a specific species (i.e. a
     * solvent) all according to the `velocityBasis` input parameter.
     *
     * @param ndim      The number of spatial dimensions (1, 2, or 3).
     * @param grad_T    The temperature gradient (ignored in this model).
     *                     (length = ndim)
     * @param ldx       Leading dimension of the grad_X array.
     *                     (usually equal to m_nsp but not always)
     * @param grad_X    Gradients of the mole fraction. Flat vector with the
     *                  m_nsp in the inner loop. length = ldx * ndim
     * @param ldf       Leading dimension of the fluxes array (usually equal to
     *                     m_nsp but not always)
     * @param grad_Phi  Gradients of the electrostatic potential. length = ndim
     * @param fluxes    Output of the diffusive mass fluxes. Flat vector with
     *                  the m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesFluxesES(size_t ndim,
                                    const doublereal* grad_T,
                                    size_t ldx,
                                    const doublereal* grad_X,
                                    size_t ldf,
                                    const doublereal* grad_Phi,
                                    doublereal* fluxes);

    //! Return the species diffusive velocities relative to the averaged
    //! velocity.
    /*!
     * This method acts similarly to getSpeciesVdiffES() but requires all
     * gradients to be preset using methods set_Grad_X(), set_Grad_V(),
     * set_Grad_T(). See the documentation of getSpeciesVdiffES() for details.
     *
     *  @param ldf  Leading dimension of the Vdiff array.
     *  @param Vdiff  Output of the diffusive velocities. Flat vector with the
     *             m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesVdiffExt(size_t ldf, doublereal* Vdiff);

    //! Return the species diffusive fluxes relative to the averaged velocity.
    /*!
     * This method acts similarly to getSpeciesFluxesES() but requires all
     * gradients to be preset using methods set_Grad_X(), set_Grad_V(),
     * set_Grad_T(). See the documentation of getSpeciesFluxesES() for details.
     *
     * units = kg/m2/s
     *
     * @param ldf  Leading dimension of the Vdiff array.
     * @param fluxes  Output of the diffusive fluxes. Flat vector with the m_nsp
     *             in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesFluxesExt(size_t ldf, doublereal* fluxes);

protected:
    //! Returns true if temperature has changed, in which case flags are set to
    //! recompute transport properties.
    /*!
     * This is called whenever a transport property is requested. The first task
     * is to check whether the temperature has changed since the last call to
     * update_T(). If it hasn't then an immediate return is carried out.
     *
     * Note this should be a lightweight function since it's part of all of the
     * interfaces.
     *
     * @returns true if the temperature has changed, and false otherwise
     */
    virtual bool update_T();

    //! Returns true if mixture composition has changed,
    //! in which case flags are set to recompute transport properties.
    /*!
     * This is called for every interface call to check whether the
     * concentrations have changed. Concentrations change whenever the pressure
     * or the mole fraction has changed. If it has changed, the recalculations
     * should be done.
     *
     * Note this should be a lightweight function since it's part of all of the
     * interfaces.
     *
     * @returns true if the mixture composition has changed, and false
     *     otherwise.
     */
    virtual bool update_C();

    //! Updates the internal value of the gradient of the logarithm of the
    //! activity, which is used in the gradient of the chemical potential.
    /*!
     * Evaluate the gradients of the activity as they alter the diffusion
     * coefficient.
     *
     * The gradient of the chemical potential can be written in terms of
     * gradient of the logarithm of the mole fraction times a correction
     * associated with the gradient of the activity coefficient relative to that
     * of the mole fraction.  Specifically, the gradients of the logarithms of
     * each are involved according to the formula
     *
     * \f[
     *     \nabla \mu_k = RT \left[ \nabla ( \ln X_k ) +
     *     \nabla ( \ln \gamma_k ) \right] = RT \left[
     *     \nabla ( \ln a_k ) \right]
     * \f]
     *
     * The gradient in the activity coefficient requires the use of ThermoPhase
     * getdlnActCoeff that calculates its change based on a change in the state
     * (i.e. temperature and composition of each species) which was first
     * implemented in MargulesVPSSTP.cpp (LiquidTransport.h doxygen)
     */
    virtual void update_Grad_lnAC();

    //! Solve the Stefan-Maxwell equations for the diffusive fluxes.
    /*!
     * The diffusive mass flux of species \e k is computed
     * using the Stefan-Maxwell equation
     * \f[
     *     X_i \nabla \mu_i
     *                     = RT \sum_i \frac{X_i X_j}{D_{ij}}
     *                           ( \vec{V}_j - \vec{V}_i )
     * \f]
     * to determine the diffusion velocity and
     * \f[
     *      \vec{N}_i = C_T X_i \vec{V}_i
     * \f]
     * to determine the diffusion flux.  Here \f$ C_T \f$ is the total
     * concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$ are the Stefa-
     * Maxwell interaction parameters in [m^2/s], \f$ \vec{V}_{i} \f$ is the
     * diffusion velocity of species \e i, \f$ \mu_i \f$ is the electrochemical
     * potential of species \e i.
     *
     * The diffusion velocity is relative to an average velocity that can be
     * computed on a mole-weighted or mass-weighted basis, or the diffusion
     * velocities may be specified as relative to a specific species (i.e. a
     * solvent) all according to the `velocityBasis` input parameter.
     *
     * The gradient in the activity coefficient requires the use of ThermoPhase
     * getdlnActCoeff that calculates its change based on a change in the state
     * i.e. temperature and composition of each species. First implemented in
     * MargulesVPSSTP.cpp.
     *
     * One of the Stefan Maxwell equations is replaced by the appropriate
     * definition of the mass-averaged velocity, the mole-averaged velocity or
     * the specification that velocities are relative to that of one species.
     */
    void stefan_maxwell_solve();

    //! Updates the array of pure species viscosities internally.
    /*!
     * The flag m_visc_ok is set to true.
     *
     * Note that for viscosity, a positive activation energy corresponds to the
     * typical case of a positive argument to the exponential so that the
     * Arrhenius expression is
     *
     * \f[
     *      \mu = A T^n \exp( + E / R T )
     * \f]
     */
    void updateViscosity_T();

    //! Update the temperature-dependent ionic conductivity terms for each
    //! species internally
    /*!
     * The flag m_ionCond_temp_ok is set to true.
     */
    void updateIonConductivity_T();

    //! Updates the array of pure species mobility ratios internally.
    /*!
     * The flag m_mobRat_ok is set to true.
     */
    void updateMobilityRatio_T();

    //! Updates the array of pure species self diffusion coeffs internally.
    /*!
     * The flag m_selfDiff_ok is set to true.
     */
    void updateSelfDiffusion_T();

    //! Update the temperature-dependent hydrodynamic radius terms for each
    //! species internally
    /*!
     * The flag m_radi_temp_ok is set to true.
     */
    void updateHydrodynamicRadius_T();

    //! Update the temperature-dependent parts of the mixture-averaged
    //! thermal conductivity internally
    void updateCond_T();

    //! Update the concentration parts of the viscosities
    /*!
     * Internal routine is run whenever the update_boolean m_visc_conc_ok is
     * false. Currently there is no concentration dependence for the pure
     * species viscosities.
     */
    void updateViscosities_C();

    //! Update the concentration parts of the ionic conductivity
    /*!
     * Internal routine is run whenever the update_boolean m_ionCond_conc_ok is
     * false. Currently there is no concentration dependence for the pure
     * species ionic conductivity.
     */
    void updateIonConductivity_C();

    //! Update the concentration parts of the mobility ratio
    /*!
     * Internal routine is run whenever the update_boolean m_mobRat_conc_ok is
     * false. Currently there is no concentration dependence for the pure
     * species mobility ratio.
     */
    void updateMobilityRatio_C();

    //! Update the concentration parts of the self diffusion
    /*!
     * Internal routine is run whenever the update_boolean m_selfDiff_conc_ok is
     * false. Currently there is no concentration dependence for the pure
     * species self diffusion.
     */
    void updateSelfDiffusion_C();

    //! Update the concentration dependence of the hydrodynamic radius
    /*!
     * Internal routine is run whenever the update_boolean m_radi_conc_ok is
     * false. Currently there is no concentration dependence for the
     * hydrodynamic radius.
     */
    void updateHydrodynamicRadius_C();

    //! Update the binary Stefan-Maxwell diffusion coefficients
    //! wrt T using calls to the appropriate LTPspecies subclass
    void updateDiff_T();

private:
    //! Number of species squared
    size_t m_nsp2;

    //! Local copy of the molecular weights of the species
    /*!
     *  Length is equal to the number of species in the phase
     */
    vector_fp m_mw;

    //! Viscosity for each species expressed as an appropriate subclass
    //! of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     */
    std::vector<LTPspecies*> m_viscTempDep_Ns;

    //! Viscosity of the mixture expressed as a subclass of
    //! LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    LiquidTranInteraction* m_viscMixModel;

    //! Ionic conductivity for each species expressed as an appropriate subclass
    //! of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     */
    std::vector<LTPspecies*> m_ionCondTempDep_Ns;

    //! Ionic Conductivity of the mixture expressed as a subclass of
    //! LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    LiquidTranInteraction* m_ionCondMixModel;

    //! Type def for LTPvector equating it with a vector of pointers to LTPspecies
    typedef std::vector<LTPspecies*> LTPvector;

    //! Mobility ratio for the binary combinations of each species in each
    //! pure phase expressed as an appropriate subclass of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     */
    std::vector<LTPvector> m_mobRatTempDep_Ns;

    //! Mobility ratio for each binary combination of mobile species in the mixture
    //! expressed as a subclass of LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    std::vector<LiquidTranInteraction*> m_mobRatMixModel;

    //! Self Diffusion for each species in each pure species phase
    //! expressed as an appropriate subclass of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     */
    std::vector<LTPvector> m_selfDiffTempDep_Ns;

    //! Self Diffusion for each species in the mixture expressed as a subclass of
    //! LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    std::vector<LiquidTranInteraction*> m_selfDiffMixModel;

    //! Thermal conductivity for each species expressed as an
    //! appropriate subclass of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     */
    std::vector<LTPspecies*> m_lambdaTempDep_Ns;

    //! Thermal conductivity of the mixture expressed as a subclass of
    //! LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    LiquidTranInteraction* m_lambdaMixModel;

    //! (NOT USED IN LiquidTransport.)
    //! Diffusion coefficient model for each species expressed as an
    //! appropriate subclass of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData().
     *
     * Since the LiquidTransport class uses the Stefan-Maxwell equation to
     * describe species diffusivity, the species-specific diffusivity is
     * irrelevant.
     */
    std::vector<LTPspecies*> m_diffTempDep_Ns;

    //! Species diffusivity of the mixture expressed as a subclass of
    //! LiquidTranInteraction.  This will return an array of Stefan-Maxwell
    //! interaction parameters for use in the Stefan-Maxwell solution.
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    LiquidTranInteraction* m_diffMixModel;

    //! Stefan-Maxwell diffusion coefficients
    DenseMatrix m_diff_Dij;

    //! Hydrodynamic radius for each species expressed as an appropriate
    //! subclass of LTPspecies
    /*!
     * These subclasses of LTPspecies evaluate the species-specific transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidSpeciesTransportData(). length = nsp
     */
    std::vector<LTPspecies*> m_radiusTempDep_Ns;

    //! (Not used in LiquidTransport) Hydrodynamic radius of the mixture
    //! expressed as a subclass of LiquidTranInteraction
    /*!
     * These subclasses of LiquidTranInteraction evaluate the mixture transport
     * properties according to the parameters parsed in
     * TransportFactory::getLiquidInteractionsTransportData().
     */
    LiquidTranInteraction* m_radiusMixModel;

    //! Species hydrodynamic radius
    vector_fp m_hydrodynamic_radius;

    //! Internal value of the gradient of the mole fraction vector
    /*!
     * Note, this is the only gradient value that can and perhaps should reflect
     * the true state of the mole fractions in the application solution vector.
     * In other words no cropping or massaging of the values to make sure they
     * are above zero should occur. - developing ....
     *
     *  m_nsp is the number of species in the fluid
     *  k is the species index
     *  n is the dimensional index (x, y, or z). It has a length
     *    equal to m_nDim
     *
     *    m_Grad_X[n*m_nsp + k]
     */
    vector_fp m_Grad_X;

    //! Gradient of the logarithm of the activity
    /*!
     * This quantity appears in the gradient of the chemical potential. It
     * replaces the gradient of the mole fraction, and in this way serves to
     * "modify" the diffusion coefficient.
     *
     * \f[
     *     m\_Grad\_lnAC[k] = \nabla ( \ln X_k ) +
     *     \nabla ( \ln \gamma_k )
     * \f]
     *
     *  k is the species index
     *  n is the dimensional index (x, y, or z). It has a length
     *    equal to m_nDim
     *
     *    m_Grad_X[n*m_nsp + k]
     */
    vector_fp m_Grad_lnAC;

    //! Internal value of the gradient of the Temperature vector
    /*!
     * Generally, if a transport property needs this in its evaluation it will
     * look to this place to get it.
     *
     * No internal property is precalculated based on gradients. Gradients are
     * assumed to be freshly updated before every property call.
     */
    vector_fp m_Grad_T;

    //! Internal value of the gradient of the Pressure vector
    /*!
     * Generally, if a transport property needs this in its evaluation it will
     * look to this place to get it.
     *
     * No internal property is precalculated based on gradients. Gradients are
     * assumed to be freshly updated before every property call.
     */
    vector_fp m_Grad_P;

    //! Internal value of the gradient of the Electric Voltage
    /*!
     * Generally, if a transport property needs this in its evaluation it will
     * look to this place to get it.
     *
     * No internal property is precalculated based on gradients. Gradients are
     * assumed to be freshly updated before every property call.
     */
    vector_fp m_Grad_V;

    //! Gradient of the electrochemical potential
    /*!
     * m_nsp is the number of species in the fluid. k is the species index. n is
     * the dimensional index (x, y, or z)
     *
     * \f[
     *    m\_Grad\_mu[n*m_nsp + k]
     * \f]
     */
    vector_fp m_Grad_mu;

    // property values

    //! Array of Binary Diffusivities
    /*!
     * These are evaluated according to the subclass of LiquidTranInteraction
     * stored in m_diffMixModel.
     *
     * This has a size equal to nsp x nsp. It is a symmetric matrix. D_ii is the
     * self diffusion coefficient. D_ii is not needed except for when there is
     * one species in the mixture.
     *
     * units m2/sec
     */
    DenseMatrix m_bdiff;

    //! Internal value of the species viscosities
    /*!
     * Viscosity of the species evaluated using subclass of LTPspecies
     * held in m_viscTempDep_Ns.
     *
     * Length = number of species
     *
     * controlling update boolean -> m_visc_temp_ok
     */
    vector_fp m_viscSpecies;

    //! Internal value of the species ionic conductivities
    /*!
     * Ionic conductivity of the species evaluated using subclass of LTPspecies
     * held in m_ionCondTempDep_Ns.
     *
     * Length = number of species
     *
     * controlling update boolean -> m_ionCond_temp_ok
     */
    vector_fp m_ionCondSpecies;

    //! Internal value of the species mobility ratios
    /*!
     * Mobility ratio of the species evaluated using subclass of LTPspecies
     * held in m_mobRatTempDep_Ns.
     *
     * Length = number of species
     *
     * controlling update boolean -> m_mobRat_temp_ok
     */
    DenseMatrix m_mobRatSpecies;

    //! Internal value of the species self diffusion coefficients
    /*!
     * Self diffusion of the species evaluated using subclass of LTPspecies
     * held in m_selfDiffTempDep_Ns.
     *
     * Length = number of species
     *
     * controlling update boolean -> m_selfDiff_temp_ok
     */
    DenseMatrix m_selfDiffSpecies;

    //! Internal value of the species individual thermal conductivities
    /*!
     * Thermal conductivities of the species evaluated using subclass
     * of LTPspecies held in m_lambdaTempDep_Ns.
     *
     * Length = number of species
     */
    vector_fp m_lambdaSpecies;

    //! State of the mole fraction vector.
    int m_iStateMF;

    //! Local copy of the mass fractions of the species in the phase
    /*!
     * The mass fraction vector comes from the ThermoPhase object.
     *
     * length = m_nsp
     */
    vector_fp m_massfracs;

    //! Local copy of the mass fractions of the species in the phase
    /**
     * This version of the mass fraction vector is adjusted to a
     * minimum lower bound of *Tiny* for use in transport calculations.
     */
    vector_fp m_massfracs_tran;

    //! Local copy of the mole fractions of the species in the phase
    /*!
     * The mole fractions here are assumed to be bounded by 0.0 and 1.0 and they
     * are assumed to add up to one exactly. This mole fraction vector comes
     * from the ThermoPhase object. Derivative quantities from this are referred
     * to as bounded.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_molefracs;

    //! Non-zero mole fraction vector used in transport property calculations
    /*!
     * The mole fractions here are assumed to be bounded by *Tiny* and 1.0 and
     * they may not be assumed to add up to one. This mole fraction vector is
     * created from the ThermoPhase object. Derivative quantities of this use
     * the _tran suffix.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_molefracs_tran;

    //! Local copy of the concentrations of the species in the phase
    /*!
     * The concentrations are consistent with the m_molefracs vector which is
     * bounded and sums to one.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_concentrations;

    //! Local copy of the total concentration. This is consistent with the
    //! m_concentrations[] and m_molefracs[] vector.
    doublereal concTot_;

    //! Local copy of the total concentration. This is consistent with the
    //! x_molefracs_tran vector and with the concTot_ number
    doublereal concTot_tran_;

    //! Mean molecular mass
    doublereal meanMolecularWeight_;

    //! Density
    doublereal dens_;

    //! Local copy of the charge of each species. Contains the charge of each
    //! species (length m_nsp)
    vector_fp m_chargeSpecies;

    //! Specific volume for each species.  Local copy from thermo object.
    vector_fp m_volume_spec;

    //! Vector of activity coefficients
    vector_fp m_actCoeff;

    //! RHS to the Stefan-Maxwell equation
    DenseMatrix m_B;

    //! Matrix for the Stefan-Maxwell equation.
    DenseMatrix m_A;

    //! Current Temperature -> locally stored. This is used to test whether new
    //! temperature computations should be performed.
    doublereal m_temp;

    //! Current value of the pressure
    doublereal m_press;

    //! Solution of the Stefan Maxwell equation in terms of flux. This is the
    //! mass flux of species k in units of kg m-3 s-1.
    Array2D m_flux;

    //! Solution of the Stefan Maxwell equation. This is the diffusion velocity
    //! of species k in units of m/s and relative to the mole-averaged velocity.
    Array2D m_Vdiff;

    //! Saved value of the mixture thermal conductivity
    doublereal m_lambda;

    //! Saved value of the mixture viscosity
    doublereal m_viscmix;

    //! Saved value of the mixture ionic conductivity
    doublereal m_ionCondmix;

    //! Saved values of the mixture mobility ratios
    vector_fp m_mobRatMix;

    //! Saved values of the mixture self diffusion coefficients
    vector_fp m_selfDiffMix;

    //! work space. Length is equal to m_nsp
    vector_fp m_spwork;

private:
    //! Boolean indicating that the top-level mixture viscosity is current. This
    //! is turned false for every change in T, P, or C.
    bool m_visc_mix_ok;

    //! Boolean indicating that weight factors wrt viscosity is current
    bool m_visc_temp_ok;

    //! Flag to indicate that the pure species viscosities are current wrt the
    //! concentration
    bool m_visc_conc_ok;

    //! Boolean indicating that the top-level mixture ionic conductivity is
    //! current. This is turned false for every change in T, P, or C.
    bool m_ionCond_mix_ok;

    //! Boolean indicating that weight factors wrt ionic conductivity is current
    bool m_ionCond_temp_ok;

    //! Flag to indicate that the pure species ionic conductivities
    //! are current wrt the concentration
    bool m_ionCond_conc_ok;

    //! Flag to indicate that the mixture conductivity is current
    bool m_cond_mix_ok;

    //! Boolean indicating that the top-level mixture mobility ratio is current.
    //! This is turned false for every change in T, P, or C.
    bool m_mobRat_mix_ok;

    //! Boolean indicating that weight factors wrt mobility ratio is current
    bool m_mobRat_temp_ok;

    //! Flag to indicate that the pure species mobility ratios
    //! are current wrt the concentration
    bool m_mobRat_conc_ok;

    //! Boolean indicating that the top-level mixture self diffusion is current.
    //! This is turned false for every change in T, P, or C.
    bool m_selfDiff_mix_ok;

    //! Boolean indicating that weight factors wrt self diffusion is current
    bool m_selfDiff_temp_ok;

    //! Flag to indicate that the pure species self diffusion are current wrt
    //! the concentration
    bool m_selfDiff_conc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_radi_mix_ok;

    //! Boolean indicating that temperature dependence of
    //! hydrodynamic radius is current
    bool m_radi_temp_ok;

    //! Flag to indicate that the hydrodynamic radius is current is current wrt
    //! the concentration
    bool m_radi_conc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_diff_mix_ok;

    //! Boolean indicating that binary diffusion coeffs are current
    bool m_diff_temp_ok;

    //! Flag to indicate that the pure species conductivities are current wrt
    //! the temperature
    bool m_lambda_temp_ok;

    //! Boolean indicating that mixture conductivity is current
    bool m_lambda_mix_ok;

    //! Mode indicator for transport models -- currently unused.
    int m_mode;

    //! Debugging flags. Turn on to get debugging information.
    bool m_debug;
};
}
#endif
