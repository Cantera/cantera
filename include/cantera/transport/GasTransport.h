/**
 * @file GasTransport.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GAS_TRANSPORT_H
#define CT_GAS_TRANSPORT_H

#include "Transport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{

class MMCollisionInt;

//! Class GasTransport implements some functions and properties that are
//! shared by the MixTransport and MultiTransport classes.
//!
//! For details, see Kee, et al. @cite kee2003 and @cite kee2017.
//!
//! @ingroup tranprops
class GasTransport : public Transport
{
public:
    //! Viscosity of the mixture  (kg /m /s)
    /*!
     * The viscosity is computed using the Wilke mixture rule (kg /m /s)
     *
     * @f[
     *     \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
     * @f]
     *
     * Here @f$ \mu_k @f$ is the viscosity of pure species @e k, and
     *
     * @f[
     *     \Phi_{k,j} = \frac{\left[1
     *                  + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
     *                  {\sqrt{8}\sqrt{1 + M_k/M_j}}
     * @f]
     *
     * @returns the viscosity of the mixture (units =  Pa s = kg /m /s)
     *
     * @see updateViscosity_T()
     */
    double viscosity() override;

    //! Get the pure-species viscosities
    void getSpeciesViscosities(double* const visc) override {
        update_T();
        updateViscosity_T();
        std::copy(m_visc.begin(), m_visc.end(), visc);
    }

    //! Returns the matrix of binary diffusion coefficients.
    /*!
     * d[ld*j + i] = rp * m_bdiff(i,j);
     *
     * @param ld   offset of rows in the storage
     * @param d    output vector of diffusion coefficients. Units of m**2 / s
     */
    void getBinaryDiffCoeffs(const size_t ld, double* const d) override;

    //! Returns the Mixture-averaged diffusion coefficients [m^2/s].
    /*!
     * Returns the mixture averaged diffusion coefficients for a gas,
     * appropriate for calculating the mass averaged diffusive flux with respect
     * to the mass averaged velocity using gradients of the mole fraction.
     * Note, for the single species case or the pure fluid case the routine
     * returns the self-diffusion coefficient. This is needed to avoid a Nan
     * result in the formula below.
     *
     * This is Eqn. 12.180 from "Chemically Reacting Flow"
     *
     * @f[
     *     D_{km}' = \frac{\left( \bar{M} - X_k M_k \right)}{ \bar{\qquad M \qquad } }  {\left( \sum_{j \ne k} \frac{X_j}{D_{kj}} \right) }^{-1}
     * @f]
     *
     * @param[out] d  Vector of mixture diffusion coefficients, @f$ D_{km}' @f$ ,
     *     for each species (m^2/s). length m_nsp
     */
    void getMixDiffCoeffs(double* const d) override;

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    //! These are the coefficients for calculating the molar diffusive fluxes
    //! from the species mole fraction gradients, computed according to
    //! Eq. 12.176 in "Chemically Reacting Flow":
    //!
    //! @f[  D_{km}^* = \frac{1-X_k}{\sum_{j \ne k}^K X_j/\mathcal{D}_{kj}} @f]
    //!
    //! @param[out] d vector of mixture-averaged diffusion coefficients for
    //!     each species, length m_nsp.
    void getMixDiffCoeffsMole(double* const d) override;

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    /*!
     * These are the coefficients for calculating the diffusive mass fluxes
     * from the species mass fraction gradients, computed according to
     * Eq. 12.178 in "Chemically Reacting Flow":
     *
     * @f[
     *     \frac{1}{D_{km}} = \sum_{j \ne k}^K \frac{X_j}{\mathcal{D}_{kj}} +
     *     \frac{X_k}{1-Y_k} \sum_{j \ne k}^K \frac{Y_j}{\mathcal{D}_{kj}}
     * @f]
     *
     * @param[out] d vector of mixture-averaged diffusion coefficients for
     *     each species, length m_nsp.
     */
    void getMixDiffCoeffsMass(double* const d) override;

    //! Return the polynomial fits to the viscosity of species i
    //! @see fitProperties()
    void getViscosityPolynomial(size_t i, double* coeffs) const override;

    //! Return the temperature fits of the heat conductivity of species i
    //! @see fitProperties()
    void getConductivityPolynomial(size_t i, double* coeffs) const override;

    //! Return the polynomial fits to the binary diffusivity of species pair (i, j)
    //! @see fitDiffCoeffs()
    void getBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs) const override;

    //! Return the polynomial fits to the collision integral of species pair (i, j)
    //! @see fitCollisionIntegrals()
    void getCollisionIntegralPolynomial(size_t i, size_t j,
                                        double* astar_coeffs,
                                        double* bstar_coeffs,
                                        double* cstar_coeffs) const override;

    //! Modify the polynomial fits to the viscosity of species i
    //! @see fitProperties()
    void setViscosityPolynomial(size_t i, double* coeffs) override;

    //! Modify the temperature fits of the heat conductivity of species i
    //! @see fitProperties()
    void setConductivityPolynomial(size_t i, double* coeffs) override;

    //! Modify the polynomial fits to the binary diffusivity of species pair (i, j)
    //! @see fitDiffCoeffs()
    void setBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs) override;

    //! Modify the polynomial fits to the collision integral of species pair (i, j)
    //! @see fitCollisionIntegrals()
    void setCollisionIntegralPolynomial(size_t i, size_t j,
                                        double* astar_coeffs,
                                        double* bstar_coeffs,
                                        double* cstar_coeffs, bool actualT) override;

    void init(ThermoPhase* thermo, int mode=0, int log_level=-7) override;

    bool CKMode() const override {
        return m_mode == CK_Mode;
    }

    void invalidateCache() override;

protected:
    GasTransport();

    virtual void update_T();
    virtual void update_C() = 0;

    //! Update the temperature-dependent viscosity terms.
    /**
     * Updates the array of pure species viscosities, and the weighting
     * functions in the viscosity mixture rule. The flag m_visc_ok is set to true.
     *
     * The formula for the weighting function is from Poling et al. @cite poling2001,
     * Eq. (9-5.14):
     *  @f[
     *      \phi_{ij} = \frac{ \left[ 1 + \left( \mu_i / \mu_j \right)^{1/2} \left( M_j / M_i \right)^{1/4} \right]^2 }
     *                    {\left[ 8 \left( 1 + M_i / M_j \right) \right]^{1/2}}
     *  @f]
     */
    virtual void updateViscosity_T();

    //! Update the pure-species viscosities. These are evaluated from the
    //! polynomial fits of the temperature and are assumed to be independent
    //! of pressure.
    virtual void updateSpeciesViscosities();

    //! Update the binary diffusion coefficients
    /*!
     * These are evaluated from the polynomial fits of the temperature at the
     * unit pressure of 1 Pa.
     */
    virtual void updateDiff_T();

    //! @name Initialization
    //! @{

    //! Setup parameters for a new kinetic-theory-based transport manager for
    //! low-density gases
    virtual void setupCollisionParameters();

    //! Setup range for polynomial fits to collision integrals of
    //! Monchick & Mason @cite monchick1961
    void setupCollisionIntegral();

    //! Read the transport database
    /*!
     * Read transport property data from a file for a list of species. Given the
     * name of a file containing transport property parameters and a list of
     * species names.
     */
    void getTransportData();

    //! Corrections for polar-nonpolar binary diffusion coefficients
    /*!
     * Calculate corrections to the well depth parameter and the diameter for use in
     * computing the binary diffusion coefficient of polar-nonpolar pairs. For more
     * information about this correction, see Dixon-Lewis @cite dixon-lewis1968.
     *
     * @param i        Species one - this is a bimolecular correction routine
     * @param j        species two - this is a bimolecular correction routine
     * @param f_eps    Multiplicative correction factor to be applied to epsilon(i,j)
     * @param f_sigma  Multiplicative correction factor to be applied to diam(i,j)
     */
    void makePolarCorrections(size_t i, size_t j, double& f_eps,
                              double& f_sigma);

    //! Generate polynomial fits to collision integrals
    /*!
     * @param integrals interpolator for the collision integrals
     */
    void fitCollisionIntegrals(MMCollisionInt& integrals);

    //! Generate polynomial fits to the viscosity @f$ \eta @f$ and conductivity
    //! @f$ \lambda @f$.
    /*!
     * If CK_mode, then the fits are of the form
     * @f[
     *      \ln \eta(i) = \sum_{n=0}^3 a_n(i) \, (\ln T)^n
     * @f]
     * and
     * @f[
     *      \ln \lambda(i) = \sum_{n=0}^3 b_n(i) \, (\ln T)^n
     * @f]
     * Otherwise the fits are of the form
     * @f[
     *      \left(\eta(i)\right)^{1/2} = T^{1/4} \sum_{n=0}^4 a_n(i) \, (\ln T)^n
     * @f]
     * and
     * @f[
     *      \lambda(i) = T^{1/2} \sum_{n=0}^4 b_n(i) \, (\ln T)^n
     * @f]
     *
     * @param integrals interpolator for the collision integrals
     */
    virtual void fitProperties(MMCollisionInt& integrals);

    //! Generate polynomial fits to the binary diffusion coefficients
    /*!
     * If CK_mode, then the fits are of the form
     * @f[
     *      \ln D(i,j) = \sum_{n=0}^3 c_n(i,j) \, (\ln T)^n
     * @f]
     * Otherwise the fits are of the form
     * @f[
     *      D(i,j) = T^{3/2} \sum_{n=0}^4 c_n(i,j) \, (\ln T)^n
     * @f]
     *
     * @param integrals interpolator for the collision integrals
     */
    virtual void fitDiffCoeffs(MMCollisionInt& integrals);

    //! Second-order correction to the binary diffusion coefficients
    /*!
     * Calculate second-order corrections to binary diffusion coefficient pair
     * (dkj, djk). At first order, the binary diffusion coefficients are
     * independent of composition, and d(k,j) = d(j,k). But at second order,
     * there is a weak dependence on composition, with the result that d(k,j) !=
     * d(j,k). This method computes the multiplier by which the first-order
     * binary diffusion coefficient should be multiplied to produce the value
     * correct to second order. The expressions here are taken from Marerro and
     * Mason @cite marrero1972.
     *
     * @param t   Temperature (K)
     * @param integrals interpolator for the collision integrals
     * @param k   index of first species
     * @param j   index of second species
     * @param xk  Mole fraction of species k
     * @param xj  Mole fraction of species j
     * @param fkj multiplier for d(k,j)
     * @param fjk multiplier for d(j,k)
     *
     * @note This method is not used currently.
     */
    void getBinDiffCorrection(double t, MMCollisionInt& integrals, size_t k,
                              size_t j, double xk, double xj,
                              double& fkj, double& fjk);

    //! @}

    //! Vector of species mole fractions. These are processed so that all mole
    //! fractions are >= *Tiny*. Length = m_kk.
    vector<double> m_molefracs;

    //! Internal storage for the viscosity of the mixture  (kg /m /s)
    double m_viscmix = 0.0;

    //! Update boolean for mixture rule for the mixture viscosity
    bool m_visc_ok = false;

    //! Update boolean for the weighting factors for the mixture viscosity
    bool m_viscwt_ok = false;

    //! Update boolean for the species viscosities
    bool m_spvisc_ok = false;

    //! Update boolean for the binary diffusivities at unit pressure
    bool m_bindiff_ok = false;

    //! Type of the polynomial fits to temperature. `CK_Mode` means Chemkin mode.
    //! Any other value means to use %Cantera's preferred fitting functions.
    int m_mode = 0;

    //! m_phi is a Viscosity Weighting Function. size = m_nsp * n_nsp
    DenseMatrix m_phi;

    //! work space length = m_kk
    vector<double> m_spwork;

    //! vector of species viscosities (kg /m /s). These are used in Wilke's
    //! rule to calculate the viscosity of the solution. length = m_kk.
    vector<double> m_visc;

    //! Polynomial fits to the viscosity of each species. m_visccoeffs[k] is
    //! the vector of polynomial coefficients for species k that fits the
    //! viscosity as a function of temperature.
    vector<vector<double>> m_visccoeffs;

    //! Local copy of the species molecular weights.
    vector<double> m_mw;

    //! Holds square roots of molecular weight ratios
    /*!
     *  @code
     *  m_wratjk(j,k)  = sqrt(mw[j]/mw[k])        j < k
     *  m_wratjk(k,j)  = sqrt(sqrt(mw[j]/mw[k]))  j < k
     *  @endcode
     */
    DenseMatrix m_wratjk;

    //! Holds square roots of molecular weight ratios
    /*!
     *  `m_wratjk1(j,k)  = sqrt(1.0 + mw[k]/mw[j])        j < k`
     */
    DenseMatrix m_wratkj1;

    //! vector of square root of species viscosities sqrt(kg /m /s). These are
    //! used in Wilke's rule to calculate the viscosity of the solution.
    //! length = m_kk.
    vector<double> m_sqvisc;

    //! Powers of the ln temperature, up to fourth order
    vector<double> m_polytempvec;

    //! Current value of the temperature at which the properties in this object
    //! are calculated (Kelvin).
    double m_temp = -1.0;

    //! Current value of Boltzmann constant times the temperature (Joules)
    double m_kbt = 0.0;

    //! current value of temperature to 1/2 power
    double m_sqrt_t = 0.0;

    //! Current value of the log of the temperature
    double m_logt = 0.0;

    //! Current value of temperature to 1/4 power
    double m_t14 = 0.0;

    //! Polynomial fits to the binary diffusivity of each species
    /*!
     * m_diffcoeff[ic] is vector of polynomial coefficients for species i
     * species j that fits the binary diffusion coefficient. The relationship
     * between i j and ic is determined from the following algorithm:
     *
     *      int ic = 0;
     *      for (i = 0; i < m_nsp; i++) {
     *         for (j = i; j < m_nsp; j++) {
     *           ic++;
     *         }
     *      }
     */
    vector<vector<double>> m_diffcoeffs;

    //! Matrix of binary diffusion coefficients at the reference pressure and
    //! the current temperature Size is nsp x nsp.
    DenseMatrix m_bdiff;

    //! temperature fits of the heat conduction
    /*!
     *  Dimensions are number of species (nsp) polynomial order of the collision
     *  integral fit (degree+1).
     */
    vector<vector<double>> m_condcoeffs;

    //! Indices for the (i,j) interaction in collision integral fits
    /*!
     *  m_poly[i][j] contains the index for (i,j) interactions in
     *  #m_omega22_poly, #m_astar_poly, #m_bstar_poly, and #m_cstar_poly.
     */
    vector<vector<int>> m_poly;

    //! Fit for omega22 collision integral
    /*!
     * m_omega22_poly[m_poly[i][j]] is the vector of polynomial coefficients
     * (length degree+1) for the collision integral fit for the species pair
     * (i,j).
     */
    vector<vector<double>> m_omega22_poly;

    //! Flag to indicate for which (i,j) interaction pairs the
    //! actual temperature is used instead of the reduced temperature
    vector<vector<int>> m_star_poly_uses_actualT;

    //! Fit for astar collision integral
    /*!
     * m_astar_poly[m_poly[i][j]] is the vector of polynomial coefficients
     * (length degree+1) for the collision integral fit for the species pair
     * (i,j).
     */
    vector<vector<double>> m_astar_poly;

    //! Fit for bstar collision integral
    /*!
     * m_bstar_poly[m_poly[i][j]] is the vector of polynomial coefficients
     * (length degree+1) for the collision integral fit for the species pair
     * (i,j).
     */
    vector<vector<double>> m_bstar_poly;

    //! Fit for cstar collision integral
    /*!
     * m_bstar_poly[m_poly[i][j]] is the vector of polynomial coefficients
     * (length degree+1) for the collision integral fit for the species pair
     * (i,j).
     */
    vector<vector<double>> m_cstar_poly;

    //! Rotational relaxation number for each species
    /*!
     * length is the number of species in the phase. units are dimensionless
     */
    vector<double> m_zrot;

    //! Dimensionless rotational heat capacity of each species
    /*!
     * These values are 0, 1 and 1.5 for single-molecule, linear, and nonlinear
     * species respectively length is the number of species in the phase.
     * Dimensionless  (Cr / R)
     */
    vector<double> m_crot;

    //! Vector of booleans indicating whether a species is a polar molecule
    /*!
     * Length is nsp
     */
    vector<bool> m_polar;

    //! Polarizability of each species in the phase
    /*!
     * Length = nsp. Units = m^3
     */
    vector<double> m_alpha;

    //! Lennard-Jones well-depth of the species in the current phase
    /*!
     * length is the number of species in the phase. Units are Joules (Note this
     * is not Joules/kmol) (note, no kmol -> this is a per molecule amount)
     */
    vector<double> m_eps;

    //! Lennard-Jones diameter of the species in the current phase
    /*!
     * length is the number of species in the phase. units are in meters.
     */
    vector<double> m_sigma;

    //! This is the reduced mass of the interaction between species i and j
    /*!
     *  reducedMass(i,j) =  mw[i] * mw[j] / (Avogadro * (mw[i] + mw[j]));
     *
     *  Units are kg (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix
     */
    DenseMatrix m_reducedMass;

    //! hard-sphere diameter for (i,j) collision
    /*!
     *  diam(i,j) = 0.5*(sigma[i] + sigma[j]);
     *  Units are m (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix m_diam;

    //! The effective well depth for (i,j) collisions
    /*!
     *     epsilon(i,j) = sqrt(eps[i]*eps[j]);
     *     Units are Joules (note, no kmol -> this is a per molecule amount)
     *
     * Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix m_epsilon;

    //! The effective dipole moment for (i,j) collisions
    /*!
     *  Given `dipoleMoment` in Debye (a Debye is 3.335e-30 C-m):
     *
     *    dipole(i,i) = 1.e-21 / lightSpeed * dipoleMoment;
     *    dipole(i,j) = sqrt(dipole(i,i) * dipole(j,j));
     *  (note, no kmol -> this is a per molecule amount)
     *
     *  Length nsp * nsp. This is a symmetric matrix.
     */
    DenseMatrix m_dipole;

    //! Reduced dipole moment of the interaction between two species
    /*!
     *  This is the reduced dipole moment of the interaction between two species
     *       0.5 * dipole(i,j)^2 / (4 * Pi * epsilon_0 * epsilon(i,j) * d^3);
     *
     *  Length nsp * nsp .This is a symmetric matrix
     */
    DenseMatrix m_delta;

    //! Pitzer acentric factor
    /*!
     * Length is the number of species in the phase. Dimensionless.
     */
    vector<double> m_w_ac;

    //! Dispersion coefficient
    vector<double> m_disp;

    //! Quadrupole polarizability
    vector<double> m_quad_polar;

    //! Level of verbose printing during initialization.
    //! @deprecated To be removed after %Cantera 3.1.
    int m_log_level = 0;
};

} // namespace Cantera

#endif
