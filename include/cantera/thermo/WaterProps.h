/**
 *  @file WaterProps.h
 *   Header for a class used to house several approximation
 *   routines for properties of water.
 *  (see \ref thermoprops
 *   and class \link Cantera::WaterProps WaterProps\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WATERPROPS_H
#define CT_WATERPROPS_H


#include "cantera/base/ct_defs.h"
namespace Cantera
{
class WaterPropsIAPWS;
class PDSS_Water;

/**
 * @defgroup relatedProps Electric Properties of Phases
 *
 * Computation of the electric properties of phases
 *
 * ### Treatment of the phase potential and the electrochemical potential of a species
 *
 * The electrochemical potential of species \f$k\f$ in a phase \f$p\f$, \f$ \zeta_k \f$,
 * is related to the chemical potential via the following equation,
 *
 *       \f[
 *            \zeta_{k}(T,P) = \mu_{k}(T,P) + z_k \phi_p
 *       \f]
 *
 * where  \f$ \nu_k \f$ is the charge of species \f$k\f$, and \f$ \phi_p \f$ is
 * the electric potential of phase \f$p\f$.
 *
 * The potential  \f$ \phi_p \f$ is tracked and internally stored within the
 * base ThermoPhase object. It constitutes a specification of the internal state
 * of the phase; it's the third state variable, the first two being temperature
 * and density (or, pressure, for incompressible equations of state). It may be
 * set with the function, ThermoPhase::setElectricPotential(), and may be
 * queried with the function ThermoPhase::electricPotential().
 *
 * Note, the overall electrochemical potential of a phase may not be changed
 * by the potential because many phases enforce charge neutrality:
 *
 *       \f[
 *            0 = \sum_k z_k X_k
 *       \f]
 *
 * Whether charge neutrality is necessary for a phase is also specified within
 * the ThermoPhase object, by the function call
 * ThermoPhase::chargeNeutralityNecessary(). Note, that it is not necessary for
 * the IdealGas phase, currently. However, it is necessary for liquid phases
 * such as DebyeHuckel and HMWSoln for the proper specification of the chemical
 * potentials.
 *
 * This equation, when applied to the \f$ \zeta_k \f$ equation described
 * above, results in a zero net change in the effective Gibbs free energy of
 * the phase. However, specific charged species in the phase may increase or
 * decrease their electrochemical potentials, which will have an effect on
 * interfacial reactions involving charged species, when there is a potential
 * drop between phases. This effect is used within the InterfaceKinetics and
 * EdgeKinetics kinetics objects classes.
 *
 * ### Electrothermochemical Properties of Phases of Matter
 *
 * The following classes are used to compute the electrical and
 * electrothermochemical properties of phases of matter. The main property
 * currently is the dielectric constant, which is an important parameter for
 * electrolyte solutions. The class WaterProps calculate the dielectric
 * constant of water as a function of temperature and pressure.
 *
 * WaterProps also calculate the constant A_debye used in the Debye Huckel and
 * Pitzer activity coefficient calculations.
 *
 * @ingroup phases
 */
//@{

//! The WaterProps class is used to house several approximation routines for
//! properties of water.
/*!
 * The class is also a wrapper around the WaterPropsIAPWS class which provides
 * the calculations for the equation of state properties for water.
 *
 * In particular, this class house routine for the calculation of the dielectric
 * constant of water
 *
 * Most if not all of the member functions are static.
 */
class WaterProps
{
public:
    //! Default constructor
    WaterProps();

    //! Constructor
    /*!
     * @param wptr Pointer to WaterPropsIAPWS object
     */
    WaterProps(WaterPropsIAPWS* wptr);

    //! Constructor with pointer to Water PDSS object
    /*!
     * @param wptr Pointer to water standard state object
     */
    WaterProps(PDSS_Water* wptr);

    // WaterProps objects are not copyable or assignable
    WaterProps(const WaterProps& b) = delete;
    WaterProps& operator=(const WaterProps& b) = delete;
    virtual ~WaterProps();

    //! Simple calculation of water density at atmospheric pressure.
    //! Valid up to boiling point.
    /*!
     * This formulation has no dependence on the pressure and shouldn't be used
     * where accuracy is needed.
     *
     * @param T temperature in kelvin
     * @param P Pressure in pascal
     * @param ifunc changes what's returned
     *
     * @return value returned depends on ifunc value:
     *   - ifunc = 0 Returns the density in kg/m^3
     *   - ifunc = 1 returns the derivative of the density wrt T.
     *   - ifunc = 2 returns the 2nd derivative of the density wrt T
     *   - ifunc = 3 returns the derivative of the density wrt P.
     *
     * Verification:
     *   Agrees with the CRC values (6-10) for up to 4 sig digits.
     *
     * units = returns density in kg m-3.
     */
    static doublereal density_T(doublereal T, doublereal P, int ifunc);

    //! Bradley-Pitzer equation for the dielectric constant
    //! of water as a function of temperature and pressure.
    /*!
     * Returns the dimensionless relative dielectric constant and its
     * derivatives.
     *
     * Range of validity: 0 to 350C, 0 to 1 kbar pressure
     *
     * @param T temperature (kelvin)
     * @param P_pascal pressure in pascal
     * @param ifunc changes what's returned from the function
     * @return Depends on the value of ifunc:
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *
     * Validation: Numerical experiments indicate that this function agrees with
     * the Archer and Wang data in the CRC p. 6-10 to all 4 significant digits
     * shown (0 to 100C).
     *
     * value at 25C and 1 atm, relEps = 78.38
     */
    doublereal relEpsilon(doublereal T, doublereal P_pascal, int ifunc = 0);

    //! ADebye calculates the value of A_Debye as a function of temperature and
    //! pressure according to relations that take into account the temperature
    //! and pressure dependence of the water density and dielectric constant.
    /*!
     * The A_Debye expression appears on the top of the ln actCoeff term in the
     * general Debye-Huckel expression It depends on temperature and pressure.
     * And, therefore, most be recalculated whenever T or P changes. The units
     * returned by this expression are sqrt(kg/gmol).
     *
     * \f[
     *   A_{Debye} = \frac{1}{8 \pi} \sqrt{\frac{2 N_{Avog} \rho_w}{1000}}
     *                     {\left(\frac{e^2}{\epsilon k_{boltz} T}\right)}^{\frac{3}{2}}
     * \f]
     *
     * Nominal value at 25C and 1atm = 1.172576 sqrt(kg/gmol).
     *
     * Based on:
     *    - epsilon/epsilon_0 = 78.54 (water at 25C)
     *    - T = 298.15 K
     *    - B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *
     * @param T  Temperature (kelvin)
     * @param P  pressure (pascal)
     * @param ifunc Changes what's returned from the routine
     * @returns a double whose meaning depends on ifunc:
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *
     * Verification: With the epsRelWater value from the Bradley-Pitzer
     * relation, and the water density from the density_IAPWS() function, The
     * A_Debye computed with this function agrees with the Pitzer table p. 99 to
     * 4 significant digits at 25C. and 20C. (Aphi = ADebye/3)
     */
    doublereal ADebye(doublereal T, doublereal P, int ifunc);

    //! Returns the saturation pressure given the temperature
    /*!
     * @param T temperature (kelvin)
     * @returns the saturation pressure (pascal)
     */
    doublereal satPressure(doublereal T);

    //! Returns the density of water
    /*!
     * This function sets the internal temperature and pressure
     * of the underlying object at the same time.
     *
     * @param T Temperature (kelvin)
     * @param P pressure (pascal)
     */
    doublereal density_IAPWS(doublereal T, doublereal P);

    //! Returns the density of water
    /*!
     * This function uses the internal state of the underlying water object
     */
    doublereal density_IAPWS() const;

    //! returns the coefficient of thermal expansion
    /*!
     * @param T Temperature (kelvin)
     * @param P pressure (pascal)
     */
    doublereal coeffThermalExp_IAPWS(doublereal T, doublereal P);

    //! Returns the isothermal compressibility of water
    /*!
     * @param T  temperature in kelvin
     * @param P  pressure in pascal
     */
    doublereal isothermalCompressibility_IAPWS(doublereal T, doublereal P);

protected:
    //! Pointer to the WaterPropsIAPWS object
    WaterPropsIAPWS* m_waterIAPWS;

    //! true if we own the WaterPropsIAPWS object
    bool m_own_sub;
};

//@}
}

#endif
