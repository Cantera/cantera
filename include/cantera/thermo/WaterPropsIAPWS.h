/**
 * @file WaterPropsIAPWS.h
 * Headers for a class for calculating the equation of state of water
 * from the IAPWS 1995 Formulation based on the steam tables thermodynamic
 * basis (See class \link Cantera::WaterPropsIAPWS WaterPropsIAPWS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef WATERPROPSIAPWS_H
#define WATERPROPSIAPWS_H

#include "WaterPropsIAPWSphi.h"

namespace Cantera
{
/**
 * @name Names for the phase regions
 *
 * These constants are defined and used in the interface to describe the
 * location of where we are in (T,rho) space.
 *
 * WATER_UNSTABLELIQUID indicates that we are in the unstable region, inside the
 * spinodal curve where dpdrho < 0.0 amongst other properties. The difference
 * between WATER_UNSTABLELIQUID and WATER_UNSTABLEGAS is that
 *     for WATER_UNSTABLELIQUID  d2pdrho2 > 0   and dpdrho < 0.0
 *     for WATER_UNSTABLEGAS     d2pdrho2 < 0   and dpdrho < 0.0
 */
//@{
#define WATER_GAS 0
#define WATER_LIQUID 1
#define WATER_SUPERCRIT 2
#define WATER_UNSTABLELIQUID 3
#define WATER_UNSTABLEGAS 4
//@}

//! Class for calculating the equation of state of water.
/*!
 * The reference is W. Wagner, A. Pruss, "The IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use," J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * This class provides a very complicated polynomial for the specific
 * Helmholtz free energy of water, as a function of temperature and density.
 *
 * \f[
 *     \frac{M\hat{f}(\rho,T)}{R T} = \phi(\delta, \tau) =
 *                     \phi^o(\delta, \tau) +  \phi^r(\delta, \tau)
 * \f]
 *
 * where
 *
 * \f[
 *     \delta = \rho / \rho_c \quad \mathrm{and} \quad \tau = T_c / T
 * \f]
 *
 * The following constants are assumed
 *
 * \f[
 *     T_c = 647.096\mathrm{\;K}
 * \f]
 * \f[
 *     \rho_c = 322 \mathrm{\;kg\,m^{-3}}
 * \f]
 * \f[
 *     R/M = 0.46151805 \mathrm{\;kJ\,kg^{-1}\,K^{-1}}
 * \f]
 *
 * The free energy is a unique single-valued function of the temperature and
 * density over its entire range.
 *
 * Note, the base thermodynamic state for this class is the one used in the
 * steam tables, i.e., the liquid at the triple point for water has the
 * following properties:
 *
 *   - u(273.16, rho)    = 0.0
 *   - s(273.16, rho)    = 0.0
 *   - psat(273.16)      = 611.655 Pascal
 *   - rho(273.16, psat) = 999.793 kg m-3
 *
 * Therefore, to use this class within %Cantera, offsets to u() and s() must
 * be used to put the water class onto the same basis as other thermodynamic
 * quantities. For example, in the WaterSSTP class, these offsets are
 * calculated in the following way. The thermodynamic base state for water is
 * set to the NIST basis here by specifying constants EW_Offset and SW_Offset.
 * These offsets are calculated on the fly so that the following properties
 * hold:
 *
 *   - Delta_Hfo_idealGas(298.15, 1bar) = -241.826 kJ/gmol
 *   - So_idealGas(298.15, 1bar)        = 188.835 J/gmolK
 *
 * The offsets are calculated by actually computing the above quantities and
 * then calculating the correction factor.
 *
 * This class provides an interface to the WaterPropsIAPWSphi class, which
 * actually calculates the \f$ \phi^o(\delta, \tau)  \f$ and the
 * \f$ \phi^r(\delta, \tau) \f$ polynomials in dimensionless form.
 *
 * All thermodynamic results from this class are returned in dimensional form.
 * This is because the gas constant (and molecular weight) used within this
 * class is allowed to be potentially different than that used elsewhere in
 * %Cantera. Therefore, everything has to be in dimensional units. Note,
 * however, the thermodynamic basis is set to that used in the steam tables.
 * (u = s = 0 for liquid water at the triple point).
 *
 * This class is not a ThermoPhase. However, it does maintain an internal
 * state of the object that is dependent on temperature and density. The
 * internal state is characterized by an internally stored \f$ \tau\f$ and a
 * \f$ \delta \f$ value, and an iState value, which indicates whether the
 * point is a liquid, a gas, or a supercritical fluid. Along with that the
 * \f$ \tau\f$ and a \f$ \delta \f$ values are polynomials of \f$ \tau\f$ and
 * a \f$ \delta \f$ that are kept by the WaterPropsIAPWSphi class. Therefore,
 * whenever  \f$ \tau\f$ or \f$ \delta \f$ is changed, the function setState()
 * must be called in order for the internal state to be kept up to date.
 *
 * The class is pretty straightforward. However, one function deserves
 * mention. The density() function calculates the density that is consistent
 * with a particular value of the temperature and pressure. It may therefore
 * be multivalued or potentially there may be no answer from this function. It
 * therefore takes a phase guess and a density guess as optional parameters.
 * If no guesses are supplied to density(), a gas phase guess is assumed. This
 * may or may not be what is wanted. Therefore, density() should usually at
 * least be supplied with a phase guess so that it may manufacture an
 * appropriate density guess. density() manufactures the initial density
 * guess, nondimensionalizes everything, and then calls
 * WaterPropsIAPWSphi::dfind(), which does the iterative calculation to find
 * the density condition that matches the desired input pressure.
 *
 * The phase guess defines are located in the .h file. they are
 *
 *   - WATER_GAS
 *   - WATER_LIQUID
 *   - WATER_SUPERCRIT
 *
 * There are only three functions which actually change the value of the
 * internal state of this object after it's been instantiated
 *
 *   - setState_TR(temperature, rho)
 *   - density(temperature, pressure, phase, rhoguess)
 *   - psat(temperature, waterState);
 *
 * The setState_TR() is the main function that sets the temperature and rho
 * value. The density() function serves as a setState_TP() function, in that
 * it sets internal state to a temperature and pressure. However, note that
 * this is potentially multivalued. Therefore, we need to supply in addition a
 * phase guess and a rho guess to the input temperature and pressure. The
 * psat() function sets the internal state to the saturated liquid or
 * saturated gas state, depending on the waterState parameter.
 *
 * Because the underlying object WaterPropsIAPWSphi is privately held, you can
 * be sure that the underlying state of this object doesn't change except due
 * to the three function calls listed above.
 *
 * @ingroup thermoprops
 */
class WaterPropsIAPWS
{
public:
    //! Base constructor
    WaterPropsIAPWS();

    WaterPropsIAPWS(const WaterPropsIAPWS& right) = delete;
    WaterPropsIAPWS& operator=(const WaterPropsIAPWS& right) = delete;

    //! Set the internal state of the object wrt temperature and density
    /*!
     * @param temperature   temperature (kelvin)
     * @param rho           density  (kg m-3)
     */
    void setState_TR(doublereal temperature, doublereal rho);

    //! Calculate the Helmholtz free energy in mks units of J kmol-1 K-1,
    //! using the last temperature and density
    doublereal helmholtzFE() const;

    //! Calculate the Gibbs free energy in mks units of J kmol-1 K-1.
    //! using the last temperature and density
    doublereal Gibbs() const;

    //! Calculate the enthalpy in mks units of J kmol-1
    //! using the last temperature and density
    doublereal enthalpy() const;

    //! Calculate the internal energy in mks units of J kmol-1
    doublereal intEnergy() const;

    //! Calculate the entropy in mks units of J kmol-1 K-1
    doublereal entropy() const;

    //! Calculate the constant volume heat capacity in mks units of J kmol-1 K-1
    //! at the last temperature and density
    doublereal cv() const;

    //! Calculate the constant pressure heat capacity in mks units of J kmol-1 K-1
    //! at the last temperature and density
    doublereal cp() const;

    //! Calculate the molar volume (kmol m-3) at the last temperature and
    //! density
    doublereal molarVolume() const;

    //! Calculates the pressure (Pascals), given the current value of the
    //! temperature and density.
    /*!
     * The density is an independent variable in the underlying equation of state
     *
     *  @returns the pressure (Pascal)
     */
    doublereal pressure() const;

    //! Calculates the density given the temperature and the pressure,
    //! and a guess at the density. Sets the internal state.
    /*!
     *  Note, below T_c, this is a multivalued function.
     *
     * The density() function calculates the density that is consistent with
     * a particular value of the temperature and pressure. It may therefore be
     * multivalued or potentially there may be no answer from this function.
     * It therefore takes a phase guess and a density guess as optional
     * parameters. If no guesses are supplied to density(), a gas phase guess
     * is assumed. This may or may not be what is wanted. Therefore, density()
     * should usually at least be supplied with a phase guess so that it may
     * manufacture an appropriate density guess. density() manufactures the
     * initial density guess, nondimensionalizes everything, and then calls
     * WaterPropsIAPWSphi::dfind(), which does the iterative calculation to
     * find the density condition that matches the desired input pressure.
     *
     * @param temperature  Kelvin
     * @param pressure     Pressure in Pascals (Newton/m**2)
     * @param phase        guessed phase of water; -1: no guessed phase
     * @param rhoguess     guessed density of the water; -1.0 no guessed density
     * @returns the density. If an error is encountered in the calculation the
     *     value of -1.0 is returned.
     */
    doublereal density(doublereal temperature, doublereal pressure,
                       int phase = -1, doublereal rhoguess = -1.0);

    //! Calculates the density given the temperature and the pressure,
    //! and a guess at the density, while not changing the internal state
    /*!
     * Note, below T_c, this is a multivalued function.
     *
     * The density() function calculates the density that is consistent with a
     * particular value of the temperature and pressure. It may therefore be
     * multivalued or potentially there may be no answer from this function.
     * It therefore takes a phase guess and a density guess as optional
     * parameters. If no guesses are supplied to density(), a gas phase guess
     * is assumed. This may or may not be what is wanted. Therefore, density()
     * should usually at least be supplied with a phase guess so that it may
     * manufacture an appropriate density guess. density() manufactures the
     * initial density guess, nondimensionalizes everything, and then calls
     * WaterPropsIAPWSphi::dfind(), which does the iterative calculation to
     * find the density condition that matches the desired input pressure.
     *
     * @param pressure    Pressure in Pascals (Newton/m**2)
     * @param phase       guessed phase of water; -1: no guessed phase
     * @param rhoguess    guessed density of the water; -1.0: no guessed density
     * @returns the density. If an error is encountered in the calculation the
     *     value of -1.0 is returned.
     */
    doublereal density_const(doublereal pressure, int phase = -1, doublereal rhoguess = -1.0) const;

    //! Returns the density (kg m-3)
    /*!
     * The density is an independent variable in the underlying equation of state
     *
     * @returns the density (kg m-3)
     */
    doublereal density() const;

    //! Returns the temperature (Kelvin)
    /*!
     * @return s the internally stored temperature
     */
    doublereal temperature() const;

    //! Returns the coefficient of thermal expansion.
    /*!
     * alpha = d (ln V) / dT at constant P.
     *
     * @returns the coefficient of thermal expansion
     */
    doublereal coeffThermExp() const;

    //! Returns the isochoric pressure derivative wrt temperature
    /*!
     * beta = M / (rho * Rgas) (d (pressure) / dT) at constant rho
     *
     * Note for ideal gases this is equal to one.
     *
     *     beta = delta (phi0_d() + phiR_d()) - tau delta (phi0_dt() + phiR_dt())
     */
    doublereal coeffPresExp() const;

    //! Returns the coefficient of isothermal compressibility for the state of
    //! the object
    /*!
     * kappa = - d (ln V) / dP at constant T.
     *
     * units - 1/Pascal
     *
     * @returns the isothermal compressibility
     */
    doublereal isothermalCompressibility() const;

    //! Returns the value of dp / drho at constant T for the state of the object
    /*!
     *  units - Joules / kg
     *
     * @returns dpdrho
     */
    doublereal dpdrho() const;

    //! This function returns an estimated value for the saturation pressure.
    /*!
     * It does this via a polynomial fit of the vapor pressure curve.
     * units = (Pascals)
     *
     * @param temperature Input temperature (Kelvin)
     *
     * @returns the estimated saturation pressure
     */
    doublereal psat_est(doublereal temperature) const;

    //! This function returns the saturation pressure given the temperature as
    //! an input parameter, and sets the internal state to the saturated
    //! conditions.
    /*!
     * Note this function will return the saturation pressure, given the
     * temperature. It will then set the state of the system to the saturation
     * condition. The input parameter waterState is used to either specify the
     * liquid state or the gas state at the desired temperature and saturated
     * pressure.
     *
     * If the input temperature, T, is above T_c, this routine will set the
     * internal state to T and the pressure to P_c. Then, return P_c.
     *
     * @param temperature   input temperature (kelvin)
     * @param waterState    integer specifying the water state
     * @returns the saturation pressure. units = Pascal
     */
    doublereal psat(doublereal temperature, int waterState = WATER_LIQUID);

    //! Return the value of the density at the water spinodal point (on the
    //! liquid side) for the current temperature.
    /*!
     * @returns the density with units of kg m-3
     */
    doublereal densSpinodalWater() const;

    //! Return the value of the density at the water spinodal point (on the gas
    //! side) for the current temperature.
    /*!
     * @returns the density with units of kg m-3
     */
    doublereal densSpinodalSteam() const;

    //! Returns the Phase State flag for the current state of the object
    /*!
     * @param checkState If true, this function does a complete check to see
     *        where in parameters space we are
     *
     *  There are three values:
     *    - WATER_GAS   below the critical temperature but below the critical density
     *    - WATER_LIQUID  below the critical temperature but above the critical density
     *    - WATER_SUPERCRIT   above the critical temperature
     */
    int phaseState(bool checkState = false) const;

    //! Returns the critical temperature of water (Kelvin)
    /*!
     * This is hard coded to the value 647.096 Kelvin
     */
    doublereal Tcrit() const {
        return 647.096;
    }

    //! Returns the critical pressure of water (22.064E6 Pa)
    /*!
     * This is hard coded to the value of 22.064E6 pascals
     */
    doublereal Pcrit() const {
        return 22.064E6;
    }

    //! Return the critical density of water (kg m-3)
    /*!
     * This is equal to 322 kg m-3.
     */
    doublereal Rhocrit() const {
        return 322.;
    }

private:
    //! Calculate the dimensionless temp and rho and store internally.
    /*!
     * @param temperature   input temperature (kelvin)
     * @param rho           density in kg m-3
     */
    void calcDim(doublereal temperature, doublereal rho);

    //! Utility routine in the calculation of the saturation pressure
    /*!
     * Calculate the Gibbs free energy in mks units of J kmol-1 K-1.
     *
     * @param temperature    temperature (kelvin)
     * @param pressure       pressure (Pascal)
     * @param densLiq        Output density of liquid
     * @param densGas        output Density of gas
     * @param delGRT         output delGRT
     */
    void corr(doublereal temperature, doublereal pressure, doublereal& densLiq,
              doublereal& densGas, doublereal& delGRT);

    //! Utility routine in the calculation of the saturation pressure
    /*!
     * @param temperature    temperature (kelvin)
     * @param pressure       pressure (Pascal)
     * @param densLiq        Output density of liquid
     * @param densGas        output Density of gas
     * @param pcorr          output corrected pressure
     */
    void corr1(doublereal temperature, doublereal pressure, doublereal& densLiq,
               doublereal& densGas, doublereal& pcorr);

    //! pointer to the underlying object that does the calculations.
    mutable WaterPropsIAPWSphi m_phi;

    //! Dimensionless temperature,  tau = T_C / T
    doublereal tau;

    //! Dimensionless density, delta = rho / rho_c
    mutable doublereal delta;

    //! Current state of the system
    mutable int iState;
};

}
#endif
