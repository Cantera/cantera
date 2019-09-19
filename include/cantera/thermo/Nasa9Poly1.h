/**
 * @file Nasa9Poly1.h Header for a single-species standard state object derived
 *     from \link Cantera::SpeciesThermoInterpType
 *     SpeciesThermoInterpType\endlink based on the NASA 9 coefficient
 *     temperature polynomial form applied to one temperature region (see \ref
 *     spthermo and class \link Cantera::Nasa9Poly1 Nasa9Poly1\endlink).
 *
 * This parameterization has one NASA temperature region.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NASA9POLY1_H
#define CT_NASA9POLY1_H

#include "SpeciesThermoInterpType.h"

namespace Cantera
{
//! The NASA 9 polynomial parameterization for one temperature range.
/*!
 * This parameterization expresses the heat capacity via a 7 coefficient
 * polynomial. Note that this is the form used in the 2002 NASA equilibrium
 * program. A reference to the form is provided below:
 *
 * "NASA Glenn Coefficients for Calculating Thermodynamic Properties of
 * Individual Species," B. J. McBride, M. J. Zehe, S. Gordon
 * NASA/TP-2002-211556, Sept. 2002
 *
 * Nine coefficients \f$(a_0,\dots,a_8)\f$ are used to represent
 * \f$ C_p^0(T)\f$, \f$ H^0(T)\f$, and \f$ S^0(T) \f$ as
 * polynomials in \f$ T \f$ :
 * \f[
 * \frac{C_p^0(T)}{R} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T
 *                  + a_4 T^2 + a_5 T^3 + a_6 T^4
 * \f]
 *
 * \f[
 * \frac{H^0(T)}{RT} = - a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2
 *     + \frac{a_3}{2} T + \frac{a_4}{3} T^2  + \frac{a_5}{4} T^3 +
 *     \frac{a_6}{5} T^4 + \frac{a_7}{T}
 * \f]
 *
 * \f[
 * \frac{s^0(T)}{R} = - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T
 *    + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3  + \frac{a_6}{4} T^4 + a_8
 * \f]
 *
 * The standard state is assumed to be an ideal gas at the standard pressure of
 * 1 bar, for gases. For condensed species, the standard state is the pure
 * crystalline or liquid substance at the standard pressure of 1 atm.
 *
 * These NASA representations may have multiple temperature regions through the
 * use of the Nasa9PolyMultiTempRegion object, which uses multiple copies of
 * this Nasa9Poly1 object to handle multiple temperature regions.
 *
 * @ingroup spthermo
 * @see Nasa9PolyMultiTempRegion
 */
class Nasa9Poly1 : public SpeciesThermoInterpType
{
public:
    Nasa9Poly1();

    //! Constructor with all input data
    /*!
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the
     *                     parameters for the standard state.
     */
    Nasa9Poly1(double tlow, double thigh, double pref, const double* coeffs);

    //! Set the array of 9 polynomial coefficients
    void setParameters(const vector_fp& coeffs);

    virtual int reportType() const;

    virtual size_t temperaturePolySize() const { return 7; }
    virtual void updateTemperaturePoly(double T, double* T_poly) const;

    /*!
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Temperature Polynomial:
     *   - tt[0] = t;
     *   - tt[1] = t*t;
     *   - tt[2] = t*t*t;
     *   - tt[3] = t*t*t*t;
     *   - tt[4] = 1.0/t;
     *   - tt[5] = 1.0/(t*t);
     *   - tt[6] = std::log(t);
     */
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const;

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const;

    //! This utility function reports back the type of parameterization and all
    //! of the parameters for the species
    /*!
     * All parameters are output variables
     *
     * @param n         Species index
     * @param type      Integer type of the standard type
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the parameters for
     *     the standard state. There are 12 of them, designed to be compatible
     *     with the multiple temperature formulation.
     *       - coeffs[0] is equal to one.
     *       - coeffs[1] is min temperature
     *       - coeffs[2] is max temperature
     *       - coeffs[3+i] from i =0,9 are the coefficients themselves
     */
    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const;

protected:
    //! array of polynomial coefficients
    vector_fp m_coeff;
};

}
#endif
