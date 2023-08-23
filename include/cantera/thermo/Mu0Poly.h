/**
 *  @file Mu0Poly.h
 *  Header for a single-species standard state object derived
 *  from @link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType@endlink  based
 *  on a piecewise constant mu0 interpolation
 *  (see @ref spthermo and class @link Cantera::Mu0Poly Mu0Poly@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MU0POLY_H
#define CT_MU0POLY_H

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/speciesThermoTypes.h"

namespace Cantera
{

//! The Mu0Poly class implements an interpolation of the Gibbs free energy based
//! on a piecewise constant heat capacity approximation.
/*!
 * The Mu0Poly class implements a piecewise constant heat capacity
 * approximation. of the standard state chemical potential of one species at a
 * single reference pressure. The chemical potential is input as a series of
 * (@f$ T @f$, @f$ \mu^o(T) @f$) values. The first temperature is assumed to be
 * equal to 298.15 K; however, this may be relaxed in the future. This
 * information, and an assumption of a constant heat capacity within each
 * interval is enough to calculate all thermodynamic functions.
 *
 * The piece-wise constant heat capacity is calculated from the change in the
 * chemical potential over each interval. Once the heat capacity is known, the
 * other thermodynamic functions may be determined. The basic equation for going
 * from temperature point 1 to temperature point 2 are as follows for @f$ T @f$,
 * @f$ T_1 <= T <= T_2 @f$
 *
 * @f[
 *      \mu^o(T_1) = h^o(T_1) - T_1 * s^o(T_1)
 * @f]
 * @f[
 *      \mu^o(T_2) - \mu^o(T_1) = Cp^o(T_1)(T_2 - T_1) - Cp^o(T_1)(T_2)ln(\frac{T_2}{T_1}) - s^o(T_1)(T_2 - T_1)
 * @f]
 * @f[
 *      s^o(T_2) = s^o(T_1) + Cp^o(T_1)ln(\frac{T_2}{T_1})
 * @f]
 * @f[
 *      h^o(T_2) = h^o(T_1) + Cp^o(T_1)(T_2 - T_1)
 * @f]
 *
 *  Within each interval the following relations are used. For @f$ T @f$, @f$
 *  T_1 <= T <= T_2 @f$
 *
 * @f[
 *      \mu^o(T) = \mu^o(T_1) + Cp^o(T_1)(T - T_1) - Cp^o(T_1)(T_2)ln(\frac{T}{T_1}) - s^o(T_1)(T - T_1)
 * @f]
 * @f[
 *      s^o(T) = s^o(T_1) + Cp^o(T_1)ln(\frac{T}{T_1})
 * @f]
 * @f[
 *      h^o(T) = h^o(T_1) + Cp^o(T_1)(T - T_1)
 * @f]
 *
 * Notes about temperature interpolation for @f$ T < T_1 @f$ and @f$ T >
 * T_{npoints} @f$: These are achieved by assuming a constant heat capacity
 * equal to the value in the closest temperature interval. No error is thrown.
 *
 * @note In the future, a better assumption about the heat capacity may be
 *       employed, so that it can be continuous.
 *
 * @ingroup spthermo
 */
class Mu0Poly: public SpeciesThermoInterpType
{
public:
    Mu0Poly();

    //! Constructor with all input data
    /*!
     * @param tlow    Minimum temperature
     * @param thigh   Maximum temperature
     * @param pref    reference pressure (Pa).
     * @param coeffs  Vector of coefficients used to set the parameters for the
     *                standard state for species n. There are @f$ 2+npoints*2
     *                @f$ coefficients, where @f$ npoints @f$ are the number of
     *                temperature points. Their identity is further broken down:
     *            - coeffs[0] = number of points (integer)
     *            - coeffs[1] = @f$ h^o(298.15 K) @f$ (J/kmol)
     *            - coeffs[2] = @f$ T_1 @f$  (Kelvin)
     *            - coeffs[3] = @f$ \mu^o(T_1) @f$ (J/kmol)
     *            - coeffs[4] = @f$ T_2 @f$  (Kelvin)
     *            - coeffs[5] = @f$ \mu^o(T_2) @f$ (J/kmol)
     *            - coeffs[6] = @f$ T_3 @f$  (Kelvin)
     *            - coeffs[7] = @f$ \mu^o(T_3) @f$ (J/kmol)
     *            - ........
     */
    Mu0Poly(double tlow, double thigh, double pref, const double* coeffs);

    //! Set parameters for @f$ \mu^o(T) @f$
    /*!
     * Calculates and stores the piecewise linear approximation to the
     * thermodynamic functions.
     *
     *  @param h0    Enthalpy at the reference temperature of 298.15 K [J/kmol]
     *  @param T_mu  Map with temperature [K] as the keys and the Gibbs free
     *               energy [J/kmol] as the values. Must contain one point at
     *               298.15 K.
     */
    void setParameters(double h0, const map<double, double>& T_mu);

    int reportType() const override {
        return MU0_INTERP;
    }

    /**
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Temperature Polynomial:
     *     tt[0] = temp (Kelvin)
     */
    void updateProperties(const double* tt, double* cp_R, double* h_RT,
                          double* s_R) const override;

    void updatePropertiesTemp(const double temp, double* cp_R, double* h_RT,
                              double* s_R) const override;

    size_t nCoeffs() const override;

    void reportParameters(size_t& n, int& type, double& tlow, double& thigh,
                          double& pref, double* const coeffs) const override;

    void getParameters(AnyMap& thermo) const override;

protected:
    //! Number of intervals in the interpolating linear approximation. Number
    //! of points is one more than the number of intervals.
    size_t m_numIntervals = 0;

    //! Value of the enthalpy at T = 298.15. This value is tied to the Heat of
    //! formation of the species at 298.15.
    double m_H298 = 0.0;

    //! Points at which the standard state chemical potential are given.
    vector<double> m_t0_int;

    //! Mu0's are primary input data. They aren't strictly needed, but are kept
    //! here for convenience.
    vector<double> m_mu0_R_int;

    //! Dimensionless Enthalpies at the temperature points
    vector<double> m_h0_R_int;

    //! Entropy at the points
    vector<double> m_s0_R_int;

    //! Heat capacity at the points
    vector<double> m_cp0_R_int;
};

}

#endif
