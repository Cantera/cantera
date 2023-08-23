/**
 *  @file ConstCpPoly.h
 * Headers for the @link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType@endlink
 * object that employs a constant heat capacity assumption (see @ref spthermo and
 * @link Cantera::ConstCpPoly ConstCpPoly@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTCPPOLY_H
#define CT_CONSTCPPOLY_H

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/speciesThermoTypes.h"

namespace Cantera
{

/**
 * A constant-heat capacity species thermodynamic property manager class. This
 * makes the assumption that the heat capacity is a constant. Then, the
 * following relations are used to complete the specification of the
 * thermodynamic functions for the species.
 *
 * @f[
 *  \hat{c}_p^\circ(T) = \hat{c}_p^\circ(T^\circ)
 * @f]
 * @f[
 *  \hat{h}^\circ(T) = \hat{h}^\circ\left(T_0\right)
 *      + \hat{c}_p^\circ \left(T-T^\circ\right)
 * @f]
 * @f[
 *  \hat{s}^\circ(T) = \hat{s}^\circ(T_0)
 *      + \hat{c}_p^\circ \ln{\left(\frac{T}{T^\circ}\right)}
 * @f]
 *
 * This parameterization takes 4 input values @f$ T^\circ @f$,
 * @f$ \hat{h}^\circ(T^\circ) @f$, @f$ \hat{s}^\circ(T^\circ) @f$ and
 * @f$ \hat{c}_p^\circ(T^\circ) @f$, see setParameters(). The default value of
 * @f$ T^\circ @f$ is 298.15 K; the default value for the other parameters is 0.0.
 * @ingroup spthermo
 */
class ConstCpPoly: public SpeciesThermoInterpType
{
public:
    ConstCpPoly();

    //! Constructor with all input data
    /*!
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the parameters for
     *                     the standard state for species n. Contains 4 parameters
     *                     in the order of setParameters() arguments.
     */
    ConstCpPoly(double tlow, double thigh, double pref, const double* coeffs);

    /**
     * Set ConstCpPoly parameters.
     * @param t0  @f$ T^\circ @f$ [K]
     * @param h0  @f$ \hat{h}^\circ(T^\circ) @f$ [J/kmol]
     * @param s0  @f$ \hat{s}^\circ(T^\circ) @f$ [J/kmol/K]
     * @param cp0  @f$ \hat{c}_p^\circ(T^\circ) @f$ [J/kmol/K]
     */
    void setParameters(double t0, double h0, double s0, double cp0);

    int reportType() const override {
        return CONSTANT_CP;
    }

    /**
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Form and Length of the temperature polynomial:
     *  - m_t[0] = tt;
     *
     */
    void updateProperties(const double* tt, double* cp_R, double* h_RT,
                          double* s_R) const override;

    void updatePropertiesTemp(const double temp, double* cp_R, double* h_RT,
                              double* s_R) const override;

    size_t nCoeffs() const override { return 4; }

    void reportParameters(size_t& n, int& type, double& tlow, double& thigh,
                          double& pref, double* const coeffs) const override;

    void getParameters(AnyMap& thermo) const override;

    double reportHf298(double* const h298=nullptr) const override;
    void modifyOneHf298(const size_t k, const double Hf298New) override;
    void resetHf298() override;

protected:
    //! Base temperature
    double m_t0 = 298.15;
    //! Dimensionless value of the heat capacity
    double m_cp0_R = 0.0;
    //! dimensionless value of the enthalpy at t0
    double m_h0_R = 0.0;
    //! Dimensionless value of the entropy at t0
    double m_s0_R = 0.0;
    //! log of the t0 value
    double m_logt0 = log(298.15);
    //! Original value of h0_R, restored by calling resetHf298()
    double m_h0_R_orig = 0.0;
};

}

#endif
