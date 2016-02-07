/**
 *  @file ConstCpPoly.h
 * Headers for the \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink
 * object that employs a constant heat capacity assumption (see \ref spthermo and
 * \link Cantera::ConstCpPoly ConstCpPoly\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_CONSTCPPOLY_H
#define CT_CONSTCPPOLY_H

#include "cantera/thermo/SpeciesThermoInterpType.h"

namespace Cantera
{

/**
 * A constant-heat capacity species thermodynamic property manager class. This
 * makes the assumption that the heat capacity is a constant. Then, the
 * following relations are used to complete the specification of the
 * thermodynamic functions for the species.
 *
 * \f[
 * \frac{c_p(T)}{R} = Cp0\_R
 * \f]
 * \f[
 * \frac{h^0(T)}{RT} = \frac{1}{T} * (h0\_R + (T - T_0) * Cp0\_R)
 * \f]
 * \f[
 * \frac{s^0(T)}{R} =  (s0\_R + (log(T) - log(T_0)) * Cp0\_R)
 * \f]
 *
 * This parameterization takes 4 input values. These are:
 *       -   c[0] = \f$ T_0 \f$(Kelvin)
 *       -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
 *       -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
 *       -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
 *
 * @ingroup spthermo
 */
class ConstCpPoly: public SpeciesThermoInterpType
{
public:
    //! empty constructor
    ConstCpPoly();

    //! Normal constructor
    /*!
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the parameters for
     *                     the standard state for species n. There are 4
     *                     coefficients for the ConstCpPoly parameterization.
     *           -   c[0] = \f$ T_0 \f$(Kelvin)
     *           -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
     *           -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
     *           -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
     */
    ConstCpPoly(double tlow, double thigh, double pref, const double* coeffs);

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const;

    virtual int reportType() const {
        return CONSTANT_CP;
    }

    /*!
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Form and Length of the temperature polynomial:
     *  - m_t[0] = tt;
     *
     */
    void updateProperties(const doublereal* tt,
                          doublereal* cp_R, doublereal* h_RT,
                          doublereal* s_R) const;

    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R, doublereal* h_RT,
                              doublereal* s_R) const;
    void reportParameters(size_t& n, int& type,
                          doublereal& tlow, doublereal& thigh,
                          doublereal& pref,
                          doublereal* const coeffs) const;
    virtual void modifyParameters(doublereal* coeffs);

    virtual doublereal reportHf298(doublereal* const h298 = 0) const;
    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New);

protected:
    //! Base temperature
    doublereal m_t0;
    //! Dimensionless value of the heat capacity
    doublereal m_cp0_R;
    //! dimensionless value of the enthaply at t0
    doublereal m_h0_R;
    //! Dimensionless value of the entropy at t0
    doublereal m_s0_R;
    //! log of the t0 value
    doublereal m_logt0;
};

}

#endif
