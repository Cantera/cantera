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
 *  A constant-heat capacity species thermodynamic property manager class.
 *  This makes the
 *  assumption that the heat capacity is a constant. Then, the following
 *  relations are used to complete the specification of the thermodynamic
 *  functions for the species.
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
 * The multispecies SimpleThermo class makes the same assumptions as
 * this class does.
 *
 * @see SimpleThermo
 * @ingroup spthermo
 */
class ConstCpPoly: public SpeciesThermoInterpType
{
public:
    //! empty constructor
    ConstCpPoly();

    //! Constructor used in templated instantiations
    /*!
     * @param n            Species index
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the
     *                     parameters for the standard state for species n.
     *                     There are 4 coefficients for the ConstCpPoly parameterization.
     *           -   c[0] = \f$ T_0 \f$(Kelvin)
     *           -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
     *           -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
     *           -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
     * @deprecated Use the constructor which does not take the species index. To
     *     be removed after Cantera 2.2.
     */
    ConstCpPoly(size_t n, doublereal tlow, doublereal thigh,
                doublereal pref,
                const doublereal* coeffs);

    //! Normal constructor
    /*!
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the
     *                     parameters for the standard state for species n.
     *                     There are 4 coefficients for the ConstCpPoly parameterization.
     *           -   c[0] = \f$ T_0 \f$(Kelvin)
     *           -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
     *           -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
     *           -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
     */
    ConstCpPoly(double tlow, double thigh, double pref, const double* coeffs);

    //! copy constructor
    ConstCpPoly(const ConstCpPoly&);

    //! Assignment operator
    ConstCpPoly& operator=(const ConstCpPoly&);

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const;

    virtual int reportType() const {
        return CONSTANT_CP;
    }

    //! Update the properties for this species, given a temperature polynomial
    /*!
     * This method is called with a pointer to an array containing the functions of
     * temperature needed by this  parameterization, and three pointers to arrays where the
     * computed property values should be written. This method updates only one value in
     * each array.
     *
     * Form and Length of the temperature polynomial:
     *  - m_t[0] = tt;
     *
     * @param tt      Vector of temperature polynomials
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
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
    //! Modify parameters for the standard state
    /*!
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     */
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
