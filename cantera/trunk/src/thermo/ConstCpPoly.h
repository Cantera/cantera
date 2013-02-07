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
 *
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
     *                     There are 4 coefficients for the %ConstCpPoly parameterization.
     *           -   c[0] = \f$ T_0 \f$(Kelvin)
     *           -   c[1] = \f$ H_k^o(T_0, p_{ref}) \f$ (J/kmol)
     *           -   c[2] = \f$ S_k^o(T_0, p_{ref}) \f$    (J/kmol K)
     *           -   c[3] = \f$ {Cp}_k^o(T_0, p_{ref}) \f$  (J(kmol K)
     *
     */
    ConstCpPoly(size_t n, doublereal tlow, doublereal thigh,
                doublereal pref,
                const doublereal* coeffs);

    //! copy constructor
    ConstCpPoly(const ConstCpPoly&);

    //! Assignment operator
    ConstCpPoly& operator=(const ConstCpPoly&);

    //! Destructor
    virtual ~ConstCpPoly();

    //! Duplicator
    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const;
    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    doublereal minTemp() const;

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    doublereal maxTemp() const;

    //! Returns the reference pressure (Pa)
    doublereal refPressure() const;

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const {
        return CONSTANT_CP;
    }

    //! Returns an integer representing the species index
    virtual size_t speciesIndex() const {
        return m_index;
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
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    void updateProperties(const doublereal* tt,
                          doublereal* cp_R, doublereal* h_RT,
                          doublereal* s_R) const;

    //! Compute the reference-state property of one species
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of one of the species. The species index is used
     * to reference into the cp_R, h_RT, and s_R arrays.
     *
     * @param temp    Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R, doublereal* h_RT,
                              doublereal* s_R) const;
    //!This utility function reports back the type of
    //! parameterization and all of the parameters for the
    //! species, index.
    /*!
     * All parameters are output variables
     *
     * @param n         Species index
     * @param type      Integer type of the standard type
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the
     *                  parameters for the standard state.
     */
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

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportHf298(doublereal* const h298 = 0) const;

    virtual void modifyOneHf298(const size_t& k, const doublereal Hf298New);

#endif

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
    //! Minimum temperature for which the parameterization is valid (Kelvin)
    doublereal m_lowT;
    //! Maximum temperature for which the parameterization is valid (Kelvin)
    doublereal m_highT;
    //! Reference pressure (Pa)
    doublereal m_Pref;
    //! Species Index
    size_t m_index;

private:

};

}

#endif
