/**
 *  @file Mu0Poly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on a piecewise constant mu0 interpolation
 *  (see \ref spthermo and class \link Cantera::Mu0Poly Mu0Poly\endlink).
 */
#ifndef CT_MU0POLY_H
#define CT_MU0POLY_H

#include "cantera/thermo/SpeciesThermoInterpType.h"

namespace Cantera
{
template<typename ValAndDerivType> class SpeciesThermo;
class XML_Node;

//!  The %Mu0Poly class implements an interpolation of the Gibbs free energy based on a
//!  piecewise constant heat capacity approximation.
/*!
 *   The %Mu0Poly class implements a piecewise constant heat capacity approximation.
 *   of the standard state chemical potential of one
 *   species at a single reference pressure.
 *   The chemical potential is input as a series of (\f$T\f$, \f$ \mu^o(T)\f$)
 *   values. The first temperature is assumed to be equal
 *   to 298.15 K; however, this may be relaxed in the future.
 *   This information, and an assumption of a constant
 *   heat capacity within each interval is enough to
 *   calculate all thermodynamic functions.
 *
 *  The piece-wise constant heat capacity is calculated from the change in the chemical potential over each interval.
 *  Once the heat capacity is known, the other thermodynamic functions may be determined.
 *  The basic equation for going from temperature point 1 to temperature point 2
 *  are as follows for \f$ T \f$,  \f$ T_1 <= T <= T_2 \f$
 *
 * \f[
 *      \mu^o(T_1) = h^o(T_1) - T_1 * s^o(T_1)
 * \f]
 * \f[
 *      \mu^o(T_2) - \mu^o(T_1) = Cp^o(T_1)(T_2 - T_1) - Cp^o(T_1)(T_2)ln(\frac{T_2}{T_1}) - s^o(T_1)(T_2 - T_1)
 * \f]
 * \f[
 *      s^o(T_2) = s^o(T_1) + Cp^o(T_1)ln(\frac{T_2}{T_1})
 * \f]
 * \f[
 *      h^o(T_2) = h^o(T_1) + Cp^o(T_1)(T_2 - T_1)
 * \f]
 *
 *  Within each interval the following relations are used. For \f$ T \f$,  \f$ T_1 <= T <= T_2 \f$
 *
 * \f[
 *      \mu^o(T) = \mu^o(T_1) + Cp^o(T_1)(T - T_1) - Cp^o(T_1)(T_2)ln(\frac{T}{T_1}) - s^o(T_1)(T - T_1)
 * \f]
 * \f[
 *      s^o(T) = s^o(T_1) + Cp^o(T_1)ln(\frac{T}{T_1})
 * \f]
 * \f[
 *      h^o(T) = h^o(T_1) + Cp^o(T_1)(T - T_1)
 * \f]
 *
 *   Notes about temperature interpolation for \f$ T < T_1 \f$ and \f$ T > T_{npoints} \f$.
 *     These are achieved by assuming a constant heat capacity
 *     equal to the value in the closest temperature interval.
 *     No error is thrown.
 *
 *   @note In the future, a better assumption about the heat
 *         capacity may be employed, so that it can be continuous.
 *
 * @ingroup spthermo
 */
template<typename ValAndDerivType>
class Mu0Poly: public SpeciesThermoInterpType<ValAndDerivType>
{

public:

    //! Constructor
    Mu0Poly();

    //! Constructor used in templated instantiations
    /*!
     *
     * In the constructor, we calculate and store the
     * piecewise linear approximation to the thermodynamic
     * functions.
     *
     * @param n            Species index
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the
     *                     parameters for the standard state for species n.
     *                     There are \f$ 2+npoints*2 \f$ coefficients, where
     *                     \f$ npoints \f$ are the number of temperature points.
     *                     Their identity is further broken down:
     *            - coeffs[0] = number of points (integer)
     *            - coeffs[1]  = \f$ h^o(298.15 K) \f$ (J/kmol)
     *            - coeffs[2]  = \f$ T_1 \f$  (Kelvin)
     *            - coeffs[3]  = \f$ \mu^o(T_1) \f$ (J/kmol)
     *            - coeffs[4]  = \f$ T_2 \f$  (Kelvin)
     *            - coeffs[5]  = \f$ \mu^o(T_2) \f$ (J/kmol)
     *            - coeffs[6]  = \f$ T_3 \f$  (Kelvin)
     *            - coeffs[7]  = \f$ \mu^o(T_3) \f$ (J/kmol)
     *            - ........
     *            .
     */
    Mu0Poly(size_t n, doublereal tlow, doublereal thigh,
            doublereal pref, const doublereal* coeffs);

    //! Copy constructor
    Mu0Poly(const Mu0Poly&);

    //! Assignment operator
    Mu0Poly& operator=(const Mu0Poly&);

    //! Destructor
    virtual ~Mu0Poly();

    //! Duplicator
    virtual SpeciesThermoInterpType<ValAndDerivType>*
    duplMyselfAsSpeciesThermoInterpType() const;

    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    virtual doublereal minTemp() const;

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    virtual doublereal maxTemp() const;

    //! Returns the reference pressure (Pa)
    virtual doublereal refPressure() const;

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const {
        return MU0_INTERP;
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
     * Temperature Polynomial:
     *
     * tPoly[0] = temp (Kelvin)
     *
     * @param tPoly  vector of temperature polynomials. Length = 1
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void updateProperties(const doublereal* tPoly,
                                  ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                  ValAndDerivType* s_R) const ;

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
    virtual void updatePropertiesTemp(const doublereal temp,
                                      ValAndDerivType* cp_R,
                                      ValAndDerivType* h_RT,
                                      ValAndDerivType* s_R) const ;

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
    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const;

    //! Modify parameters for the standard state
    /*!
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     */
    virtual void modifyParameters(doublereal* coeffs);

protected:

    /**
     * Number of intervals in the interpolating linear
     * approximation. Number of points is one more than the
     * number of intervals.
     */
    size_t m_numIntervals;

    /**
     * Value of the enthalpy at T = 298.15.
     *  This value is tied to the Heat of formation of
     *  the species at 298.15.
     */
    doublereal m_H298;

    /**
     * Points at which the standard state chemical potential
     * are given.
     */
    vector_fp m_t0_int;

    /**
     * Mu0's are primary input data. They aren't strictly
     * needed, but are kept here for convenience.
     */
    vector_fp m_mu0_R_int;

    //! Dimensionless Enthalpies at the temperature points
    vector_fp m_h0_R_int;

    //! Entropy at the points
    vector_fp m_s0_R_int;

    //! Heat capacity at the points
    vector_fp m_cp0_R_int;
    //! Limiting low temperature
    doublereal m_lowT;
    //! Limiting high temperature
    doublereal m_highT;

    //! Reference pressure
    doublereal m_Pref;

    //! Species index
    size_t m_index;

private:

    //! process the coefficients
    /*!
     * Mu0Poly():
     *
     * In the constructor, we calculate and store the
     * piecewise linear approximation to the thermodynamic
     * functions.
     *
     * @param coeffs coefficients. These are defined as follows:
     *
     *  coeffs[0] = number of points (integer)
     *         1  = H298(J/kmol)
     *         2  = T1  (Kelvin)
     *         3  = mu1 (J/kmol)
     *         4  = T2  (Kelvin)
     *         5  = mu2 (J/kmol)
     *         6  = T3  (Kelvin)
     *         7  = mu3 (J/kmol)
     *         ........
     */
    void processCoeffs(const doublereal* coeffs);

};

//!  Install a Mu0 polynomial thermodynamic reference state
/*!
 * Install a Mu0 polynomial thermodynamic reference state property
 * parameterization for species k into a SpeciesThermo instance,
 * getting the information from an XML database.
 *
 * @param speciesName  Name of the species
 * @param sp           Owning SpeciesThermo object
 * @param k            Species index
 * @param Mu0Node_ptr  Pointer to the XML element containing the
 *                     Mu0 information.
 *
 *  @ingroup spthermo
 */
template<typename ValAndDerivType>
void installMu0ThermoFromXML(const std::string& speciesName,
                             SpeciesThermo<ValAndDerivType>& sp, size_t k,
                             const XML_Node* Mu0Node_ptr);
}

#endif


