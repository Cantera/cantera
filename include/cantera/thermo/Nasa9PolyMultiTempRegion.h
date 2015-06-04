/**
 *  @file Nasa9PolyMultiTempRegion.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType
 *  SpeciesThermoInterpType\endlink  based
 *  on the NASA 9 coefficient temperature polynomial form
 *   applied to multiple temperature regions
 *  (see \ref spthermo and class \link Cantera::Nasa9PolyMultiTempRegion Nasa9PolyMultiTempRegion\endlink).
 *
 *  This parameterization has multiple NASA temperature regions.
 */

#ifndef CT_NASA9POLYMULTITEMPREGION_H
#define CT_NASA9POLYMULTITEMPREGION_H
// Copyright 2007  Sandia National Laboratories

#include "cantera/thermo/Nasa9Poly1.h"

namespace Cantera
{
//! The NASA 9 polynomial parameterization for a single species
//! encompassing multiple temperature regions.
/*!
 * The parameterization used in each temperature range is described in the
 * documentation for class Nasa9Poly1.
 *
 * These NASA representations may have multiple temperature regions
 * through the use of this Nasa9PolyMultiTempRegion object, which uses
 * multiple copies of the Nasa9Poly1 object to handle multiple temperature
 * regions.
 *
 * @ingroup spthermo
 * @see Nasa9Poly1
 */
class Nasa9PolyMultiTempRegion : public SpeciesThermoInterpType
{
public:
    //! Empty constructor
    Nasa9PolyMultiTempRegion();

    //! Constructor used in templated instantiations
    /*!
     * @param regionPts Vector of pointers to Nasa9Poly1 objects. These
     *                  objects all refer to the temperature regions for the
     *                  same species. The vector must be in increasing
     *                  temperature region format.  Together they
     *                  represent the reference temperature parameterization
     *                  for a single species.
     *
     *  Note, after the constructor, we will own the underlying
     *  Nasa9Poly1 objects and be responsible for owning them.
     */
    Nasa9PolyMultiTempRegion(std::vector<Cantera::Nasa9Poly1*> &regionPts);

    //! Copy constructor
    /*!
     * @param b object to be copied
     */
    Nasa9PolyMultiTempRegion(const Nasa9PolyMultiTempRegion& b);

    //! Assignment operator
    /*!
     * @param b object to be copied
     */
    Nasa9PolyMultiTempRegion& operator=(const Nasa9PolyMultiTempRegion& b);

    //! Destructor
    virtual ~Nasa9PolyMultiTempRegion();

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const;

    virtual int reportType() const;

    virtual void setIndex(size_t index);

    virtual size_t temperaturePolySize() const { return 7; }
    virtual void updateTemperaturePoly(double T, double* T_poly) const;

    //! Update the properties for this species, given a temperature polynomial
    /*!
     * This method is called with a pointer to an array containing the
     * functions of temperature needed by this  parameterization, and three
     * pointers to arrays where the computed property values should be
     * written. This method updates only one value in each array.
     *
     * Temperature Polynomial:
     *  tt[0] = t;
     *  tt[1] = t*t;
     *  tt[2] = t*t*t;
     *  tt[3] = t*t*t*t;
     *  tt[4] = 1.0/t;
     *  tt[5] = 1.0/(t*t);
     *  tt[6] = std::log(t);
     *
     * @param tt      vector of temperature polynomials
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updateProperties(const doublereal* tt,
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
     * Temperature Polynomial:
     *  tt[0] = t;
     *  tt[1] = t*t;
     *  tt[2] = t*t*t;
     *  tt[3] = t*t*t*t;
     *  tt[4] = 1.0/t;
     *  tt[5] = 1.0/(t*t);
     *  tt[6] = std::log(t);
     *
     * @param temp    Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updatePropertiesTemp(const doublereal temp,
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
     *      There are 1 + 11*nzones coefficients
     *      coeffs[0] is equal to nTempZones.
     *      index = 1
     *      for each zone:
     *        coeffs[index] = minTempZone
     *        coeffs[index+1] = maxTempZone
     *        coeffs[index+2+i] from i =0,9 are the coefficients themselves
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
    //! Number of temperature regions
    size_t m_numTempRegions;

    //! Lower boundaries of each temperature regions
    vector_fp m_lowerTempBounds;

    //! pointers to the objects
    /*!
     * This object will now own these pointers and delete
     * them when the current object is deleted.
     */
    std::vector<Nasa9Poly1*>m_regionPts;

    //! current region
    mutable int m_currRegion;
};

}
#endif
