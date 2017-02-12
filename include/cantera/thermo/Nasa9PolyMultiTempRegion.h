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

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_NASA9POLYMULTITEMPREGION_H
#define CT_NASA9POLYMULTITEMPREGION_H

#include "cantera/thermo/Nasa9Poly1.h"

namespace Cantera
{
//! The NASA 9 polynomial parameterization for a single species encompassing
//! multiple temperature regions.
/*!
 * The parameterization used in each temperature range is described in the
 * documentation for class Nasa9Poly1.
 *
 * These NASA representations may have multiple temperature regions through the
 * use of this Nasa9PolyMultiTempRegion object, which uses multiple copies of
 * the Nasa9Poly1 object to handle multiple temperature regions.
 *
 * @ingroup spthermo
 * @see Nasa9Poly1
 */
class Nasa9PolyMultiTempRegion : public SpeciesThermoInterpType
{
public:
    //! Constructor used in templated instantiations
    /*!
     * @param regionPts Vector of pointers to Nasa9Poly1 objects. These objects
     *     all refer to the temperature regions for the same species. The vector
     *     must be in increasing temperature region format.  Together they
     *     represent the reference temperature parameterization for a single
     *     species.
     *
     * Note, after the constructor, we will own the underlying Nasa9Poly1
     * objects and be responsible for owning them.
     */
    Nasa9PolyMultiTempRegion(std::vector<Nasa9Poly1*> &regionPts);

    virtual ~Nasa9PolyMultiTempRegion();

    virtual int reportType() const;

    virtual size_t temperaturePolySize() const { return 7; }
    virtual void updateTemperaturePoly(double T, double* T_poly) const;

    //! @copydoc Nasa9Poly1::updateProperties
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const;

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const;

    //! This utility function reports back the type of parameterization and all
    //! of the parameters for the species, index.
    /*!
     * All parameters are output variables
     *
     * @param n         Species index
     * @param type      Integer type of the standard type
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the parameters for
     *     the standard state. There are 1 + 11*nzones coefficients.
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

protected:
    //! Lower boundaries of each temperature regions
    vector_fp m_lowerTempBounds;

    //! Individual temperature region objects
    std::vector<std::unique_ptr<Nasa9Poly1>> m_regionPts;

    //! current region
    mutable int m_currRegion;
};

}
#endif
