/**
 *  @file PDSS_SSVol.h
 *    Declarations for the class PDSS_SSVol (pressure dependent standard state)
 *    which handles calculations for a single species with an expression for the standard state molar volume in a phase
 *    given by an enumerated data type
 *    (see class \ref pdssthermo and \link Cantera::PDSS_SSVol PDSS_SSVol\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_SSVOL_H
#define CT_PDSS_SSVOL_H

#include "PDSS.h"

namespace Cantera
{
//! Class for pressure dependent standard states that uses a standard state
//! volume model of some sort.
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * Class PDSS_SSVol is an implementation class that compute the properties of a
 * single species in a phase at its standard states, for a range of temperatures
 * and pressures. This particular class assumes that the calculation of the
 * thermodynamics functions can be separated into a temperature polynomial
 * representation for thermo functions that can be handled by a SimpleThermo
 * object and a separate calculation for the standard state volume. The Models
 * include a cubic polynomial in temperature for either the standard state
 * volume or the standard state density. The manager uses a SimpleThermo object
 * to handle the calculation of the reference state. This object then adds the
 * pressure dependencies and the volume terms to these thermo functions to
 * complete the representation.
 *
 * The class includes the following models for the representation of the
 * standard state volume:
 *
 * - Temperature polynomial for the standard state volume
 *   - This standard state model is invoked with the keyword "temperature_polynomial".
 *     The standard state volume is considered a function of temperature only.
 *     \f[
 *       V^o_k(T,P) = a_0 + a_1 T + a_2 T^2 + a_3 T^3
 *     \f]
 *
 * - Temperature polynomial for the standard state density
 *   - This standard state model is invoked with the keyword "density_temperature_polynomial".
 *     The standard state density, which is the inverse of the volume,
 *     is considered a function of temperature only.
 *    \f[
 *       {\rho}^o_k(T,P) = \frac{M_k}{V^o_k(T,P)} = a_0 + a_1 T + a_2 T^2 + a_3 T^3
 *    \f]
 *
 * ## Specification of Species Standard State Properties
 *
 * The standard molar Gibbs free energy for species *k* is determined from
 * the enthalpy and entropy expressions
 *
 *       \f[
 *            G^o_k(T,P) = H^o_k(T,P) - S^o_k(T,P)
 *       \f]
 *
 * The enthalpy is calculated mostly from the MultiSpeciesThermo object's enthalpy
 * evaluator. The dependence on pressure originates from the Maxwell relation
 *
 *       \f[
 *            {\left(\frac{dH^o_k}{dP}\right)}_T = T  {\left(\frac{dS^o_k}{dP}\right)}_T + V^o_k
 *       \f]
 * which is equal to
 *
 *       \f[
 *            {\left(\frac{dH^o_k}{dP}\right)}_T =  V^o_k -  T  {\left(\frac{dV^o_k}{dT}\right)}_P
 *       \f]
 *
 * The entropy is calculated mostly from the MultiSpeciesThermo objects entropy
 * evaluator. The dependence on pressure originates from the Maxwell relation:
 *
 *       \f[
 *              {\left(\frac{dS^o_k}{dP}\right)}_T =  - {\left(\frac{dV^o_k}{dT}\right)}_P
 *       \f]
 *
 * The standard state constant-pressure heat capacity expression is obtained
 * from taking the temperature derivative of the Maxwell relation involving the
 * enthalpy given above to yield an expression for the pressure dependence of
 * the heat capacity.
 *
 *       \f[
 *            {\left(\frac{d{C}^o_{p,k}}{dP}\right)}_T =  - T  {\left(\frac{{d}^2{V}^o_k}{{dT}^2}\right)}_T
 *       \f]
 *
 * The standard molar Internal Energy for species *k* is determined from the
 * following relation.
 *
 *       \f[
 *            U^o_k(T,P) = H^o_k(T,P) - p V^o_k
 *       \f]
 *
 * ## XML Example
 *
 * An example of the specification of a standard state for the LiCl molten salt
 * which has a temperature dependent standard state volume.
 *
 * @code
 * <speciesData id="species_MoltenSalt">
 * <species name="LiCl(L)">
 *    <atomArray> Li:1 Cl:1 </atomArray>
 *    <standardState  model="density_temperature_polynomial">
 *       <densityTemperaturePolynomial units="gm/cm3" >
 *          1.98715, -5.890906E-4, 0.0, 0.0
 *       </densityTemperaturePolynomial>
 *    </standardState>
 *    <thermo>
 *      <Shomate Pref="1 bar" Tmax="2000.0" Tmin="700.0">
 *        <floatArray size="7">
 *          73.18025, -9.047232, -0.316390,
 *          0.079587, 0.013594, -417.1314,
 *          157.6711
 *        </floatArray>
 *      </Shomate>
 *    </thermo>
 *  </species>
 *  </speciesData>
 * @endcode
 *
 * @ingroup pdssthermo
 */
class PDSS_SSVol : public PDSS_Nondimensional
{
public:
    //! Default Constructor
    PDSS_SSVol();

    //! @name  Molar Thermodynamic Properties of the Species Standard State in the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal intEnergy_mole() const;
    virtual doublereal cv_mole() const;

    //! @}

    //! @name Mechanical Equation of State Properties
    //! @{

    virtual void setPressure(doublereal pres);
    virtual void setTemperature(doublereal temp);
    virtual void setState_TP(doublereal temp, doublereal pres);
    //! @deprecated To be removed after Cantera 2.5.
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //! @name Miscellaneous properties of the standard state
    //! @{

    virtual doublereal satPressure(doublereal t);

    //! @}
    //! @name Initialization of the Object
    //! @{

    virtual void initThermo();

    //! Set polynomial coefficients for the standard state molar volume as a
    //! function of temperature. Cubic polynomial (4 coefficients). Leading
    //! coefficient is the constant (temperature-independent) term [m^3/kmol].
    void setVolumePolynomial(double* coeffs);

    //! Set polynomial coefficients for the standard state density as a function
    //! of temperature. Cubic polynomial (4 coefficients). Leading coefficient
    //! is the constant (temperature-independent) term [kg/m^3].
    void setDensityPolynomial(double* coeffs);

    virtual void setParametersFromXML(const XML_Node& speciesNode);
    //@}

private:
    //! Does the internal calculation of the volume
    void calcMolarVolume();

    //! Types of general formulations for the specification of the standard
    //! state volume
    enum class SSVolume_Model {
        //! This approximation is for a species with a cubic polynomial in
        //! temperature
        /*!
         *       V^ss = a_0 + a_1 T + a_2 T^2 + a_3 T^3
         */
        tpoly,
        //! This approximation is for a species where the density is expressed
        //! as a cubic polynomial in temperature
        /*!
         *       V^ss = M / (a_0 + a_1 T + a_2 T^2 + a_3 T^3)
         */
        density_tpoly
    };

    //! Enumerated data type describing the type of volume model
    //! used to calculate the standard state volume of the species
    SSVolume_Model volumeModel_;

    //! coefficients for the temperature representation
    vector_fp TCoeff_;

    //! Derivative of the volume wrt temperature
    mutable doublereal dVdT_;

    //! 2nd derivative of the volume wrt temperature
    mutable doublereal d2VdT2_;
};

}

#endif
