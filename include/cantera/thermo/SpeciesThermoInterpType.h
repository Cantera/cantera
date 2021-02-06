/**
 *  @file SpeciesThermoInterpType.h
 *
 * Pure Virtual Base class for individual species reference state thermodynamic
 * managers and text for the spthermo module (see \ref spthermo and class \link
 * Cantera::SpeciesThermoInterpType SpeciesThermoInterpType \endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SPECIESTHERMOINTERPTYPE_H
#define CT_SPECIESTHERMOINTERPTYPE_H

#include "cantera/base/ct_defs.h"
#include "speciesThermoTypes.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

class PDSS;

/**
  * @defgroup spthermo Species Reference-State Thermodynamic Properties
  *
  * To compute the thermodynamic properties of multicomponent solutions, it is
  * necessary to know something about the thermodynamic properties of the
  * individual species present in the solution. Exactly what sort of species
  * properties are required depends on the thermodynamic model for the solution.
  * For a gaseous solution (i.e., a gas mixture), the species properties
  * required are usually ideal gas properties at the mixture temperature and at
  * a reference pressure (almost always at 1 bar). For other types of solutions,
  * however, it may not be possible to isolate the species in a "pure" state.
  * For example, the thermodynamic properties of, say, Na+ and Cl- in saltwater
  * are not easily determined from data on the properties of solid NaCl, or
  * solid Na metal, or chlorine gas. In this case, the solvation in water is
  * fundamental to the identity of the species, and some other reference state
  * must be used. One common convention for liquid solutions is to use
  * thermodynamic data for the solutes in the limit of infinite dilution within
  * the pure solvent; another convention is to reference all properties to unit
  * molality.
  *
  * In defining these standard states for species in a phase, we make the
  * following definition. A reference state is a standard state of a species in
  * a phase limited to one particular pressure, the reference pressure. The
  * reference state specifies the dependence of all thermodynamic functions as a
  * function of the temperature, in between a minimum temperature and a maximum
  * temperature. The reference state also specifies the molar volume of the
  * species as a function of temperature. The molar volume is a thermodynamic
  * function. A full standard state does the same thing as a reference state,
  * but specifies the thermodynamics functions at all pressures.
  *
  * The class SpeciesThermoInterpType is an abstract base class for calculation
  * of thermodynamic functions for a single species in its reference state. The
  * following classes inherit from SpeciesThermoInterpType.
  *
  *   - NasaPoly1          in file NasaPoly1.h
  *      - This is a one zone model, consisting of a 7
  *        coefficient NASA Polynomial format.
  *      .
  *   - NasaPoly2          in file NasaPoly2.h
  *      - This is a two zone model, with each zone consisting of a 7
  *        coefficient NASA Polynomial format.
  *      .
  *   - ShomatePoly        in file ShomatePoly.h
  *      - This is a one zone model, consisting of a 7
  *        coefficient Shomate Polynomial format.
  *      .
  *   - ShomatePoly2       in file ShomatePoly.h
  *      - This is a two zone model, with each zone consisting of a 7
  *        coefficient Shomate Polynomial format.
  *      .
  *   - ConstCpPoly        in file ConstCpPoly.h
  *      - This is a one-zone constant heat capacity model.
  *      .
  *   - Mu0Poly            in file Mu0Poly.h
  *      - This is a multi-zone model. The chemical potential is given
  *        at a set number of temperatures. Between each temperature
  *        the heat capacity is treated as a constant.
  *      .
  *   - Nasa9Poly1          in file Nasa9Poly1.h
  *      - This is a one zone model, consisting of the 9
  *        coefficient NASA Polynomial format.
  *      .
  *   - Nasa9PolyMultiTempRegion       in file Nasa9PolyMultiTempRegion.h
  *      - This is a multiple zone model, consisting of the 9
  *        coefficient NASA Polynomial format in each zone.
  *      .
  * The most important member function for the SpeciesThermoInterpType class is
  * the member function SpeciesThermoInterpType::updatePropertiesTemp(). The
  * function calculates the values of Cp, H, and S for the specific species
  * pertaining to this class.
  *
  * A key concept for reference states is that there is a maximum and a minimum
  * temperature beyond which the thermodynamic formulation isn't valid. Calls
  * for temperatures outside this range will cause the object to throw a
  * CanteraError.
  *
  * @ingroup thermoprops
  */

//! Abstract Base class for the thermodynamic manager for an individual
//! species' reference state
/*!
 * One key feature is that the update routines use the same form as the update
 * routines in the MultiSpeciesThermo class. They update values of cp_R,
 * s_R, and H_R.
 *
 * @ingroup spthermo
 */
class SpeciesThermoInterpType
{
public:
    SpeciesThermoInterpType();

    SpeciesThermoInterpType(double tlow, double thigh, double pref);

    // SpeciesThermoInterpType objects are not copyable or assignable
    SpeciesThermoInterpType(const SpeciesThermoInterpType& b) = delete;
    SpeciesThermoInterpType& operator=(const SpeciesThermoInterpType& b) = delete;

    virtual ~SpeciesThermoInterpType() {}

    //! Returns the minimum temperature that the thermo parameterization is
    //! valid
    virtual doublereal minTemp() const {
        return m_lowT;
    }

    //! Set the minimum temperature at which the thermo parameterization is valid
    virtual void setMinTemp(double Tmin) {
        m_lowT = Tmin;
    }

    //! Returns the maximum temperature that the thermo parameterization is
    //! valid
    virtual doublereal maxTemp() const {
        return m_highT;
    }

    //! Set the maximum temperature at which the thermo parameterization is valid
    virtual void setMaxTemp(double Tmax) {
        m_highT = Tmax;
    }

    //! Returns the reference pressure (Pa)
    virtual doublereal refPressure() const {
        return m_Pref;
    }

    //! Set the reference pressure [Pa]
    virtual void setRefPressure(double Pref) {
        m_Pref = Pref;
    }

    //! Check for problems with the parameterization, and generate warnings or
    //! throw and exception if any are found.
    virtual void validate(const std::string& name) {}

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const { return 0; };

    //! Number of terms in the temperature polynomial for this parameterization
    virtual size_t temperaturePolySize() const { return 1; }

    //! Given the temperature *T*, compute the terms of the temperature
    //! polynomial *T_poly*.
    virtual void updateTemperaturePoly(double T, double* T_poly) const {
        T_poly[0] = T;
    }

    //! Update the properties for this species, given a temperature polynomial
    /*!
     * This method is called with a pointer to an array containing the functions
     * of temperature needed by this parameterization, and three pointers to
     * arrays where the computed property values should be written. This method
     * updates only one value in each array.
     *
     * The form and length of the Temperature Polynomial may vary depending on
     * the parameterization.
     *
     * @param tt      vector of evaluated temperature functions
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const;

    //! Compute the reference-state property of one species
    /*!
     * Given temperature T in K, this method updates the values of the non-
     * dimensional heat capacity at constant pressure, enthalpy, and entropy, at
     * the reference pressure, of the species.
     *
     * @param temp    Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R,
                                      doublereal* h_RT,
                                      doublereal* s_R) const;

    //! This utility function returns the number of coefficients
    //! for a given type of species parameterization
    virtual size_t nCoeffs() const;

    //! This utility function returns the type of parameterization and all
    //! of the parameters for the species.
    /*!
     * All parameters are output variables
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the
     *                  parameters for the standard state.
     */
    virtual void reportParameters(size_t& index, int& type,
                                  doublereal& minTemp, doublereal& maxTemp,
                                  doublereal& refPressure,
                                  doublereal* const coeffs) const;

    //! Report the 298 K Heat of Formation of the standard state of one species
    //! (J kmol-1)
    /*!
     * The 298K Heat of Formation is defined as the enthalpy change to create
     * the standard state of the species from its constituent elements in their
     * standard states at 298 K and 1 bar.
     *
     * @param h298 If this is nonnull, the current value of the Heat of
     *             Formation at 298K and 1 bar for species m_speciesIndex is
     *             returned in h298[m_speciesIndex].
     * @return the current value of the Heat of Formation at 298K and 1 bar for
     *               species m_speciesIndex.
     */
    virtual doublereal reportHf298(doublereal* const h298 = 0) const;

    //! Modify the value of the 298 K Heat of Formation of one species in the
    //! phase (J kmol-1)
    /*!
     * The 298K heat of formation is defined as the enthalpy change to create
     * the standard state of the species from its constituent elements in their
     * standard states at 298 K and 1 bar.
     *
     * @param  k           Species k
     * @param  Hf298New    Specify the new value of the Heat of Formation at
     *                     298K and 1 bar
     */
    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New);

    //! Restore the original heat of formation for this species
    /*!
     *  Resets changes made by modifyOneHf298().
     */
    virtual void resetHf298() {
        throw NotImplementedError("SpeciesThermoInterpType::resetHf298");
    }

protected:
    //!  lowest valid temperature
    doublereal m_lowT;
    //! Highest valid temperature
    doublereal m_highT;
    //! Reference state pressure
    doublereal m_Pref;
};

}

#endif
