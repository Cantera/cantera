/**
 *  @file SpeciesThermoInterpType.h
 * Pure Virtual Base class for individual species reference state
 * thermodynamic managers and text for the spthermo module
 * (see \ref spthermo and class \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType \endlink).
 */

// Copyright 2001  California Institute of Technology

#include "cantera/base/ct_defs.h"
#include "speciesThermoTypes.h"

#ifndef CT_SPECIESTHERMOINTERPTYPE_H
#define CT_SPECIESTHERMOINTERPTYPE_H

namespace Cantera
{

class PDSS;
class VPSSMgr;

/**
  * @defgroup spthermo Species Reference-State Thermodynamic Properties
  *
  *  The %ThermoPhase object relies on classes to calculate the thermodynamic
  *  properties of the reference state for all of the species in the phase.
  *  This group describes the types and functionality of the classes that
  *  calculate the reference state thermodynamic functions within %Cantera.
  *
  * To compute the thermodynamic properties of multicomponent
  * solutions, it is necessary to know something about the
  * thermodynamic properties of the individual species present in
  * the solution. Exactly what sort of species properties are
  * required depends on the thermodynamic model for the
  * solution. For a gaseous solution (i.e., a gas mixture), the
  * species properties required are usually ideal gas properties at
  * the mixture temperature and at a reference pressure (almost always at
  * 1 bar). For other types of solutions, however, it may
  * not be possible to isolate the species in a "pure" state. For
  * example, the thermodynamic properties of, say, Na+ and Cl- in
  * saltwater are not easily determined from data on the properties
  * of solid NaCl, or solid Na metal, or chlorine gas. In this
  * case, the solvation in water is fundamental to the identity of
  * the species, and some other reference state must be used. One
  * common convention for liquid solutions is to use thermodynamic
  * data for the solutes in the limit of infinite dilution within the
  * pure solvent; another convention is to reference all properties
  * to unit molality.
  *
  * In defining these standard states for species in a phase, we make
  * the following definition. A reference state is a standard state
  * of a species in a phase limited to one particular pressure, the reference
  * pressure. The reference state specifies the dependence of all
  * thermodynamic functions as a function of the temperature, in
  * between a minimum temperature and a maximum temperature. The
  * reference state also specifies the molar volume of the species
  * as a function of temperature. The molar volume is a thermodynamic
  * function.
  * A full standard state does the same thing as a reference state,
  * but specifies the thermodynamics functions at all pressures.
  *
  * Whatever the conventions used by a particular solution model,
  * means need to be provided to compute the species properties in
  * the reference state. Class SpeciesThermo is the base class
  * for a family of classes that compute properties of all
  * species in a phase in their reference states, for a range of temperatures.
  * Note, the pressure dependence of the species thermodynamic functions is not
  * handled by this particular species thermodynamic model. %SpeciesThermo
  * calculates the reference-state thermodynamic values of all species in a single
  * phase during each call.
  *
  * The class SpeciesThermoInterpType is a pure virtual base class for
  * calculation of thermodynamic functions for a single species
  * in its reference state.
  * The following classes inherit from %SpeciesThermoInterpType.
  *
  *   - NasaPoly1          in file NasaPoly1.h
  *      - This is a one zone model,  consisting of a 7
  *        coefficient Nasa Polynomial format.
  *      .
  *   - NasaPoly2          in file NasaPoly2.h
  *      - This is a two zone model, with each zone consisting of a 7
  *        coefficient Nasa Polynomial format.
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
  *      - This is a multizoned model. The chemical potential is given
  *        at a set number of temperatures. Between each temperature
  *        the heat capacity is treated as a constant.
  *      .
  *   - Nasa9Poly1          in file Nasa9Poly1.h
  *      - This is a one zone model,  consisting of the 9
  *        coefficient Nasa Polynomial format.
  *      .
  *   - Nasa9PolyMultiTempRegion       in file Nasa9PolyMultiTempRegion.h
  *      - This is a multiple zone model, consisting of the 9
  *        coefficient Nasa Polynomial format in each zone.
  *      .
  *   - STITbyPDSS          in file SpeciesThermoInterpType.h
  *      - This is an object that calculates reference state thermodynamic
  *        functions by relying on a pressure dependent
  *        standard state object (i.e., a PDSS object) to calculate
  *        the thermodynamic functions.
  *      .
  *
  *  The most important member function for the SpeciesThermoInterpType class
  *  is the member function
  *  \link SpeciesThermoInterpType::updatePropertiesTemp() updatePropertiesTemp()\endlink.
  *  The function calculates the values of Cp, H, and S for the specific
  *  species pertaining to this class. It takes as its arguments the
  *  base pointer for the vector of  Cp, H, and S values for all species
  *  in the phase. The offset for the species is known within the
  *  object.
  *
  *  A key concept for reference states is that there is a maximum and a minimum
  *  temperature beyond which the thermodynamic formulation isn't valid.
  *  Calls for temperatures outside this range will cause the
  *  object to throw a CanteraError.
  *
  * @ingroup thermoprops
  */

//!  Pure Virtual Base class for the thermodynamic manager for
//!  an individual species' reference state
/*!
 * This differs from the SpeciesThermo virtual
 * base class in the sense that this class is meant to handle only
 * one species. The speciesThermo class is meant to handle the
 * calculation of all the species (or a large subset) in a phase.
 *
 * One key feature is that the update routines use the same
 * form as the update routines in the speciesThermo class. They update
 * into a vector of cp_R, s_R, and H_R that spans all of the species in
 * a phase. Therefore, this class must carry along a species index into that
 * vector.
 *
 * These routine may be templated. A key requirement of the template is that
 * there is a constructor with the following form:
 *
 *  @code
 *   SpeciesThermoInterpType(int index, doublereal tlow, doublereal thigh,
 *                           doublereal pref,  const doublereal* coeffs)
 *  @endcode
 *
 * @ingroup spthermo
 */
class SpeciesThermoInterpType
{
public:

    //! Constructor
    SpeciesThermoInterpType();

    //! Constructor
    SpeciesThermoInterpType(size_t n, doublereal tlow,
                            doublereal thigh, doublereal pref) :
        m_lowT(tlow),
        m_highT(thigh),
        m_Pref(pref),
        m_index(n) {}

    //! Destructor
    virtual ~SpeciesThermoInterpType();

    //! duplicator
    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const = 0;

    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    virtual doublereal minTemp() const {
        return m_lowT;
    }

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    virtual doublereal maxTemp() const {
        return m_highT;
    }

    //! Returns the reference pressure (Pa)
    virtual doublereal refPressure() const {
        return m_Pref;
    }

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const = 0;

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
     * The form and length of the Temperature Polynomial may vary depending on the
     * parameterization.
     *
     * @param tempPoly  vector of temperature polynomials
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updateProperties(const doublereal* tempPoly,
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
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R,
                                      doublereal* h_RT,
                                      doublereal* s_R) const = 0;

    //!This utility function reports back the type of
    //! parameterization and all of the parameters for the
    //! species, index.
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
     * @deprecated
     */
    virtual void reportParameters(size_t& index, int& type,
                                  doublereal& minTemp, doublereal& maxTemp,
                                  doublereal& refPressure,
                                  doublereal* const coeffs) const = 0;

    //! Modify parameters for the standard state
    /*!
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     * @deprecated
     */
    virtual void modifyParameters(doublereal* coeffs) {}

#ifdef H298MODIFY_CAPABILITY

    //! Report the 298 K Heat of Formation of the standard state of one species (J kmol-1)
    /*!
     *   The 298K Heat of Formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param h298 If this is nonnull,  the current value of the Heat of Formation at 298K and 1 bar for
     *               species m_speciesIndex is returned in h298[m_speciesIndex].
     *   @return     Returns the current value of the Heat of Formation at 298K and 1 bar for
     *               species m_speciesIndex.
     */
    virtual doublereal reportHf298(doublereal* const h298 = 0) const;

    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar
     */
    virtual void modifyOneHf298(const int k, const doublereal Hf298New);

#endif

protected:
    //!  lowest valid temperature
    doublereal m_lowT;
    //! Highest valid temperature
    doublereal m_highT;
    //! Reference state pressure
    doublereal m_Pref;
    //! species index
    size_t m_index;
};

//!  Class for the thermodynamic manager for an individual species' reference state
//!  which uses the PDSS base class to satisfy the requests.
/*!
 *  This class is a pass-through class for handling thermodynamics calls
 *  for reference state thermo to an pressure dependent standard state (PDSS)
 *  class. For some situations, it makes no sense to have a reference state
 *  at all. One example of this is the real water standard state.
 *
 *  What this class does is just to pass through the calls for thermo at (T , p0)
 *  to the PDSS class, which evaluates the calls at (T, p0).
 *
 * @ingroup spthermo
 */
class STITbyPDSS : public SpeciesThermoInterpType
{
public:
    //! Constructor
    STITbyPDSS();

    //! Main Constructor
    /*!
     *  @param speciesIndex species index for this object. Note, this must
     *         agree with what was internally set before.
     *
     *  @param vpssmgr_ptr  Pointer to the Variable pressure standard state manager
     *                      that owns the PDSS object that will handle calls for this object
     *
     *  @param PDSS_ptr     Pointer to the PDSS object that handles calls for this object
     */
    STITbyPDSS(size_t speciesIndex, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr);

    //! copy constructor
    /*!
     *  @param b Object to be copied
     */
    STITbyPDSS(const STITbyPDSS& b);

    virtual SpeciesThermoInterpType* duplMyselfAsSpeciesThermoInterpType() const;

    //! Initialize and/or Reinitialize all the pointers for this object
    /*!
     *  This routine is needed because the STITbyPDSS object doesn't own the
     *  underlying objects. Therefore, shallow copies during duplication operations
     *  may fail.
     *
     *  @param speciesIndex species index for this object. Note, this must
     *         agree with what was internally set before.
     *
     *  @param vpssmgr_ptr  Pointer to the Variable pressure standard state manager
     *                      that owns the PDSS object that will handle calls for this object
     *
     *  @param PDSS_ptr     Pointer to the PDSS object that handles calls for this object
     */
    void initAllPtrs(size_t speciesIndex, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr);

    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    virtual doublereal minTemp() const;

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    virtual doublereal maxTemp() const;

    //! Returns the reference pressure (Pa)
    virtual doublereal refPressure() const;

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const;

    virtual void updateProperties(const doublereal* tempPoly,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const;

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R,
                                      doublereal* h_RT,
                                      doublereal* s_R) const;

    //! @deprecated
    virtual void reportParameters(size_t& index, int& type,
                                  doublereal& minTemp, doublereal& maxTemp,
                                  doublereal& refPressure,
                                  doublereal* const coeffs) const;

    //! @deprecated
    virtual void modifyParameters(doublereal* coeffs);

private:
    //! Pointer to the Variable pressure standard state manager
    //! that owns the PDSS object that will handle calls for this object
    VPSSMgr* m_vpssmgr_ptr;

    //! Pointer to the PDSS object that handles calls for this object
    /*!
     * This object is not owned by the current one.
     */
    PDSS* m_PDSS_ptr;
};

}

#endif
