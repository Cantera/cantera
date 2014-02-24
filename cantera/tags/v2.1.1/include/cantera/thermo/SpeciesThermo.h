/**
 *  @file SpeciesThermo.h
 *  Virtual base class for the calculation of multiple-species thermodynamic
 *  reference-state property managers and text for the mgrsrefcalc module (see \ref mgrsrefcalc
 *  and class \link Cantera::SpeciesThermo SpeciesThermo\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_SPECIESTHERMO_H
#define CT_SPECIESTHERMO_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{
class SpeciesThermoInterpType;

/**
 * @defgroup mgrsrefcalc Managers for Calculating Reference-State Thermodynamics
 *
 *  The ThermoPhase object relies on a set of manager classes to calculate
 *  the thermodynamic properties of the reference state for all
 *  of the species in the phase. This may be a computationally
 *  significant cost, so efficiency is important.
 *  This group describes how this is done efficiently within Cantera.
 *
 * To compute the thermodynamic properties of multicomponent
 * solutions, it is necessary to know something about the
 * thermodynamic properties of the individual species present in
 * the solution. Exactly what sort of species properties are
 * required depends on the thermodynamic model for the
 * solution. For a gaseous solution (i.e., a gas mixture), the
 * species properties required are usually ideal gas properties at
 * the mixture temperature and at a reference pressure (almost always at
 * 1 bar).
 *
 * In defining these standard states for species in a phase, we make
 * the following definition. A reference state is a standard state
 * of a species in a phase limited to one particular pressure, the reference
 * pressure. The reference state specifies the dependence of all
 * thermodynamic functions as a function of the temperature, in
 * between a minimum temperature and a maximum temperature. The
 * reference state also specifies the molar volume of the species
 * as a function of temperature. The molar volume is a thermodynamic
 * function. By contrast, a full standard state does the same thing
 * as a reference state, but specifies the thermodynamics functions
 * at all pressures.
 *
 *  Whatever the conventions used by a particular solution model,
 *  means need to be provided to compute the species properties in
 *  the reference state. Class SpeciesThermo is the base class
 *  for a family of classes that compute properties of all
 *  species in a phase in their reference states, for a range of temperatures.
 *  Note, the pressure dependence of the species thermodynamic functions is not
 *  handled by this particular species thermodynamic model. %SpeciesThermo
 *  calculates the reference-state thermodynamic values of all species in a single
 *  phase during each call. The vector nature of the operation leads to
 *  a lower operation count and better efficiency, especially if the
 *  individual reference state classes are known to the reference-state
 *  manager class so that common operations may be grouped together.
 *
 *  The most important member function for the %SpeciesThermo class
 *  is the member function \link SpeciesThermo::update() update()\endlink.
 *  The function calculates the values of Cp, H, and S for all of the
 *  species at once at the specified temperature.
 *
 *  Usually, all of the species in a phase are installed into a %SpeciesThermo
 *  class. However, there is no requirement that a %SpeciesThermo
 *  object handles all of the species in a phase. There are
 *  two member functions that are called to install each species into
 *  the %SpeciesThermo.
 *  One routine is called \link SpeciesThermo::install() install()\endlink.
 *  It is called with the index of the species in the phase,
 *  an integer type delineating
 *  the SpeciesThermoInterpType object, and a listing of the
 *  parameters for that parameterization. A factory routine is called based
 *  on the integer type.  The other routine is called
 *  \link SpeciesThermo::install_STIT() install_STIT()\endlink.
 *  It accepts as an argument a pointer to an already formed
 *  SpeciesThermoInterpType object.
 *
 *  The following classes inherit from %SpeciesThermo. Each of these classes
 *  handle multiple species, usually all of the species in a phase. However,
 *  there is no requirement that a %SpeciesThermo object handles all of the
 *  species in a phase.
 *
 *   - NasaThermo          in file NasaThermo.h
 *      - This is a two zone model, with each zone consisting of a 7
 *        coefficient Nasa Polynomial format.
 *   - ShomateThermo       in file ShomateThermo.h
 *      - This is a two zone model, with each zone consisting of a 7
 *        coefficient Shomate Polynomial format.
 *   - SimpleThermo        in file SimpleThermo.h
 *      - This is a one-zone constant heat capacity model.
 *   - GeneralSpeciesThermo in file GeneralSpeciesThermo.h
 *      - This is a general model. Each species is handled separately
 *        via a vector over SpeciesThermoInterpType classes.
 *   - SpeciesThermoDuo      in file SpeciesThermoMgr.h
 *      - This is a combination of two SpeciesThermo types.
 *
 * The class SpeciesThermoInterpType is a pure virtual base class for
 * calculation of thermodynamic functions for a single species
 * in its reference state.
 * The following classes inherit from %SpeciesThermoInterpType.
 *
 *   - NasaPoly1          in file NasaPoly1.h
 *      - This is a one zone model,  consisting of a 7
 *        coefficient Nasa Polynomial format.
 *   - NasaPoly2          in file NasaPoly2.h
 *      - This is a two zone model, with each zone consisting of a 7
 *        coefficient Nasa Polynomial format.
 *   - ShomatePoly        in file ShomatePoly.h
 *      - This is a one zone model, consisting of a 7
 *        coefficient Shomate Polynomial format.
 *   - ShomatePoly2       in file ShomatePoly.h
 *      - This is a two zone model, with each zone consisting of a 7
 *        coefficient Shomate Polynomial format.
 *   - ConstCpPoly        in file ConstCpPoly.h
 *      - This is a one-zone constant heat capacity model.
 *   - Mu0Poly            in file Mu0Poly.h
 *      - This is a multizoned model. The chemical potential is given
 *        at a set number of temperatures. Between each temperature
 *        the heat capacity is treated as a constant.
 *   - Nasa9Poly1          in file Nasa9Poly1.h
 *      - This is a one zone model,  consisting of the 9
 *        coefficient Nasa Polynomial format.
 *   - Nasa9PolyMultiTempRegion       in file Nasa9PolyMultiTempRegion.h
 *      - This is a multiple zone model, consisting of the 9
 *        coefficient Nasa Polynomial format in each zone.
 *
 * In particular the NasaThermo %SpeciesThermo-derived model has been
 * optimized for execution speed. It's the main-stay of gas phase computations
 * involving large numbers of species in a phase. It combines the calculation
 * of each species, which individually have NasaPoly2 representations, to
 * minimize the computational time.
 *
 * The GeneralSpeciesThermo %SpeciesThermo object is completely general. It
 * does not try to coordinate the individual species calculations at all and
 * therefore is the slowest but most general implementation.
 *
 * @ingroup thermoprops
 */
//@{

//! Pure Virtual base class for the species thermo manager classes.
/*!
 *  This class defines the interface which all subclasses must implement.
 *
 * Class SpeciesThermo is the base class for a family of classes that compute
 * properties of a set of species in their reference state at a range of
 * temperatures. Note, the pressure dependence of the reference state is not
 * handled by this particular species standard state model.
 */
class SpeciesThermo
{
public:

    //! Constructor
    SpeciesThermo() : m_allow_discontinuities(false) {}

    //! Destructor
    virtual ~SpeciesThermo() {}

    //! Copy Constructor for the %SpeciesThermo object.
    /*!
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    SpeciesThermo(const SpeciesThermo& right) {}

    //! Assignment operator for the %SpeciesThermo object
    /*!
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    SpeciesThermo& operator=(const SpeciesThermo& right) {
        return *this;
    }

    //! Duplication routine for objects derived from SpeciesThermo
    /*!
     *  This function can be used to duplicate objects derived from
     *  SpeciesThermo even if the application only has a pointer to
     *  SpeciesThermo to work with.
     */
    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const = 0;

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @see speciesThermoTypes.h
     *
     * @param name      Name of the species
     * @param index     The 'update' method will update the property
     *                  values for this species
     *                  at position i index in the property arrays.
     * @param type      int flag specifying the type of parameterization to be
     *                 installed.
     * @param c        vector of coefficients for the parameterization.
     *                 This vector is simply passed through to the
     *                 parameterization constructor.
     * @param minTemp  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure standard-state pressure for this
     *                    parameterization.
     */
    virtual void install(const std::string& name, size_t index, int type,
                         const doublereal* c,
                         doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure)=0;

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType* stit_ptr) = 0;

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of the non-
     * dimensional heat capacity at constant pressure, enthalpy, and entropy,
     * at the reference pressure, Pref of each of the standard states.
     *
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void update(doublereal T, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const=0;

    //! Like update(), but only updates the single species k.
    /*!
     *  The default treatment is to just call update() which means that
     *  potentially the operation takes a m_kk*m_kk hit.
     *
     * @param k       species index
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void update_one(size_t k, doublereal T,
                            doublereal* cp_R,
                            doublereal* h_RT,
                            doublereal* s_R) const {
        update(T, cp_R, h_RT, s_R);
    }

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this method returns the minimum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum temperature for
     * species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k=npos) const =0;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this method returns the maximum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum temperature for
     * parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(size_t k=npos) const =0;

    //! The reference-state pressure for species k.
    /*!
     * Returns the reference state pressure in Pascals for species k. If k is
     * left out of the argument list, it returns the reference state pressure
     * for the first species. Note that some SpeciesThermo implementations,
     * such as those for ideal gases, require that all species in the same
     * phase have the same reference state pressures.
     *
     * @param k Species Index
     */
    virtual doublereal refPressure(size_t k=npos) const =0;

    //! This utility function reports the type of parameterization
    //! used for the species with index number *index*.
    /*!
     * @param index  Species index
     */
    virtual int reportType(size_t index=npos) const = 0;

    //! This utility function reports back the type of parameterization and
    //! all of the parameters for the species with index number *index*.
    /*!
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     * @deprecated
     */
    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const =0;

#ifdef H298MODIFY_CAPABILITY
    //! Report the 298 K Heat of Formation of the standard state of one species (J kmol-1)
    /*!
     *   The 298K Heat of Formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param k    species index
     *   @return     Returns the current value of the Heat of Formation at 298K and 1 bar
     */
    virtual doublereal reportOneHf298(int k) const = 0;

    //!  Modify the value of the 298 K Heat of Formation of the standard state of
    //!  one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Index of the species
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar.
     *                       units = J/kmol.
     */
    virtual void modifyOneHf298(const int k, const doublereal Hf298New) = 0;
#endif

    bool m_allow_discontinuities;
};
//@}
}

#endif
