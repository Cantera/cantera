/**
 * @file MultiSpeciesThermo.h
 *  Header for a general species thermodynamic property manager for a phase (see
 * \link Cantera::MultiSpeciesThermo MultiSpeciesThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTISPECIESTHERMO_H
#define CT_MULTISPECIESTHERMO_H

#include "SpeciesThermoInterpType.h"

namespace Cantera
{

//! A species thermodynamic property manager for a phase.
/*!
 * This is a general manager that can handle a wide variety of species
 * thermodynamic polynomials for individual species and compute their
 * nondimensional, reference-state thermodynamic properties (i.e. as a function
 * of temperature only).
 *
 * The ThermoPhase object relies on MultiSpeciesThermo to calculate the
 * thermodynamic properties of the reference state for all of the species in the
 * phase, for a range of temperatures. Note, the pressure dependence of the
 * species thermodynamic functions is not handled at this level. Species using
 * the same parameterization are grouped together in order to minimize the
 * operation count and achieve better efficiency.
 *
 * The most important member function for the MultiSpeciesThermo class is the
 * member function MultiSpeciesThermo::update(). The function calculates the
 * values of Cp/R, H/RT, and S/R for all of the species at once at the specified
 * temperature.
 *
 * Usually, all of the species in a phase are installed into a
 * MultiSpeciesThermo object. However, there is no requirement that a
 * MultiSpeciesThermo object handles all of the species in a phase. The member
 * function
 * \link MultiSpeciesThermo::install_STIT() install_STIT()\endlink
 * is called to install each species into the MultiSpeciesThermo object.
 *
 * @ingroup spthermo
 */
class MultiSpeciesThermo
{
public:
    //! Constructor
    MultiSpeciesThermo();

    // MultiSpeciesThermo objects are not copyable or assignable
    MultiSpeciesThermo(const MultiSpeciesThermo& b) = delete;
    MultiSpeciesThermo& operator=(const MultiSpeciesThermo& b) = delete;
    virtual ~MultiSpeciesThermo() {}

    //! Install a new species thermodynamic property parameterization for one
    //! species.
    /*!
     * @param index Index of the species being installed
     * @param stit Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(size_t index,
                              shared_ptr<SpeciesThermoInterpType> stit);

    //! Modify the species thermodynamic property parameterization for a species
    /*!
     * @param index Index of the species being installed
     * @param spec  Pointer to the SpeciesThermoInterpType object
     */
    virtual void modifySpecies(size_t index,
                               shared_ptr<SpeciesThermoInterpType> spec);

    //! Like update_one, but without applying offsets to the output pointers
    /*!
     * @param k       species index
     * @param T       Temperature (Kelvin)
     * @param cp_R    Dimensionless heat capacity
     * @param h_RT    Dimensionless enthalpy
     * @param s_R     Dimensionless entropy
     */
    virtual void update_single(size_t k, double T, double* cp_R,
                               double* h_RT, double* s_R) const;

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
                        doublereal* h_RT, doublereal* s_R) const;

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this method returns the minimum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum temperature for
     * species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k=npos) const;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this method returns the maximum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum temperature for
     * parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(size_t k=npos) const;

    //! The reference-state pressure for species k.
    /*!
     * Returns the reference state pressure in Pascals for species k. If k is
     * left out of the argument list, it returns the reference state pressure
     * for the first species.
     *
     * @param k Species Index
     */
    virtual doublereal refPressure(size_t k=npos) const;

    //! This utility function reports the type of parameterization used for the
    //! species with index number *index*.
    /*!
     * @param index  Species index
     */
    virtual int reportType(size_t index) const;

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
     */
    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const;

    //! Report the 298 K Heat of Formation of the standard state of one species
    //! (J kmol-1)
    /*!
     * The 298K Heat of Formation is defined as the enthalpy change to create
     * the standard state of the species from its constituent elements in their
     * standard states at 298 K and 1 bar.
     *
     * @param k    species index
     * @returns the current value of the Heat of Formation at 298K and 1 bar
     */
    virtual doublereal reportOneHf298(const size_t k) const;

    //! Modify the value of the 298 K Heat of Formation of the standard state of
    //! one species in the phase (J kmol-1)
    /*!
     * The 298K heat of formation is defined as the enthalpy change to create
     * the standard state of the species from its constituent elements in their
     * standard states at 298 K and 1 bar.
     *
     * @param  k           Index of the species
     * @param  Hf298New    Specify the new value of the Heat of Formation at
     *                     298K and 1 bar. units = J/kmol.
     */
    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New);

    //! Restore the original heat of formation of one or more species
    /*!
     *  Resets changes made by modifyOneHf298(). If the species index is not
     *  specified, the heats of formation for all species are restored.
     */
    virtual void resetHf298(const size_t k);

    //! Check if data for all species (0 through nSpecies-1) has been installed.
    bool ready(size_t nSpecies);

private:
    //! Provide the SpeciesThermoInterpType object
    /*!
     * @param k  species index
     * @return pointer to the SpeciesThermoInterpType object.
     */
    SpeciesThermoInterpType* provideSTIT(size_t k);
    const SpeciesThermoInterpType* provideSTIT(size_t k) const;

protected:
    //! Mark species *k* as having its thermodynamic data installed
    void markInstalled(size_t k);

    typedef std::pair<size_t, shared_ptr<SpeciesThermoInterpType> > index_STIT;
    typedef std::map<int, std::vector<index_STIT> > STIT_map;
    typedef std::map<int, vector_fp> tpoly_map;

    //! This is the main data structure, which contains the
    //! SpeciesThermoInterpType objects, sorted by the parameterization type.
    //! `m_sp[i]` is the vector of [species index, STIT] pairs which use
    //! parameterization `i`.
    STIT_map m_sp;

    //! Temperature polynomials for each thermo parameterization
    mutable tpoly_map m_tpoly;

    //! Map from species index to location within #m_sp, such that
    //! `m_sp[m_speciesLoc[k].first][m_speciesLoc[k].second]` is the
    //! SpeciesThermoInterpType object for species `k`.
    std::map<size_t, std::pair<int, size_t> > m_speciesLoc;

    //! Maximum value of the lowest temperature
    doublereal m_tlow_max;

    //! Minimum value of the highest temperature
    doublereal m_thigh_min;

    //! reference pressure (Pa)
    doublereal m_p0;

    //! indicates if data for species has been installed
    std::vector<bool> m_installed;
};

}

#endif
