/**
 * @file NasaThermo.h
 *   Header for the 2 regime 7 coefficient NASA thermodynamic
 *   polynomials for multiple species in a phase, derived from the
 *   \link Cantera::SpeciesThermo SpeciesThermo\endlink base class (see \ref mgrsrefcalc and
 *   \link Cantera::NasaThermo NasaThermo\endlink).
 */
// Copyright 2003 California Institute of Technology

#ifndef CT_NASATHERMO_H
#define CT_NASATHERMO_H

#include "cantera/thermo/SpeciesThermoMgr.h"
#include "cantera/thermo/NasaPoly1.h"

namespace Cantera
{
/**
 * A species thermodynamic property manager for the NASA
 * polynomial parameterization with two temperature ranges.
 *
 * This class is designed to efficiently evaluate the properties
 * of a large number of species with the NASA parameterization.
 *
 * The original NASA polynomial parameterization expressed the
 * heat capacity as a fourth-order polynomial in temperature, with
 * separate coefficients for each of two temperature ranges. (The
 * newer NASA format adds coefficients for 1/T and 1/T^2, and
 * allows multiple temperature ranges.) This class is designed for
 * use with the original parameterization, which is used, for
 * example, by the Chemkin software package.
 *
 * In many cases, the midpoint temperature is the same for many
 * species.  To take advantage of this, class NasaThermo groups
 * species with a common midpoint temperature, so that checking
 * which range the desired temperature is in need be done only
 * once for each group.
 *
 * @note There is a special CTML element for entering the
 * coefficients of this parameterization.
 * @see importCTML
 *
 * @ingroup mgrsrefcalc
 * @deprecated To be removed after Cantera 2.2. Use GeneralSpeciesThermo instead.
 */
class NasaThermo : public SpeciesThermo
{
public:
    NasaThermo();

    NasaThermo(const NasaThermo& right);

    NasaThermo& operator=(const NasaThermo& right);

    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const {
        NasaThermo* nt = new NasaThermo(*this);
        return (SpeciesThermo*) nt;
    }

    //! install a new species thermodynamic property
    //!  parameterization for one species.
    /*!
     * @param name     Name of the species
     * @param index    The 'update' method will update the property values for
     *                 this species at position i index in the property
     *                 arrays.
     * @param type     int flag specifying the type of parameterization to be
     *                 installed.
     * @param c        vector of coefficients for the parameterization.
     * - c[0]          midpoint temperature
     * - c[1] - c[7]   coefficients for low T range
     * - c[8] - c[14]  coefficients for high T range
     * @param min_temp minimum temperature for which this parameterization
     *                 is valid.
     * @param max_temp maximum temperature for which this parameterization
     *                 is valid.
     * @param ref_pressure standard-state pressure for this parameterization.
     * @see speciesThermoTypes.h
     */
    virtual void install(const std::string& name, size_t index, int type,
                         const doublereal* c,
                         doublereal min_temp, doublereal max_temp,
                         doublereal ref_pressure);

    virtual void install_STIT(size_t index, shared_ptr<SpeciesThermoInterpType> stit_ptr) {
        throw CanteraError("install_STIT", "not implemented");
    }

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void update_one(size_t k, doublereal t, doublereal* cp_R,
                            doublereal* h_RT, doublereal* s_R) const;

    virtual void update(doublereal t, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const;

    virtual doublereal minTemp(size_t k=npos) const {
        if (k == npos) {
            return m_tlow_max;
        } else {
            return m_tlow[k];
        }
    }

    virtual doublereal maxTemp(size_t k=npos) const {
        if (k == npos) {
            return m_thigh_min;
        } else {
            return m_thigh[k];
        }
    }

    virtual doublereal refPressure(size_t k=npos) const {
        return m_p0;
    }

    virtual int reportType(size_t index) const {
        return NASA;
    }

    /*!
     * This utility function reports back the type of
     * parameterization and all of the parameters for the
     * species, index.
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * For the NASA object, there are 15 coefficients.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const;

    virtual doublereal reportOneHf298(const size_t k) const;
    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New);

    //! Initialized to the type of parameterization
    /*!
     * Note, this value is used in some template functions
     */
    const int ID;

protected:
    //! Vector of vector of NasaPoly1's for the high temp region.
    /*!
     * This is the high temp region representation.
     * The first Length is equal to the number of groups.
     * The second vector is equal to the number of species
     * in that particular group.
     */
    std::vector<std::vector<NasaPoly1> > m_high;

    //! Vector of vector of NasaPoly1's for the low temp region.
    /*!
     * This is the low temp region representation.
     * The first Length is equal to the number of groups.
     * The second vector is equal to the number of species
     * in that particular group.
     */
    std::vector<std::vector<NasaPoly1> > m_low;

    //! Map between the midpoint temperature, as an int, to the group number
    /*!
     * Length is equal to the number of groups. Only used in the setup.
     */
    std::map<int, int>                 m_index;

    //! Vector of log temperature limits
    /*!
     * Length is equal to the number of groups.
     */
    vector_fp                          m_tmid;

    //! Maximum value of the low temperature limit
    doublereal                         m_tlow_max;

    //! Minimum value of the high temperature limit
    doublereal                         m_thigh_min;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp                          m_tlow;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp                          m_thigh;

    //! Reference pressure (Pa)
    /*!
     * all species must have the same reference pressure.
     */
    doublereal                         m_p0;

    //! number of groups
    int                                m_ngroups;

    //! Vector of temperature polynomials
    mutable vector_fp                  m_t;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the group index, where the temperature polynomials
     * for that species are stored. group indices start at 1,
     * so a decrement is always performed to access vectors.
     */
    std::map<size_t, size_t> m_group_map;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the position index within the group, where the
     * temperature polynomials for that species are stored.
     */
    std::map<size_t, size_t> m_posInGroup_map;

    //! Species name as a function of the species index
    std::map<size_t, std::string> m_name;

protected:
    //! Compute the nondimensional heat capacity using the given NASA polynomial
    /*!
     * @param t temperature
     * @param c coefficient array
     */
     doublereal cp_R(double t, const doublereal* c);

    //! Compute the nondimensional enthalpy using the given NASA polynomial
    /*!
     * @param t temperature
     * @param c coefficient array
     */
    doublereal enthalpy_RT(double t, const doublereal* c);

    //! Compute the nondimensional entropy using the given NASA polynomial
    /*!
     * @param t temperature
     * @param c coefficient array
     */
    doublereal entropy_R(double t, const doublereal* c);

    //! Adjust polynomials to be continuous at the midpoint temperature.
    /*!
     * Check to see if the provided coefficients are nearly continuous. Adjust
     * the values to get more precise continuity to avoid convergence
     * issues with algorithms that expect these quantities to be continuous.
     *
     * @param name string name of species
     * @param tmid  Mid temperature, between the two temperature regions
     * @param clow  coefficients for lower temperature region
     * @param chigh coefficients for higher temperature region
     */
    double checkContinuity(const std::string& name, double tmid,
                           doublereal* clow, doublereal* chigh);

    //! Adjust polynomials to be continuous at the midpoint temperature.
    /*!
     * We seek a set of coefficients for the low- and high-temperature
     * polynomials which are continuous in Cp, H, and S at the midpoint while
     * minimizing the difference between the values in Cp, H, and S over the
     * entire valid temperature range. To do this, we formulate a linear
     * least-squares problem to be solved for 11 of the 14 coefficients, with
     * the remaining 3 coefficients eliminated in the process of satisfying
     * the continuity constraints.
     *
     * @param Tlow  Minimum temperature at which the low-T polynomial is valid
     * @param Tmid  Mid temperature, between the two temperature regions
     * @param Thigh Maximum temperature at which the high-T polynomial is valid
     * @param clow  coefficients for lower temperature region
     * @param chigh coefficients for higher temperature region
     */
    void fixDiscontinuities(doublereal Tlow, doublereal Tmid, doublereal Thigh,
                            doublereal* clow, doublereal* chigh);
};

}

#endif
