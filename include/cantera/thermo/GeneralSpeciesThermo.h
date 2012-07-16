/**
 * @file GeneralSpeciesThermo.h
 *  Headers for a completely general species thermodynamic property
 *  manager for a phase (see \ref mgrsrefcalc and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 *
 *  Because it is general, it is slow.
 */
#ifndef CT_GENERALSPECIESTHERMO_H
#define CT_GENERALSPECIESTHERMO_H
#include <string>
#include "cantera/base/ct_defs.h"
#include "SpeciesThermoMgr.h"
#include "NasaPoly1.h"
#include "Nasa9Poly1.h"
#include "speciesThermoTypes.h"


namespace Cantera
{

//! A species thermodynamic property manager for a phase.
/*!
 * This is a general manager that can handle a wide variety
 * of species thermodynamic polynomials for individual species.
 * It is slow, however, because it recomputes the functions of
 * temperature needed for each species. What it does is to create
 * a vector of SpeciesThermoInterpType objects.
 *
 * @ingroup mgrsrefcalc
 */
class GeneralSpeciesThermo : public SpeciesThermo
{

public:

    //! Constructor
    GeneralSpeciesThermo();

    //! Copy constructor
    /*!
     * @param b   Object to be copied
     */
    GeneralSpeciesThermo(const GeneralSpeciesThermo& b);

    //! Assignment operator
    /*!
     * @param b   Object to be copied
     */
    GeneralSpeciesThermo& operator=(const GeneralSpeciesThermo& b);

    //! Destructor
    virtual ~GeneralSpeciesThermo();

    //! Duplicator
    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const ;

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * Install a SpeciesThermoInterpType object for the species, index.
     * This routine contains an internal list of  SpeciesThermoInterpType
     * objects that it knows about. A factory-type lookup is done
     * to create the object.
     *
     * @param name      Name of the species
     * @param index     The 'update' method will update the property
     *                  values for this species
     *                  at position i index in the property arrays.
     * @param type      int flag specifying the type of parameterization to be
     *                 installed.
     * @param c        vector of coefficients for the parameterization.
     *                 This vector is simply passed through to the
     *                 parameterization constructor. Its length depends upon
     *                 the parameterization.
     * @param minTemp  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure standard-state pressure for this
     *                    parameterization.
     * @see speciesThermoTypes.h
     *
     * @todo Create a factory method for SpeciesThermoInterpType.
     *       That's basically what we are doing here.
     */
    virtual void install(std::string name, size_t index, int type,
                         const doublereal* c,
                         doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure);

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType* stit_ptr);

    //! Install a PDSS object to handle the reference state thermodynamics
    //! calculation
    /*!
     * @param k           species index
     * @param PDSS_ptr    Pressure dependent standard state (PDSS) object
     *                    that will handle the reference state calc
     * @param vpssmgr_ptr Pointer to the variable pressure standard state
     *                    manager that handles the PDSS object.
     */
    void installPDSShandler(size_t k, PDSS* PDSS_ptr, VPSSMgr* vpssmgr_ptr);

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update_one(size_t k, doublereal T, doublereal* cp_R,
                            doublereal* h_RT,
                            doublereal* s_R) const;

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of each of the standard states.
     *
     * @param T       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update(doublereal T, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const;

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the minimum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum
     * temperature for species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k=npos) const;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the maximum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum
     * temperature for parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(size_t k=npos) const;

    //! The reference-state pressure for species k.
    /*!
     *
     * returns the reference state pressure in Pascals for
     * species k. If k is left out of the argument list,
     * it returns the reference state pressure for the first
     * species.
     * Note that some SpeciesThermo implementations, such
     * as those for ideal gases, require that all species
     * in the same phase have the same reference state pressures.
     *
     * @param k Species Index
     */
    virtual doublereal refPressure(size_t k=npos) const;

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual int reportType(size_t index) const;

    //! This utility function reports back the type of
    //! parameterization and all of the parameters for the species, index.
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

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportOneHf298(int k) const;

    virtual void modifyOneHf298(const int k, const doublereal Hf298New);

#endif

private:
    //! Provide the SpeciesthermoInterpType object
    /*!
     *  provide access to the SpeciesThermoInterpType object.
     *  This
     *
     *  @param k  integer parameter
     *
     * @return pointer to the SpeciesThermoInterpType object.
     */
    SpeciesThermoInterpType* provideSTIT(size_t k);

protected:

    /**
     * This is the main unknown in the object. It is
     * a list of pointers to type SpeciesThermoInterpType.
     * Note, this object owns the objects, so they are deleted
     * in the destructor of this object.
     *   Note, that in some instances, m_sp[k] = 0, e.g., no
     * SpeciesThermoInterpType is installed for one or more
     * species. These cases must be handled by the calling
     * routine.
     */
    std::vector<SpeciesThermoInterpType*> m_sp;

    //! Maximum value of the lowest temperature
    doublereal                         m_tlow_max;

    //! Minimum value of the highest temperature
    doublereal                         m_thigh_min;

    //! reference pressure (Pa)
    doublereal                         m_p0;

    /**
     * Internal variable indicating the length of the
     * number of species in the phase.
     */
    size_t m_kk;


    //! Make the class VPSSMgr a friend because we need to access
    //! the function provideSTIT()
    friend class VPSSMgr;


};

}

#endif

