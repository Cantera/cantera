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

#include "SpeciesThermoMgr.h"
#include "SpeciesThermoInterpType.h"

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
     * @param refPressure standard-state pressure for this parameterization.
     * @see speciesThermoTypes.h
     *
     * @deprecated Use newSpeciesThermoInterpType and
     *     GeneralSpeciesThermo::install_STIT. To be removed after Cantera 2.2.
     */
    virtual void install(const std::string& name, size_t index, int type,
                         const doublereal* c,
                         doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure);

    virtual void install_STIT(size_t index, SpeciesThermoInterpType* stit_ptr);

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
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void update_one(size_t k, doublereal T, doublereal* cp_R,
                            doublereal* h_RT,
                            doublereal* s_R) const;

    virtual void update(doublereal T, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const;

    virtual doublereal minTemp(size_t k=npos) const;
    virtual doublereal maxTemp(size_t k=npos) const;
    virtual doublereal refPressure(size_t k=npos) const;
    virtual int reportType(size_t index) const;

    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const;

    virtual doublereal reportOneHf298(const size_t k) const;

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New);

private:
    //! Provide the SpeciesthermoInterpType object
    /*!
     * @param k  species index
     *
     * @return pointer to the SpeciesThermoInterpType object.
     */
    SpeciesThermoInterpType* provideSTIT(size_t k);
    const SpeciesThermoInterpType* provideSTIT(size_t k) const;

    void clear(); //<! Delete owned SpeciesThermoInterpType objects.

protected:
    typedef std::map<int, std::vector<SpeciesThermoInterpType*> > STIT_map;
    typedef std::map<int, std::vector<double> > tpoly_map;
    /**
     * This is the main unknown in the object. It contains pointers to
     * SpeciesThermoInterpType objects, sorted by the parameterization type.
     * This object owns the SpeciesThermoInterpType objects, so they are deleted
     * in the destructor of this object.
     */
    STIT_map m_sp;

    //! Temperature polynomials for each thermo parameterization
    mutable tpoly_map m_tpoly;

    std::map<size_t, std::pair<int, size_t> > m_speciesLoc;

    //! Maximum value of the lowest temperature
    doublereal m_tlow_max;

    //! Minimum value of the highest temperature
    doublereal m_thigh_min;

    //! reference pressure (Pa)
    doublereal m_p0;

    //! Make the class VPSSMgr a friend because we need to access
    //! the function provideSTIT()
    friend class VPSSMgr;
};

}

#endif
