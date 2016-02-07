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

#include "SpeciesThermo.h"
#include "SpeciesThermoInterpType.h"

namespace Cantera
{

//! A species thermodynamic property manager for a phase.
/*!
 * This is a general manager that can handle a wide variety of species
 * thermodynamic polynomials for individual species. It is slow, however,
 * because it recomputes the functions of temperature needed for each species.
 * What it does is to create a vector of SpeciesThermoInterpType objects.
 *
 * @ingroup mgrsrefcalc
 */
class GeneralSpeciesThermo : public SpeciesThermo
{
public:
    //! Constructor
    GeneralSpeciesThermo();

    GeneralSpeciesThermo(const GeneralSpeciesThermo& b);
    GeneralSpeciesThermo& operator=(const GeneralSpeciesThermo& b);
    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const;

    virtual void install_STIT(size_t index,
                              shared_ptr<SpeciesThermoInterpType> stit_ptr);

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
     * @return pointer to the SpeciesThermoInterpType object.
     */
    SpeciesThermoInterpType* provideSTIT(size_t k);
    const SpeciesThermoInterpType* provideSTIT(size_t k) const;

protected:
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

    std::map<size_t, std::pair<int, size_t> > m_speciesLoc;

    //! Maximum value of the lowest temperature
    doublereal m_tlow_max;

    //! Minimum value of the highest temperature
    doublereal m_thigh_min;

    //! reference pressure (Pa)
    doublereal m_p0;

    //! Make the class VPSSMgr a friend because we need to access the function
    //! provideSTIT()
    friend class VPSSMgr;
};

}

#endif
