/**
 *  @file SpeciesThermoMgr.h
 *  This file contains descriptions of templated subclasses of
 *  the virtual base class, SpeciesThermo, which includes SpeciesThermoDuo
 *  (see \ref mgrsrefcalc and class
 *   \link Cantera::SpeciesThermoDuo SpeciesThermoDuo\endlink)
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_SPECIESTHERMO_MGR_H
#define CT_SPECIESTHERMO_MGR_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "SpeciesThermo.h"
#include "cantera/base/utilities.h"

namespace Cantera
{
//! Unknown species thermo manager string error
/*!
 * @ingroup mgrsrefcalc
 */
class UnknownSpeciesThermo : public CanteraError
{
public:
    //! constructor
    /*!
     * @param proc   name of the procecdure
     * @param type   unknown type
     */
    UnknownSpeciesThermo(const std::string& proc, int type) :
        CanteraError(proc, "Specified species parameterization type (" + int2str(type)
                     + ") does not match any known type.") {}

    //! Alternate constructor
    /*!
     * @param proc   name of the procecdure
     * @param stype   String name for the unknown type
     */
    UnknownSpeciesThermo(const std::string& proc, const std::string& stype) :
        CanteraError(proc, "Specified species parameterization type (" + stype
                     + ") does not match any known type.") {}
};

/**
 *  This species thermo manager requires that all species have one
 *  of two parameterizations.
 *
 * Note this seems to be a slow way to do things, and it may be on its way out.
 *
 * @deprecated To be removed after Cantera 2.2. Use GeneralSpeciesThermo instead.
 * @ingroup mgrsrefcalc
 */
template<class T1, class T2>
class SpeciesThermoDuo : public SpeciesThermo
{
public:
    //! Constructor
    SpeciesThermoDuo() {
        warn_deprecated("class SpeciesThermoDuo", "To be removed after "
            "Cantera 2.2. Use GeneralSpeciesThermo instead.");
    };

    //! copy constructor
    /*!
     * @param right Object to be copied
     */
    SpeciesThermoDuo(const SpeciesThermoDuo& right) {
        *this = right;
    }

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    SpeciesThermoDuo& operator=(const SpeciesThermoDuo& right);

    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const;

    virtual void install(const std::string& name, size_t sp, int type,
                         const doublereal* c,
                         doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure);

    virtual void install_STIT(size_t index,
                              shared_ptr<SpeciesThermoInterpType> stit_ptr) {
        throw CanteraError("install_STIT", "not implemented");
    }

    virtual void update(doublereal t, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const;

    virtual doublereal minTemp(size_t k = npos) const {
        return std::max(m_thermo1.minTemp(),m_thermo2.minTemp());
    }

    virtual doublereal maxTemp(size_t k = npos) const {
        return std::min(m_thermo1.maxTemp(), m_thermo2.maxTemp());
    }

    virtual doublereal refPressure(size_t k = npos) const {
        return m_p0;
    }

    virtual int reportType(size_t k) const;

    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const;

    virtual doublereal reportOneHf298(const size_t k) const {
        throw CanteraError("reportHF298", "unimplemented");
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New) {
        throw CanteraError("reportHF298", "unimplemented");
    }

private:

    //! Thermo Type 1
    T1 m_thermo1;
    //! Thermo Type 2
    T2 m_thermo2;
    //! Reference pressure
    doublereal m_p0;
    //! map from species to type
    std::map<size_t, int> speciesToType;
};


// ------------------------- cpp part of file -------------------------------------

// Definitions for the SpeciesThermoDuo<T1,T2> templated class

template<class T1, class T2>
SpeciesThermoDuo<T1, T2> &
SpeciesThermoDuo<T1, T2>::operator=(const SpeciesThermoDuo& right)
{
    if (&right == this) {
        return *this;
    }

    SpeciesThermo::operator=(right);
    m_thermo1 = right.m_thermo1;
    m_thermo2 = right.m_thermo2;
    m_p0      = right.m_p0;
    speciesToType = right.speciesToType;

    return *this;
}

template<class T1, class T2>
SpeciesThermo*
SpeciesThermoDuo<T1, T2>::duplMyselfAsSpeciesThermo() const
{
    return new SpeciesThermoDuo<T1,T2>(*this);
}

template<class T1, class T2>
void
SpeciesThermoDuo<T1, T2>::install(const std::string& name, size_t sp, int type,
                                  const doublereal* c,
                                  doublereal minTemp_,
                                  doublereal maxTemp_,
                                  doublereal refPressure_)
{
    m_p0 = refPressure_;
    if (type == m_thermo1.ID) {
        m_thermo1.install(name, sp, type, c, minTemp_, maxTemp_,
                          refPressure_);
        speciesToType[sp] = m_thermo1.ID;
    } else if (type == m_thermo2.ID) {
        m_thermo2.install(name, sp, type, c, minTemp_, maxTemp_,
                          refPressure_);
        speciesToType[sp] = m_thermo2.ID;
    } else {
        throw UnknownSpeciesThermo("SpeciesThermoDuo:install",type);
    }
    markInstalled(sp);
}

template<class T1, class T2>
void
SpeciesThermoDuo<T1, T2>::update(doublereal t, doublereal* cp_R,
                                 doublereal* h_RT, doublereal* s_R) const
{
    m_thermo1.update(t, cp_R, h_RT, s_R);
    m_thermo2.update(t, cp_R, h_RT, s_R);
}

template<class T1, class T2>
int
SpeciesThermoDuo<T1, T2>::reportType(size_t k) const
{
    return getValue(speciesToType, k, -1);
}

template<class T1, class T2>
void
SpeciesThermoDuo<T1, T2>::reportParams(size_t index, int& type,
                                       doublereal* const c,
                                       doublereal& minTemp_,
                                       doublereal& maxTemp_,
                                       doublereal& refPressure_) const
{
    int ctype = reportType(index);
    if (ctype == m_thermo1.ID) {
        m_thermo1.reportParams(index, type, c, minTemp_, maxTemp_,
                               refPressure_);
    } else if (ctype == m_thermo2.ID) {
        m_thermo2.reportParams(index, type, c, minTemp_, maxTemp_,
                               refPressure_);
    } else {
        throw CanteraError("SpeciesThermoDuo", "mismatched SpeciesThermoInterpType");
    }
}

}
#endif
