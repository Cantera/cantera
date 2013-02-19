/**
 *  @file GeneralSpeciesThermo.cpp
 *  Declarations for a completely general species thermodynamic property
 *  manager for a phase (see \ref spthermo and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 */
// Copyright 2001-2004  California Institute of Technology
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/NasaPoly1.h"
#include "NasaPoly2.h"
#include "ShomatePoly.h"
#include "ConstCpPoly.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/thermo/AdsorbateThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include <iostream>

using namespace std;

namespace Cantera {
//=====================================================================================================================
/*
 * Constructors
 */

template<typename ValAndDerivType>
GeneralSpeciesThermo<ValAndDerivType>::GeneralSpeciesThermo() :
        SpeciesThermo<ValAndDerivType>(),
        m_tlow_max(0.0),
        m_thigh_min(1.0E30),
        m_p0(OneAtm),
        m_kk(0)
{
    m_tlow_max = 0.0;
    m_thigh_min = 1.0E30;
}
//=====================================================================================================================

template<typename ValAndDerivType>
GeneralSpeciesThermo<ValAndDerivType>::GeneralSpeciesThermo(const GeneralSpeciesThermo<ValAndDerivType> & b) :
        m_tlow_max(b.m_tlow_max),
        m_thigh_min(b.m_thigh_min),
        m_p0(OneAtm),
        m_kk(b.m_kk)
{
    m_sp.resize(m_kk, 0);
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType<ValAndDerivType>* bk = b.m_sp[k];
        if (bk) {
            m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
        }
    }
}

template<typename ValAndDerivType>
template<typename ValAndDerivType2>
GeneralSpeciesThermo<ValAndDerivType>::GeneralSpeciesThermo(const GeneralSpeciesThermo<ValAndDerivType2> & b) :
        m_tlow_max(b.m_tlow_max),
        m_thigh_min(b.m_thigh_min),
        m_p0(OneAtm),
        m_kk(b.m_kk)
{
    m_sp.resize(m_kk, 0);
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType<ValAndDerivType2>* bk = b.m_sp[k];
        if (bk) {
            m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
        }
    }
}

template<>
template<>
GeneralSpeciesThermo<doublereal>::GeneralSpeciesThermo(const GeneralSpeciesThermo<doubleFAD> & b) :
        m_tlow_max(b.m_tlow_max),
        m_thigh_min(b.m_thigh_min),
        m_p0(OneAtm),
        m_kk(b.m_kk)
{
    m_sp.resize(m_kk, 0);
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType<doubleFAD>* bk = b.m_sp[k];
        if (bk) {
            m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpTypeDouble();
        }
    }
}

//=====================================================================================================================
template<typename ValAndDerivType>
template<typename ValAndDerivType2>
GeneralSpeciesThermo<ValAndDerivType>&
GeneralSpeciesThermo<ValAndDerivType>::operator=(const GeneralSpeciesThermo<ValAndDerivType2> & b)
{
    if ((GeneralSpeciesThermo<ValAndDerivType> *) &b != this) {
        m_tlow_max = b.m_tlow_max;
        m_thigh_min = b.m_thigh_min;

        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[k];
            if (sp) {
                delete sp;
                m_sp[k] = 0;
            }
        }
        m_kk = b.m_kk;
        m_sp.resize(m_kk, 0);
        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<ValAndDerivType2>* bk = b.m_sp[k];
            if (bk) {
                m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
            }
        }
    }
    return *this;
}

template<>
template<>
GeneralSpeciesThermo<doublereal>&
GeneralSpeciesThermo<doublereal>::operator=(const GeneralSpeciesThermo<doublereal> & b)
{
    if (&b != this) {
        m_tlow_max = b.m_tlow_max;
        m_thigh_min = b.m_thigh_min;

        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<doublereal>* sp = m_sp[k];
            if (sp) {
                delete sp;
                m_sp[k] = 0;
            }
        }
        m_kk = b.m_kk;
        m_sp.resize(m_kk, 0);
        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<doublereal>* bk = b.m_sp[k];
            if (bk) {
                m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
            }
        }
    }
    return *this;
}

template<>
template<>
GeneralSpeciesThermo<doublereal>&
GeneralSpeciesThermo<doublereal>::operator=(const GeneralSpeciesThermo<doubleFAD> & b)
{
    if ((GeneralSpeciesThermo<doublereal> *) &b != this) {
        m_tlow_max = b.m_tlow_max;
        m_thigh_min = b.m_thigh_min;

        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<doublereal>* sp = m_sp[k];
            if (sp) {
                delete sp;
                m_sp[k] = 0;
            }
        }
        m_kk = b.m_kk;
        m_sp.resize(m_kk, 0);
        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType<doubleFAD>* bk = b.m_sp[k];
            if (bk) {
                m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpTypeDouble();
            }
        }
    }
    return *this;
}

template<typename ValAndDerivType>
GeneralSpeciesThermo<ValAndDerivType>::~GeneralSpeciesThermo()
{
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[k];
        if (sp) {
            delete sp;
            m_sp[k] = 0;
        }
    }
}

template<typename ValAndDerivType>
SpeciesThermo<ValAndDerivType> *
GeneralSpeciesThermo<ValAndDerivType>::duplMyselfAsSpeciesThermo() const
{
    return new GeneralSpeciesThermo<ValAndDerivType>(*this);
}
//=====================================================================================================================
// Duplication routine for objects which inherit from %SpeciesThermo
/*
 *  This virtual routine can be used to duplicate %SpeciesThermo  objects
 *  inherited from %SpeciesThermo even if the application only has
 *  a pointer to %SpeciesThermo to work with.
 *
 *  This routine returns a doublereal templated version of SpeciesThermo no matter
 *  what templated version the underlying class is.
 *
 *  @return Duplicated <double> version of the SpeciesThermo
 */
template<typename ValAndDerivType>
SpeciesThermo<doublereal>* GeneralSpeciesThermo<ValAndDerivType>::duplMyselfAsSpeciesThermoDouble() const
{
    GeneralSpeciesThermo<doublereal>* nn = new GeneralSpeciesThermo<doublereal>(*this);
    return (SpeciesThermo<doublereal>*) nn;
}
//=====================================================================================================================
/*
 * Install parameterization for a species.
 * @param index    Species index
 * @param type     parameterization type
 * @param c        coefficients. The meaning of these depends on
 *                 the parameterization.
 */
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::install(const std::string& name, size_t index, int type, const doublereal* c,
                                                    doublereal minTemp, doublereal maxTemp, doublereal refPressure)
{
    /*
     * Resize the arrays if necessary, filling the empty
     * slots with the zero pointer.
     */

    if (minTemp <= 0.0) {
        throw CanteraError("Error in GeneralSpeciesThermo.cpp", " Cannot take 0 tmin as input. \n\n");
    }

    if (index >= m_kk) {
        m_sp.resize(index + 1, 0);
        m_kk = index + 1;
    }

    //int nfreq = 3;
    /*
     * Create the necessary object
     */

    switch (type) {
    case NASA1:
        m_sp[index] = new NasaPoly1<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    case SHOMATE1:
        m_sp[index] = new ShomatePoly<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    case CONSTANT_CP:
    case SIMPLE:
        m_sp[index] = new ConstCpPoly<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    case MU0_INTERP:
        m_sp[index] = new Mu0Poly<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    case SHOMATE2:
        m_sp[index] = new ShomatePoly2<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    case NASA2:
        m_sp[index] = new NasaPoly2<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;

    case STAT:
        m_sp[index] = new StatMech<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c, name);
        break;

    case ADSORBATE:
        m_sp[index] = new Adsorbate<ValAndDerivType>(index, minTemp, maxTemp, refPressure, c);
        break;
    default:
        throw UnknownSpeciesThermoModel("GeneralSpeciesThermo::install", "unknown species type", int2str(type));
        break;
    }
    if (!m_sp[index]) {
        cout << "Null m_sp... index = " << index << endl;
        cout << "type = " << type << endl;
    }
    m_tlow_max = max(minTemp, m_tlow_max);
    m_thigh_min = min(maxTemp, m_thigh_min);
}
//=====================================================================================================================
// Install a new species thermodynamic property
// parameterization for one species.
/*
 * @param stit_ptr Pointer to the SpeciesThermoInterpType object
 *          This will set up the thermo for one species
 */
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::install_STIT(SpeciesThermoInterpType<ValAndDerivType> * stit_ptr)
{
    /*
     * Resize the arrays if necessary, filling the empty
     * slots with the zero pointer.
     */
    if (!stit_ptr) {
        throw CanteraError("GeneralSpeciesThermo::install_STIT", "zero pointer");
    }
    size_t index = stit_ptr->speciesIndex();
    if (index >= m_kk) {
        m_sp.resize(index + 1, 0);
        m_kk = index + 1;
    }
    AssertThrow(m_sp[index] == 0, "Index position isn't null, duplication of assignment: " + int2str(index));
    /*
     *  Now, simply assign the position
     */
    m_sp[index] = stit_ptr;

    /*
     * Calculate max and min
     */
    double minTemp = stit_ptr->minTemp();
    double maxTemp = stit_ptr->maxTemp();

    m_tlow_max = max(minTemp, m_tlow_max);
    m_thigh_min = min(maxTemp, m_thigh_min);
}
//=====================================================================================================================
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::installPDSShandler(size_t k, PDSS* PDSS_ptr, VPSSMgr* vpssmgr_ptr)
{
    STITbyPDSS<ValAndDerivType>* stit_ptr = new STITbyPDSS<ValAndDerivType>(k, vpssmgr_ptr, PDSS_ptr);
    install_STIT(stit_ptr);
}
//=====================================================================================================================
/**
 *  Update the properties for one species.
 */
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::update_one(size_t k, doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                       ValAndDerivType* s_R) const
{
    SpeciesThermoInterpType<ValAndDerivType>* sp_ptr = m_sp[k];
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
}
//=====================================================================================================================
/*
 *  Update the properties for all species.
 */
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::update(doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                   ValAndDerivType* s_R) const
{
    typename std::vector<SpeciesThermoInterpType<ValAndDerivType>*>::const_iterator _begin, _end;
    _begin = m_sp.begin();
    _end = m_sp.end();
    SpeciesThermoInterpType<ValAndDerivType>* sp_ptr = 0;
    for (; _begin != _end; ++_begin) {
        sp_ptr = * (_begin);
        if (sp_ptr) {
            sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
        }
        // else {
        //   writelog("General::update: sp_ptr is NULL!\n");
        //}
    }
}
//=====================================================================================================================
/**
 * This utility function reports the type of parameterization
 * used for the species, index.
 */
template<typename ValAndDerivType>
int GeneralSpeciesThermo<ValAndDerivType>::reportType(size_t index) const
{
    SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[index];
    if (sp) {
        return sp->reportType();
    }
    return -1;
}
//=====================================================================================================================
/**
 * This utility function reports back the type of
 * parameterization and all of the parameters for the
 * species, index.
 *  For the NASA object, there are 15 coefficients.
 */
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::reportParams(size_t index, int& type, doublereal* const c, doublereal& minTemp,
                                                         doublereal& maxTemp, doublereal& refPressure) const
{
    SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[index];
    size_t n;
    if (sp) {
        sp->reportParameters(n, type, minTemp, maxTemp, refPressure, c);
        if (n != index) {
            throw CanteraError("GeneralSpeciesThermo::reportParams", "Internal error encountered");
        }
    } else {
        type = -1;
    }
}
//=====================================================================================================================
// //! Modify parameters for the standard state
// /*!
//  * @param index Species index
//  * @param c     Vector of coefficients used to set the
//  *              parameters for the standard state.
//  */
// void GeneralSpeciesThermo::
// modifyParams(size_t index, doublereal* c)
// {
//   SpeciesThermoInterpType* sp = m_sp[index];
//   if (sp) {
//     sp->modifyParameters(c);
//   }
// }
//=====================================================================================================================
/**
 * Return the lowest temperature at which the thermodynamic
 * parameterization is valid.  If no argument is supplied, the
 * value is the one for which all species parameterizations
 * are valid. Otherwise, if an integer argument is given, the
 * value applies only to the species with that index.
 */
template<typename ValAndDerivType>
doublereal GeneralSpeciesThermo<ValAndDerivType>::minTemp(size_t k) const
{
    if (k == npos) {
        return m_tlow_max;
    } else {
        SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[k];
        if (sp) {
            return sp->minTemp();
        }
    }
    return m_tlow_max;
}
//=====================================================================================================================
template<typename ValAndDerivType>
doublereal GeneralSpeciesThermo<ValAndDerivType>::maxTemp(size_t k) const
{
    if (k == npos) {
        return m_thigh_min;
    } else {
        SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[k];
        if (sp) {
            return sp->maxTemp();
        }
    }
    return m_thigh_min;
}
//=====================================================================================================================
template<typename ValAndDerivType>
doublereal GeneralSpeciesThermo<ValAndDerivType>::refPressure(size_t k) const
{
    if (k == npos) {
        return m_p0;
    } else {
        SpeciesThermoInterpType<ValAndDerivType>* sp = m_sp[k];
        if (sp) {
            return sp->refPressure();
        }
    }
    return m_p0;
}
//=====================================================================================================================
template<typename ValAndDerivType>
SpeciesThermoInterpType<ValAndDerivType> * GeneralSpeciesThermo<ValAndDerivType>::provideSTIT(size_t k)
{
    return (m_sp[k]);
}
//=====================================================================================================================
#ifdef H298MODIFY_CAPABILITY

template<typename ValAndDerivType>
doublereal GeneralSpeciesThermo<ValAndDerivType>::reportOneHf298(size_t k) const
{
    SpeciesThermoInterpType<ValAndDerivType>* sp_ptr = m_sp[k];
    doublereal h = -1.0;
    if (sp_ptr) {
        h = sp_ptr->reportHf298(0);
    }
    return h;
}
//=====================================================================================================================
template<typename ValAndDerivType>
void GeneralSpeciesThermo<ValAndDerivType>::modifyOneHf298(const size_t k, const doublereal Hf298New)
{
    SpeciesThermoInterpType<ValAndDerivType>* sp_ptr = m_sp[k];
    if (sp_ptr) {
        sp_ptr->modifyOneHf298(k, Hf298New);
    }
}

#endif
//=====================================================================================================================
// Explicit Instantiation Section

template class GeneralSpeciesThermo<doublereal> ;

//template  GeneralSpeciesThermo<doublereal>::GeneralSpeciesThermo(const GeneralSpeciesThermo<doublereal> & b);
template GeneralSpeciesThermo<doublereal>& GeneralSpeciesThermo<doublereal>::operator=(const GeneralSpeciesThermo<doublereal> & b);

#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class GeneralSpeciesThermo<doubleFAD> ;
//template  GeneralSpeciesThermo<doubleFAD>::GeneralSpeciesThermo(const GeneralSpeciesThermo<doubleFAD> & b);
template GeneralSpeciesThermo<doubleFAD>& GeneralSpeciesThermo<doubleFAD>::operator=(const GeneralSpeciesThermo<doubleFAD> & b);

template GeneralSpeciesThermo<doublereal>& GeneralSpeciesThermo<doublereal>::operator=(const GeneralSpeciesThermo<doubleFAD> & b);

#endif
#endif

}
//=====================================================================================================================
