/**
 *  @file ShomatePoly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the Shomate temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::ShomatePoly ShomatePoly\endlink and
 *   \link Cantera::ShomatePoly2 ShomatePoly2\endlink).
 *    Shomate polynomial expressions.
 */
#include "cantera/thermo/SpeciesThermoMgr.h"
#include "cantera/Cantera.h"
#include "ShomateThermo.h"
#include "NasaThermo.h"
#include "cantera/thermo/SimpleThermo.h"

namespace Cantera {

//==============================================================================================================================

template<typename ValAndDerivType, template<typename > class T1, template<typename > class T2>
template<typename ValAndDerivType2>
SpeciesThermoDuo< ValAndDerivType, T1, T2> &
SpeciesThermoDuo<ValAndDerivType, T1, T2>::operator=( const SpeciesThermoDuo< ValAndDerivType2, T1, T2>& right)
{
    if ((SpeciesThermoDuo< ValAndDerivType, T1, T2> *) &right == this) {
        return *this;
    }
    /*
     * We are counting on their being the appropriate copy constructors defined here.
     */
    m_thermo1 = right.m_thermo1;
    m_thermo2 = right.m_thermo2;
    m_p0 = right.m_p0;
    speciesToType = right.speciesToType;

    return *this;
}

//============================================================================================================================================

template<typename ValAndDerivType, template<typename > class T1, template<typename > class T2>
SpeciesThermo<doublereal>*
SpeciesThermoDuo<ValAndDerivType, T1, T2>::duplMyselfAsSpeciesThermoDouble() const
{
    /*
     *  This is actually sufficient I believe for all cases
     */
    SpeciesThermoDuo< doublereal, T1, T2 >* nn = new SpeciesThermoDuo<doublereal, T1, T2>(*this);
    return (SpeciesThermo<doublereal>*) nn;
}

//============================================================================================================================================

template class SpeciesThermoDuo<doublereal, NasaThermo, ShomateThermo> ;
template class SpeciesThermoDuo<doublereal, NasaThermo, SimpleThermo> ;
template class SpeciesThermoDuo<doublereal, ShomateThermo, SimpleThermo> ;

#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class SpeciesThermoDuo<doubleFAD, NasaThermo, ShomateThermo> ;
template class SpeciesThermoDuo<doubleFAD, NasaThermo, SimpleThermo> ;
template class SpeciesThermoDuo<doubleFAD, ShomateThermo, SimpleThermo>;
#endif
#endif

//============================================================================================================================================

}
