/**
 *  @file ThermoPhase.cpp
 * Definition file for class ThermoPhase, the base class for phases with
 * thermodynamic properties
 * (see class \link Cantera::ThermoPhase ThermoPhase\endlink).
 */

//  Copyright 2002 California Institute of Technology
#include "cantera/thermo/ThermoDerivInfo.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/VarTypes.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera {
//====================================================================================================================

template<typename ValAndDerivType>
ThermoDerivInfo<ValAndDerivType>::ThermoDerivInfo()
{
}


ThermoDerivInfo<doubleFAD>::ThermoDerivInfo() :
        indVar_Method_(INDVAR_TP_MOLEFRACTION_VECTOR),
        hasVoltage_(false),
        hasSurfaceTension_(false),
        tp_double_ptr_(0),
        numIndVars_(0),
        IndVarList_()
{
    IndVar_ProblemSpecification& id = DefaultIndVar();
    indVar_Method_ = id.indVar_Method_;
    hasVoltage_ = id.hasVoltage_;
    hasSurfaceTension_ = id.hasSurfaceTension_;
}

//====================================================================================================================

template<typename ValAndDerivType>
ThermoDerivInfo<ValAndDerivType>::~ThermoDerivInfo()
{
}


ThermoDerivInfo<doubleFAD>::~ThermoDerivInfo()
{
    if (tp_double_ptr_) {
        delete tp_double_ptr_;
        tp_double_ptr_ = 0;
    }
}

//====================================================================================================================
/*
 * Copy Constructor for the ThermoPhase object.
 *
 * Currently, this is implemented, but not tested. If called it will
 * throw an exception until fully tested.
 */
template<typename ValAndDerivType>
ThermoDerivInfo<ValAndDerivType>::ThermoDerivInfo(const ThermoDerivInfo<ValAndDerivType>& right)
{
    /*
     * Call the assignment operator
     */
    operator=(right);
}
//====================================================================================================================

template<typename ValAndDerivType>
ThermoDerivInfo<ValAndDerivType>& ThermoDerivInfo<ValAndDerivType>::operator=(const ThermoDerivInfo<ValAndDerivType>& right)
{
    return *this;
}

ThermoDerivInfo<doubleFAD>& ThermoDerivInfo<doubleFAD>::operator=(const ThermoDerivInfo<doubleFAD>& right)
{
    if (this == &right) {
        return *this;
    }
    delete tp_double_ptr_;
    tp_double_ptr_ = 0;
    if (right.tp_double_ptr_) {
        tp_double_ptr_ = (right.tp_double_ptr_)->duplMyselfAsThermoPhase();
    }
    return *this;
}

//====================================================================================================================

template<typename ValAndDerivType>
template<typename ValAndDerivType2>
ThermoDerivInfo<ValAndDerivType>& ThermoDerivInfo<ValAndDerivType>::operator=(const ThermoDerivInfo<ValAndDerivType2>& right)
{
    return *this;
}

template<typename ValAndDerivType2>
ThermoDerivInfo<doubleFAD>& ThermoDerivInfo<doubleFAD>::operator=(const ThermoDerivInfo<ValAndDerivType2>& right)
{
    return *this;
}

//====================================================================================================================

template<typename ValAndDerivType>
void ThermoDerivInfo<ValAndDerivType>::generateThermoPhaseCopyDouble(const ThermoPhase<ValAndDerivType> *tpFAD)
{
}


void ThermoDerivInfo<doubleFAD>::generateThermoPhaseCopyDouble(const ThermoPhase<doubleFAD> *tpFAD)
{
    tp_double_ptr_ = tpFAD->duplMyselfAsThermoPhaseDouble();

    size_t nsp = tp_double_ptr_->nSpecies();

}
//====================================================================================================================

template<typename ValAndDerivType>
void ThermoDerivInfo<ValAndDerivType>::generateDerivVector()
{
}


void ThermoDerivInfo<doubleFAD>::generateDerivVector()
{

    size_t nsp = tp_double_ptr_->nSpecies();
    IndVarList_.push_back(VarType(VARIABLE_TYPE_TEMPERATURE));
    IndVarList_.push_back(VarType(VARIABLE_TYPE_PRESSURE));
    numIndVars_ = 2;
    for (size_t k = 0; k < nsp; k++) {
        string sn = tp_double_ptr_->speciesName(k);
        IndVarList_.push_back(VarType(VARIABLE_TYPE_MOLEFRACTION_SPECIES, k, sn.c_str()));
        numIndVars_++;
    }
}



//====================================================================================================================

//! Explicit Instantiations
template class ThermoDerivInfo<doublereal> ;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class ThermoDerivInfo<doubleFAD>;
template ThermoDerivInfo<doublereal>& ThermoDerivInfo<doublereal>::operator=(const ThermoDerivInfo<doubleFAD>& right);
#endif
#endif

}
