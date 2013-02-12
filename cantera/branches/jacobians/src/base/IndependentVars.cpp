/**
 * @file VarTypes.cpp
 *
 */

/*
 *  $Id: m1d_VarTypes.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "cantera/base/IndependentVars.h"

#include <cstring>
#include <cstdio>
#include <string>

using namespace std;

namespace Cantera {

//============================================================================================================
IndVar_Phase_Specification::IndVar_Phase_Specification() :
        indVar_Method_(INDVAR_TP_MOLEFRACTION_VECTOR),
        hasVoltage_(false),
        hasSurfaceTension_(false),
        tp_ptr_(0)
{

}
//============================================================================================================
IndVar_Phase_Specification::~IndVar_Phase_Specification()
{
}
//============================================================================================================
// Copy Constructor
/*
 * @param b object to be copied
 */
IndVar_Phase_Specification::IndVar_Phase_Specification(const IndVar_Phase_Specification &b) :
        indVar_Method_(INDVAR_TP_MOLEFRACTION_VECTOR),
        hasVoltage_(false),
        hasSurfaceTension_(false),
        tp_ptr_(0)
{
    operator=(b);
}
//============================================================================================================

IndVar_Phase_Specification & IndVar_Phase_Specification::operator=(const IndVar_Phase_Specification &b)
{
    if (&b == this) {
        return *this;
    }
    indVar_Method_ = b.indVar_Method_;
    hasVoltage_ = b.hasVoltage_;
    hasSurfaceTension_ = b.hasSurfaceTension_;
    tp_ptr_ = b.tp_ptr_;

    return *this;
}
//=====================================================================================================================
IndVar_ProblemSpecification::IndVar_ProblemSpecification() :
        indVar_Method_(INDVAR_TP_MOLEFRACTION_VECTOR),
        hasVoltage_(false),
        hasSurfaceTension_(false)
{
}
//=====================================================================================================================
IndVar_ProblemSpecification::IndVar_ProblemSpecification(const IndVar_ProblemSpecification &b) :
        indVar_Method_(INDVAR_TP_MOLEFRACTION_VECTOR),
        hasVoltage_(false),
        hasSurfaceTension_(false)
{
    operator=(b);
}
//=====================================================================================================================
// Assignment operator
/*
 * @param b object to be copied
 *
 * @return  returns a changeable reference to the current object
 */
IndVar_ProblemSpecification & IndVar_ProblemSpecification::operator=(const IndVar_ProblemSpecification &b)
{
    if (&b == this) {
        return *this;
    }
    indVar_Method_ = b.indVar_Method_;
    hasVoltage_ = b.hasVoltage_;
    hasSurfaceTension_ = b.hasSurfaceTension_;

    return *this;
}
//=====================================================================================================================
IndVar_ProblemSpecification::~IndVar_ProblemSpecification()
{
}
//=====================================================================================================================
/* end of namespace */
}
//=====================================================================================================================
