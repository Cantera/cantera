/**
 * @file VarTypes.cpp
 *
 */

/*
 *  $Id: m1d_VarTypes.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "cantera/base/VarTypes.h"
#include "cantera/base/ctexceptions.h"

#include <cstring>
#include <cstdio>
#include <string>

using namespace std;

namespace Cantera {
//=====================================================================================================================
// Returns the variable main name given the variable type
/*
 *  Static function
 *
 *  @param variableType  The variable type
 */
std::string VarType::VarMainName(const VAR_TYPE variableType)
{
    switch (variableType) {

    case VARIABLE_TYPE_TEMPERATURE:
        return string("Temperature");
        break;
    case VARIABLE_TYPE_PRESSURE:
        return string("Pres_axial");
        break;
    case VARIABLE_TYPE_PRESSUREGRADIENT_RADIAL:
        return string("PresGrad_radial");
        break;
    case VARIABLE_TYPE_MOLEFRACTION_SPECIES:
        return string("MF_sp");
        break;
    case VARIABLE_TYPE_MASSFRACTION_SPECIES:
        return string("WF_sp");
        break;
    case VARIABLE_TYPE_VOLUMEFRACTION_SPECIES:
        return string("VF_sp");
        break;
    case VARIABLE_TYPE_VOLUMEFRACTION_PHASE:
        return string("VF_ph");
        break;
    case VARIABLE_TYPE_CONCENTRATION_SPECIES:
        return string("Conc_sp");
        break;
    case VARIABLE_TYPE_VOLTAGE:
        return string("Volt");
        break;
    case VARIABLE_TYPE_SURFACETENSION_PHASE:
        return string("SurfTn_ph");
        break;
    default:
        throw CanteraError("VarMainName", "unknown");
        break;
    }
    return string("");
}
//=====================================================================================================================
// Constructor
VarType::VarType() :
        VariableType(VARIABLE_TYPE_UNSPECIFIED),
        VariableSubType(npos)
{
    VariableMainName[0] = '\0';
    VariableSubTypeName[0] = '\0';
}
//===========================================================================
VarType::VarType(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType, const char *subName) :
        VariableType(variableType),
        VariableSubType(variableSubType)
{
    VariableMainName[0] = '\0';
    std::string h = VarMainName(VariableType);
    strncpy(VariableMainName, h.c_str(), 23);
    VariableSubTypeName[0] = '\0';
    if (subName) {
        strncpy(VariableSubTypeName, subName, 23);
    }
    VariableSubTypeName[23] = '\0';
}
//=====================================================================================================================
//! Destructor
VarType::~VarType()
{
}
//=====================================================================================================================
// Copy Constructor
/*
 * @param r Object to be copied
 */
VarType::VarType(const VarType &r) :
        VariableType(VARIABLE_TYPE_UNSPECIFIED),
        VariableSubType(npos)
{
    VariableMainName[0] = '\0';
    VariableSubTypeName[0] = '\0';
    operator=(r);
}
//=====================================================================================================================
// Assignment Operator
/*
 * @param r Object to be copied.
 * @return Returns a variable reference to the current object
 */
VarType &
VarType::operator=(const VarType &r)
{
    if (this == &r) {
        return *this;
    }
    VariableType = r.VariableType;
    VariableSubType = r.VariableSubType;
    strncpy(VariableMainName, r.VariableMainName, 23);
    strncpy(VariableSubTypeName, r.VariableSubTypeName, 23);
    return *this;
}
//=====================================================================================================================
void VarType::setID(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType, const char *subName)
{
    VariableType = variableType;
    VariableSubType = variableSubType;
    if (!subName) {
        strncpy(VariableSubTypeName, subName, 23);
    } else {
        VariableSubTypeName[0] = '\0';
    }
    VariableMainName[0] = '\0';
    VariableSubTypeName[23] = '\0';
}
//=====================================================================================================================
std::string VarType::VariableName(const int len) const
{
    char buf[128];
    if (strlen(VariableMainName) > 0) {
        strncpy(buf, VariableMainName, 63);
    } else {
        sprintf(buf, "Var_%d ", (int) VariableSubType);
    }
    int ll = strlen(buf);
    if (len > (ll + 3)) {
        int left = len - (ll);
        if (strlen(VariableSubTypeName) > 0) {
            snprintf(buf + ll, left, "(%s", VariableSubTypeName);
        } else {
            snprintf(buf + ll, left, "(%d", (int) VariableSubType);
        }
        int ll = strlen(buf);
        sprintf(buf + ll, ")");
    }
    return std::string(buf);
}
//=====================================================================================================================
bool operator==(const VarType &a, const VarType &b)
{
    if (a.VariableType == b.VariableType) {
        if (a.VariableSubType == b.VariableSubType) {
            return true;
        }
    }
    return false;
}
//=====================================================================================================================
bool operator!=(const VarType &a, const VarType &b)
{
    if (a.VariableType == b.VariableType) {
        if (a.VariableSubType == b.VariableSubType) {
            return false;
        }
    }
    return true;
}
//=====================================================================================================================
}
/* end of namespace */
//=====================================================================================================================
