/**
 * @file m1d_EqnVarTypes.h
 *
 */

/*
 *  $Id: EqnVarTypes.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef CT_VARTYPES_H
#define CT_VARTYPES_H

#include "IndependentVars.h"

#include <string>

namespace Cantera {

//! This class contains a name of a variable
/*!
 *   The class is used to identify variables.
 */
class VarType
{
public:

    //! Constructor
    VarType();

    //! Constructor with specification
    /*!
     * @param variableType
     * @param variableSubType
     * @param subName
     */
    VarType(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = 0);

    //! Destructor
    ~VarType();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    VarType(const VarType &r);

    //! Assignment Operator
    /*!
     * @param r Object to be copied.
     * @return Returns a variable reference to the current object
     */
    VarType & operator=(const VarType &r);

    //! Returns the variable main name given the variable type
    /*!
     *  Static function
     *
     *  @param variableType  The variable type
     */
    static std::string VarMainName(const VAR_TYPE variableType);

    //! Set the variable type of the object
    /*!
     *
     * @param variableType
     * @param variableSubType
     * @param subName
     */
    void setID(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = 0);

    std::string VariableName(const int len) const;

    //! Variable type
    VAR_TYPE VariableType;

    //! Variable subtype
    VAR_TYPE_SUBNUM VariableSubType;

    //!  Main Part of the variable Name
    char VariableMainName[24];

    //! Sub type name of the variable
    char VariableSubTypeName[24];
};

//! Equality boolean operator for VarType Objects
/*!
 *
 * @param a Object 1
 * @param b Object 2
 * @return Returns whether the two are the same
 */
bool
operator==(const VarType &a, const VarType &b);

//! Inequality boolean operator for VarType Objects
/*!
 *
 * @param a Object 1
 * @param b Object 2
 * @return Returns whether the two are the same
 */
bool
operator!=(const VarType &a, const VarType &b);

}
#endif

