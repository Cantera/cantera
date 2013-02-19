/**
 * @file m1d_EqnVarTypes.h
 *
 */

/*
 *  $Id: EqnVarTypes.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef CT_DOUBLEVALANDDERIV_H
#define CT_DOUBLEVALANDDERIV_H

#include "IndependentVars.h"
#include "ct_defs.h"
#include "VarTypes.h"

#include <string>
#include <vector>

namespace Cantera {

//! This class contains a name of a variable
/*!
 *   The class is used to identify variables.
 */
class doubleValAndDeriv : public doubleFAD
{
public:

    //! Constructor
    doubleValAndDeriv();

    //! Destructor: virtual
    virtual ~doubleValAndDeriv();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    doubleValAndDeriv(const doubleValAndDeriv &r);

    //! Assignment Operator
    /*!
     * @param r Object to be copied.
     * @return Returns a variable reference to the current object
     */
    doubleValAndDeriv & operator=(const doubleValAndDeriv &r);


    //! Returns a vector of variable types
    /*!
     *
     */
    std::vector<VarType>& getIndVarTypesList();


};


}
#endif

