/**
 *  @file ThermoDerivInfo.h
 *   ThermoPhase object
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_THERMODERIVINFO_H
#define CT_THERMODERIVINFO_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/VarTypes.h"

namespace Cantera {

template<typename ValAndDerivType> class ThermoPhase;

/**  \addtogroup thermoprops */
/* @{
 */

template<typename ValAndDerivType>
class ThermoDerivInfo
{

public:

    /// Constructors
    ThermoDerivInfo();

    //! Copy Constructor
    ThermoDerivInfo(const ThermoDerivInfo<ValAndDerivType>&);

    //! Assignment operator
    ThermoDerivInfo<ValAndDerivType>& operator=(const ThermoDerivInfo<ValAndDerivType>&);

    //! Assignment operator
    template<typename ValAndDerivType2>
    ThermoDerivInfo<ValAndDerivType>& operator=(const ThermoDerivInfo<ValAndDerivType2>&);

    /// Destructor.
    virtual ~ThermoDerivInfo();

    //! Generate copy
    virtual void generateThermoPhaseCopyDouble(const ThermoPhase<ValAndDerivType> *tpFAD);

    //! Generate the vector of independent unknowns as a vector of VarTypes.
    virtual void generateDerivVector();

};

template<>
class ThermoDerivInfo<doubleFAD>
{

public:

    /// Constructors
    ThermoDerivInfo();

    //! Copy Constructor
    ThermoDerivInfo(const ThermoDerivInfo<doubleFAD>&);

    //! Assignment operator
    ThermoDerivInfo<doubleFAD>& operator=(const ThermoDerivInfo<doubleFAD>&);

    //! Assignment operator
    template<typename ValAndDerivType2>
    ThermoDerivInfo<doubleFAD>& operator=(const ThermoDerivInfo<ValAndDerivType2>&);

    //! Destructor.
    virtual ~ThermoDerivInfo();

    //! Generate copy
    virtual void generateThermoPhaseCopyDouble(const ThermoPhase<doubleFAD> *tpFAD);

    //! Generate the vector of independent unknowns as a vector of VarTypes.
    virtual void generateDerivVector();

protected:

    //! Specification of the main way to identify independent variables.
    INDVAR_FORM indVar_Method_;

    //! Has a voltage
    bool hasVoltage_;

    //! Has a surface tension variable
    bool hasSurfaceTension_;

    //! Pointer to a slave ThermoPhase object that gets used
    ThermoPhase<doublereal>* tp_double_ptr_;

    //! Number of independent variables
    size_t numIndVars_;

    //! Vector of VarTypes that represents the independent variables in the problem
    std::vector<VarType> IndVarList_;

};

/* @} */
}

#endif

