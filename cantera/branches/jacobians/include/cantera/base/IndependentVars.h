/**
 * @file IndependentVars.h
 * Header for specification of independent Variables within Cantera
 */

#ifndef CT_INDEPENDENTVARS_H
#define CT_INDEPENDENTVARS_H

#include "ct_defs.h"

namespace Cantera {

// Forward declaration
template<typename ValAndDerivType> class ThermoPhase;

//!  Specification of the types of independent variable combinations allowed within Cantera
/*!
 *   Independent variables are defined by the calling application as the variables that are
 *   varied by the application in order to find the solution to the application's problem. 
 *   Typically the application needs the derivatives of various quantities supplied by Cantera
 *   with respect to these independent variables. 
 *
 *   They differ from Cantera's State variables, which are used to define the intrinsic state
 *   of an object. Once the intrinsic state of an object has been defined, the value of the
 *   object becomes independent of how it was originally defined.
 *
 *   For gradient quantities calculated by Cantera, Cantera requires the State variable value of a phase at
 *   two locations given by certain distance apart, or the value of one state variable specified state, 
 *   and a spacial gradient of the phase given by variations of some quantities that specifiy a valid
 *   internal state of the phase. After this specification is done, the value of the gradients become
 *   indendent of how it was originally defined.
 *
 *   The Jacobian of the cantera's state variables and gradients (defined above) are dependent on the
 *   specification of the indepdent variables. They can't be separated from this specification. Therefore,
 *   they are independently defined from the State variables. 
 *
 *   Cantera will accept what the application needs for independent variables
 *
 *   Cantera will communicate the Jacobian values with the application using a vector to  pass results back and forth.
 *   Some testing has shown that it is important to avoid unnecessary scatters and gathers during this operation.
 *          
 *   The form of the vector will depend on the independent variables. However, it is generally as shown  
 *   
 *             0          temperature
 *             1          pressure
 *             X_0_p0     Mole fraction of species 0 in the first phase 0
 *             X_i_p0     Mole fraction of species i in the first phase 0
 *             X_N_p0     Mole fraction of species N in the first phase 0
 *             X_0_p1     Mole fraction of species N in the next phase 1
 *             X_i_p1     Mole fraction of species i in the next phase 1
 *             X_N_p1     Mole fraction of species N in the next phase 1
 *             V          Voltage (if needed)
 *             Pi_ps1     Surface tension (if needed) of first surface phase 0
 *             Pi_psn     Surface tension (if needed) of next surface phase 1
 *             Z_extern   Other unspecified external independent variables or parameters
 *
 *   The enum below sets the type of independent variable to be assocated as the default within the 
 *   application and within each phase.
 *
 */
enum IndependentVar_Formulation
{

    //!  Temperature, pressure, and complete mole fraction vector
    /*!
     *     Composition specified in terms of the temperature, pressure, and complete mole fraction
     *     vector, 1 to N. Note mole fractions are defined as independent of each other.
     *     Pressure units are in pascals and temperature is in kelvin.
     */
    INDVAR_TP_MOLEFRACTION_VECTOR = 0,

    //!  Temperature, pressure, and incomplete mole fraction vector
    /*!
     *     Composition specified in terms of the temperature, pressure, and incomplete mole fraction
     *     vector, 2 to N. Note mole fractions are defined as independent of each other except for
     *     a special species whose mole fraction moves in a compensating manner.
     *     Pressure units are in pascals and temperature is in kelvin.
     *
     *     The special species is defines initial as the first species, initially. It may be overridden. For example
     *     it is common to specify the species with largest mole fraction as the special species.
     */
    INDVAR_TP_MOLEFRACTIONNM1_VECTOR = 1,

    //!  Temperature and complete concentration vector
    /*!
     *     Composition specified in terms of the temperature and complete concentration
     *     vector, 1 to N. Note concentrations are defined as independent of each other.
     *     Concentrations are defined in units of kmol m-3.
     *     For surfaces and edge phases they are defined in terms of kmol m-2 and kmol-1.
     */
    INDVAR_T_CONCENTRATION_VECTOR = 2,

    //!  Temperature, pressure, and complete mole number vector
    /*!
     *     Composition specified in terms of the temperature, pressure, and complete mole number
     *     vector, 1 to N. Note mole numbers are defined as independent of each other.
     *     Pressure units are in pascals and temperature is in kelvin.
     */
    INDVAR_TP_MOLENUMBER_VECTOR = 3,

    //!  Temperature, pressure, and complete mass fraction vector
    /*!
     *     Composition specified in terms of the temperature, pressure, and complete mass fraction
     *     vector, 1 to N. Note mass fractions are defined as independent of each other.
     *     Pressure units are in pascals and temperature is in kelvin.
     */
    INDVAR_TP_MASSFRACTION_VECTOR = 4,

    //!  Temperature, pressure, and incomplete mass fraction vector
    /*!
     *     Composition specified in terms of the temperature, pressure, and incomplete mass fraction
     *     vector, 2 to N. Note mass fractions are defined as independent of each other except for
     *     a special species whose mass fraction moves in a compensating manner.
     *     Pressure units are in pascals and temperature is in kelvin.
     *
     *     The special species is defines initial as the first species, initially. It may be overridden. For example
     *     it is common to specify the species with largest mass fraction as the special species.
     */
    INDVAR_TP_MASSFRACTIONNM1_VECTOR = 5,

    //!  Temperature and complete density vector
    /*!
     *     Composition specified in terms of the temperature and complete density
     *     vector, 1 to N. Note concentrations are defined as independent of each other.
     *     Densities are defined in units of kg m-3.
     *     For surfaces and edge phases they are defined in terms of kg m-2 and kg m-1.
     */
    INDVAR_T_DENSITY_VECTOR = 6,

    //!  Undefined independent variable type. The problem starts out with this value
    INDVAR_UNDEFINED = -1

};

//! Short name for the enum type that specifies how unknowns are treated
typedef enum IndependentVar_Formulation INDVAR_FORM;

//! Variable Names
enum Var_Name_Enum
{
    VARIABLE_TYPE_UNSPECIFIED = -2, //!< Variable_Type_Unspecified
    VARIABLE_TYPE_ANY, //!< Variable_Type_Any
    TEMPERATURE, //!< Temperature
    PRESSURE, //!< Isotropic 3D Pressure of the phase
    PRESSUREGRADIENT_RADIAL, //!< PressureGradient_radial
    MOLEFRACTION_SPECIES, // 3
    MASSFRACTION_SPECIES, // 4
    VOLUMEFRACTION_SPECIES, //!< VolumeFraction_species
    VOLUMEFRACTION_PHASE, //!< VolumeFraction_phase
    CONCENTRATION_SPECIES, // adding in extra variable types
    VOLTAGE, // 7
    SURFACETENSION_PHASE, //!< Surface tension of the phase

    //! must be last in the list
    MAX_VAR_NAME //!< Max_Var_Name
};

//! simplifying name for the enum
typedef enum Var_Name_Enum VAR_TYPE;

//!  Typedef for the subtype of the variable
typedef size_t VAR_TYPE_SUBNUM;

//!  Specification of independent variables
/*!
 *  There exists one of these for each phase in the problem that needs a specification of an independent variables vector.
 */
class IndVar_Phase_Specification
{
public:
    //! Constructor
    IndVar_Phase_Specification();

    //! Destructor
    ~IndVar_Phase_Specification();

    //! Copy Constructor
    /*!
     * @param b object to be copied
     */
    IndVar_Phase_Specification(const IndVar_Phase_Specification &b);

    //! Assignment operator
    /*!
     * @param b object to be copied
     *
     * @return  returns a changeable reference to the current object
     */
    IndVar_Phase_Specification & operator=(const IndVar_Phase_Specification &b);

    //! Specification of the main way to identify independent variables.
    INDVAR_FORM indVar_Method_;

    //! Has a voltage
    bool hasVoltage_;

    //! Has a surface tension variable
    bool hasSurfaceTension_;

    //! Const pointer to the ThermoPhase object whose independent variables we are specifying.
    const ThermoPhase<doubleFAD> *tp_ptr_;

};

//!  Specification of independent variables
/*!
 *  There exists one of these per problem
 */
class IndVar_ProblemSpecification
{
public:
    //! Constructor
    IndVar_ProblemSpecification();

    //! Destructor
    ~IndVar_ProblemSpecification();

    //! Copy Constructor
    /*!
     * @param b object to be copied
     */
    IndVar_ProblemSpecification(const IndVar_ProblemSpecification &b);

    //! Assignment operator
    /*!
     * @param b object to be copied
     *
     * @return  returns a changeable reference to the current object
     */
    IndVar_ProblemSpecification & operator=(const IndVar_ProblemSpecification &b);

    //! Specification of the main way to identify independent variables.
    INDVAR_FORM indVar_Method_;

    //! Has a voltage
    bool hasVoltage_;

    //! Has a surface tension variable
    bool hasSurfaceTension_;

};

} // End of namespace Cantera

#endif

