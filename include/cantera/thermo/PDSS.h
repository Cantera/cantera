/**
 *  @file PDSS.h
 *    Declarations for the virtual base class PDSS (pressure dependent standard state)
 *    which handles calculations for a single species in a phase
 *    (see \ref pdssthermo and class \link Cantera::PDSS PDSS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_H
#define CT_PDSS_H
#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
/**
 * @defgroup pdssthermo Species Standard-State Thermodynamic Properties
 *
 * In this module we describe %Cantera's treatment of pressure dependent
 * standard states (PDSS) objects. These are objects that calculate the standard
 * state of a single species that depends on both temperature and pressure.
 *
 * To compute the thermodynamic properties of multicomponent solutions, it is
 * necessary to know something about the thermodynamic properties of the
 * individual species present in the solution. Exactly what sort of species
 * properties are required depends on the thermodynamic model for the solution.
 * For a gaseous solution (i.e., a gas mixture), the species properties required
 * are usually ideal gas properties at the mixture temperature and at a
 * reference pressure (almost always at 1 bar). For other types of solutions,
 * however, it may not be possible to isolate the species in a "pure" state. For
 * example, the thermodynamic properties of, say, Na+ and Cl- in saltwater are
 * not easily determined from data on the properties of solid NaCl, or solid Na
 * metal, or chlorine gas. In this case, the solvation in water is fundamental
 * to the identity of the species, and some other reference state must be used.
 * One common convention for liquid solutions is to use thermodynamic data for
 * the solutes in the limit of infinite dilution within the pure solvent;
 * another convention is to reference all properties to unit molality.
 *
 * In defining these standard states for species in a phase, we make the
 * following definition. A reference state is a standard state of a species in a
 * phase limited to one particular pressure, the reference pressure. The
 * reference state specifies the dependence of all thermodynamic functions as a
 * function of the temperature, in between a minimum temperature and a maximum
 * temperature. The reference state also specifies the molar volume of the
 * species as a function of temperature. The molar volume is a thermodynamic
 * function. A full standard state does the same thing as a reference state, but
 * specifies the thermodynamics functions at all pressures.
 *
 * Class PDSS is the base class for a family of classes that compute properties
 * of a single species in a phase at its standard states, for a range of
 * temperatures and pressures. PDSS objects are used by derivatives of the
 * VPStandardState class. These classes assume that there exists a standard
 * state for each species in the phase, where the thermodynamic functions are
 * specified as a function of temperature and pressure.  Standard state objects
 * for each species in the phase are all derived from the PDSS virtual base
 * class.
 *
 * The following classes inherit from PDSS. Each of these classes handles just
 * one species.
 *
 * - PDSS_IdealGas
 *   - standardState model = "IdealGas"
 *   - This model assumes that the species in the phase obeys the ideal gas law
 *     for their pressure dependence. The manager uses a SimpleThermo object to
 *     handle the calculation of the reference state. This object adds the
 *     pressure dependencies to the thermo functions.
 *
 * - PDSS_ConstVol
 *    - standardState model = "ConstVol" or "constant_incompressible"
 *    - This model assumes that the species in the phase obeys the constant
 *      partial molar volume pressure dependence. The manager uses a
 *      SimpleThermo object to handle the calculation of the reference state.
 *      This object adds the pressure dependencies to these thermo functions.
 *
 * - PDSS_SSVol
 *   - standardState model = "temperature_polynomial"
 *   - standardState model = "density_temperature_polynomial"
 *   - This model assumes that the species in the phase obey a fairly general
 *     equation of state, but one that separates out the calculation of the
 *     standard state density and/or volume. Models include a cubic polynomial
 *     in temperature for either the standard state volume or the standard state
 *     density. The manager uses a SimpleThermo object to handle the calculation
 *     of the reference state. This object then adds the pressure dependencies
 *     and the volume terms to these thermo functions to complete the
 *     representation.
 *
 * - PDSS_Water
 *   - standardState model = "Water"
 *   - This model assumes that Species 0 is assumed to be water, and a real
 *     equation of state is used to model the T, P behavior. Note, the model
 *     assumes that the species is liquid water, and not steam.
 *
 * - PDSS_HKFT
 *   - standardState model = "HKFT"
 *   - This model assumes that the species follows the HKFT pressure dependent
 *     equation of state
 *
 * Normally the PDSS object is not called directly. Instead the
 * VPStandardStateTP object manages the calls to the PDSS object for the entire
 * set of species that comprise a phase.
 *
 * The PDSS objects may or may not utilize a SpeciesThermoInterpType reference
 * state manager class to calculate the reference state thermodynamics functions
 * in their own calculation. There are some classes, such as PDSS_IdealGas and
 * PDSS+_ConstVol, which utilize the SpeciesThermoInterpType object because the
 * calculation is very similar to the reference state calculation, while there
 * are other classes, PDSS_Water and PDSS_HKFT, which don't utilize the
 * reference state calculation at all, because it wouldn't make sense to. For
 * example, using the PDSS_Water module, there isn't anything special about the
 * reference pressure of 1 bar, so the reference state calculation would
 * represent a duplication of work. Additionally, when evaluating thermodynamic
 * properties at higher pressures and temperatures, near the critical point,
 * evaluation of the thermodynamics at a pressure of 1 bar may lead to
 * situations where the liquid is unstable, i.e., beyond the spinodal curve
 * leading to potentially wrong evaluation results.
 *
 * @ingroup thermoprops
 */

class XML_Node;
class SpeciesThermoInterpType;
class VPStandardStateTP;

//! Virtual base class for a species with a pressure dependent standard state
/*!
 * Virtual base class for calculation of the pressure dependent standard state
 * for a single species
 *
 * Class PDSS is the base class for a family of classes that compute properties
 * of a set of species in their standard states at a range of temperatures and
 * pressures. The independent variables for this object are temperature and
 * pressure. The class may have a reference to a SpeciesThermoInterpType object
 * which handles the calculation of the reference state temperature behavior of
 * the species.
 *
 * This class is analogous to the SpeciesThermoInterpType class, except that
 * the standard state inherently incorporates the pressure dependence.
 *
 * The class operates on a setState temperature and pressure basis. It only
 * recalculates the standard state when the setState functions for temperature
 * and pressure are called.
 *
 * @ingroup pdssthermo
 */
class PDSS
{
public:
    //! @name Constructors
    //! @{

    //! Default Constructor
    PDSS();

    // PDSS objects are not copyable or assignable
    PDSS(const PDSS& b) = delete;
    PDSS& operator=(const PDSS& b) = delete;
    virtual ~PDSS() {}

    //! @}
    //! @name  Utilities
    //! @{

     //! @}
     //! @name Molar Thermodynamic Properties of the Species Standard State in
     //!     the Solution
     //! @{

    //! Return the molar enthalpy in units of J kmol-1
    /*!
     * @return the species standard state enthalpy in J kmol-1 at the current
     *     temperature and pressure.
     */
    virtual doublereal enthalpy_mole() const;

    //! Return the standard state molar enthalpy divided by RT
    /*!
     * @return The dimensionless species standard state enthalpy divided at
     *     the current temperature and pressure.
     */
    virtual doublereal enthalpy_RT() const;

    //! Return the molar internal Energy in units of J kmol-1
    /*!
     * @return The species standard state internal Energy in J kmol-1 at the
     *     current temperature and pressure.
     */
    virtual doublereal intEnergy_mole() const;

    //! Return the molar entropy in units of J kmol-1 K-1
    /*!
     * @return The species standard state entropy in J kmol-1 K-1 at the
     *     current temperature and pressure.
     */
    virtual doublereal entropy_mole() const;

    //! Return the standard state entropy divided by RT
    /*!
     * @return The species standard state entropy divided by RT at the current
     *     temperature and pressure.
     */
    virtual doublereal entropy_R() const;

    //! Return the molar Gibbs free energy in units of J kmol-1
    /*!
     * @return The species standard state Gibbs free energy in J kmol-1 at the
     * current temperature and pressure.
     */
    virtual doublereal gibbs_mole() const;

    //! Return the molar Gibbs free energy divided by RT
    /*!
     * @return The species standard state Gibbs free energy divided by RT at
     *     the current temperature and pressure.
     */
    virtual doublereal gibbs_RT() const;

    //! Return the molar const pressure heat capacity in units of J kmol-1 K-1
    /*!
     * @return The species standard state Cp in J kmol-1 K-1 at the current
     *     temperature and pressure.
     */
    virtual doublereal cp_mole() const;

    //! Return the molar const pressure heat capacity divided by RT
    /*!
     * @return The species standard state Cp divided by RT at the current
     *     temperature and pressure.
     */
    virtual doublereal cp_R() const;

    //! Return the molar const volume heat capacity in units of J kmol-1 K-1
    /*!
     * @return The species standard state Cv in J kmol-1 K-1 at the
     * current temperature and pressure.
     */
    virtual doublereal cv_mole() const;

    //! Return the molar volume at standard state
    /*!
     * @return The standard state molar volume at the current temperature and
     *     pressure. Units are m**3 kmol-1.
     */
    virtual doublereal molarVolume() const;

    //! Return the standard state density at standard state
    /*!
     * @return The standard state density at the current temperature and
     *     pressure. units are kg m-3
     */
    virtual doublereal density() const;

    //! Get the difference in the standard state enthalpy
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal enthalpyDelp_mole() const;

    //! Get the difference in the standard state entropy between
    //! the current pressure and the reference pressure, p0
    virtual doublereal entropyDelp_mole() const;

    //! Get the difference in the standard state Gibbs free energy
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal gibbsDelp_mole() const;

    //! Get the difference in standard state heat capacity
    //! between the current pressure and the reference pressure, p0.
    virtual doublereal cpDelp_mole() const;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    //! Return the reference pressure for this phase.
    doublereal refPressure() const {
        return m_p0;
    }

    //! return the minimum temperature
    doublereal minTemp() const {
        return m_minTemp;
    }

    //! return the minimum temperature
    doublereal maxTemp() const {
        return m_maxTemp;
    }

    //! Return the molar Gibbs free energy divided by RT at reference pressure
    /*!
     * @return The reference state Gibbs free energy at the current
     *     temperature, divided by RT.
     */
    virtual doublereal gibbs_RT_ref() const;

    //! Return the molar enthalpy divided by RT at reference pressure
    /*!
     * @return The species reference state enthalpy at the current
     *     temperature, divided by RT.
     */
    virtual doublereal enthalpy_RT_ref() const;

    //! Return the molar entropy divided by R at reference pressure
    /*!
     * @return The species reference state entropy at the current
     *     temperature, divided by R.
     */
    virtual doublereal entropy_R_ref() const;

    //! Return the molar heat capacity divided by R at reference pressure
    /*!
     * @return The species reference state heat capacity divided by R at the
     * current temperature.
     */
    virtual doublereal cp_R_ref() const;

    //! Return the molar volume at reference pressure
    /*!
     * @return The reference state molar volume. units are m**3 kmol-1.
     */
    virtual doublereal molarVolume_ref() const;

    //! @}
    //!  @name Mechanical Equation of State Properties
    //! @{

    //! Returns the pressure (Pa)
    virtual doublereal pressure() const;

    //! Sets the pressure in the object
    /*!
     * Currently, this sets the pressure in the PDSS object. It is indeterminant
     * what happens to the owning VPStandardStateTP object.
     *
     * @param   pres   Pressure to be set (Pascal)
     */
    virtual void setPressure(doublereal pres);

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     *     \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Set the internal temperature
    /*!
     * @param temp Temperature (Kelvin)
     */
    virtual void setTemperature(doublereal temp);

    //! Return the current stored temperature
    virtual doublereal temperature() const;

    //! Set the internal temperature and pressure
    /*!
     * @param  temp     Temperature (Kelvin)
     * @param  pres     pressure (Pascals)
     */
    virtual void setState_TP(doublereal temp, doublereal pres);

    //! Set the internal temperature and density
    /*!
     * @param  temp     Temperature (Kelvin)
     * @param  rho      Density (kg m-3)
     */
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //! @name  Miscellaneous properties of the standard state
    //! @{

    //! critical temperature
    virtual doublereal critTemperature() const;

    //! critical pressure
    virtual doublereal critPressure() const;

    //! critical density
    virtual doublereal critDensity() const;

    //! saturation pressure
    /*!
     *  @param T Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal T);

    //! Return the molecular weight of the species
    //! in units of kg kmol-1
    doublereal molecularWeight() const;

    //! Set the molecular weight of the species
    /*!
     * @param mw Molecular Weight in kg kmol-1
     */
    void setMolecularWeight(doublereal mw);

    //! @}
    //! @name Initialization of the Object
    //! @{

    //! Set the SpeciesThermoInterpType object used to calculate reference
    //! state properties
    void setReferenceThermo(shared_ptr<SpeciesThermoInterpType> stit) {
        m_spthermo = stit;
    }

    //! Set the parent VPStandardStateTP object of this PDSS object
    /*!
     * This information is only used by certain PDSS subclasses
     * @param phase   Pointer to the parent phase
     * @param k       Index of this species in the phase
     */
    virtual void setParent(VPStandardStateTP* phase, size_t k) {}

    //! Initialization routine
    /*!
     * This is a cascading call, where each level should call the the parent
     * level.
     */
    virtual void initThermo() {}

    //! Set model parameters from an AnyMap phase description, e.g. from the
    //! `equation-of-state` field of a species definition.
    void setParameters(const AnyMap& node) {
        m_input = node;
    }

    //! Initialization routine for the PDSS object based on the speciesNode
    /*!
     * This is a cascading call, where each level should call the the parent
     * level. This function is called before initThermo()
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setParametersFromXML(const XML_Node& speciesNode) {}

    //! This utility function reports back the type of parameterization and
    //! all of the parameters for the species, index.
    /*!
     * @param kindex    Species index (unused)
     * @param type      Integer type of the standard type (unused)
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(size_t& kindex, int& type, doublereal* const c,
                              doublereal& minTemp, doublereal& maxTemp,
                              doublereal& refPressure) const;

    //@}

protected:
    //! Current temperature used by the PDSS object
    mutable doublereal m_temp;

    //! State of the system - pressure
    mutable doublereal m_pres;

    //! Reference state pressure of the species.
    doublereal m_p0;

    //! Minimum temperature
    doublereal m_minTemp;

    //! Maximum temperature
    doublereal m_maxTemp;

    //! Molecular Weight of the species
    doublereal m_mw;

    //! Input data supplied via setParameters. This may include parameters for
    //! different phase models, which will be used when initThermo() is called.
    AnyMap m_input;

    //! Pointer to the species thermodynamic property manager. Not used in all
    //! PDSS models.
    shared_ptr<SpeciesThermoInterpType> m_spthermo;
};

//! Base class for PDSS classes which compute molar properties directly
class PDSS_Molar : public virtual PDSS
{
public:
    virtual doublereal enthalpy_RT() const;
    virtual doublereal entropy_R() const;
    virtual doublereal gibbs_RT() const;
    virtual doublereal cp_R() const;
};

//! Base class for PDSS classes which compute nondimensional properties directly
class PDSS_Nondimensional : public virtual PDSS
{
public:
    PDSS_Nondimensional();

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;

    virtual double enthalpy_RT_ref() const;
    virtual double entropy_R_ref() const;
    virtual double gibbs_RT_ref() const;
    virtual double cp_R_ref() const;
    virtual double molarVolume_ref() const;
    virtual double enthalpy_RT() const;
    virtual double entropy_R() const;
    virtual double gibbs_RT() const;
    virtual double cp_R() const;
    virtual double molarVolume() const;
    virtual double density() const;

protected:
    double m_h0_RT; //!< Reference state enthalpy divided by RT
    double m_cp0_R; //!< Reference state heat capacity divided by R
    double m_s0_R; //!< Reference state entropy divided by R
    double m_g0_RT; //!< Reference state Gibbs free energy divided by RT
    double m_V0; //!< Reference state molar volume (m^3/kmol)
    double m_hss_RT; //!< Standard state enthalpy divided by RT
    double m_cpss_R; //!< Standard state heat capacity divided by R
    double m_sss_R; //!< Standard state entropy divided by R
    double m_gss_RT; //!< Standard state Gibbs free energy divided by RT
    double m_Vss; //!< Standard State molar volume (m^3/kmol)
};

}

#endif
