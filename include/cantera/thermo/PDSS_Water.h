/**
 *  @file PDSS_Water.h
 * Implementation of a pressure dependent standard state
 * virtual function for a Pure Water Phase
 * (see \ref pdssthermo and class \link Cantera::PDSS_Water PDSS_Water\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_PDSS_WATER_H
#define CT_PDSS_WATER_H

#include "cantera/base/ct_defs.h"
#include "PDSS.h"
#include "VPStandardStateTP.h"

namespace Cantera
{
class WaterPropsIAPWS;
class WaterProps;

//!  Class for the liquid water pressure dependent
//!  standard state
/*!
 * Notes:
 *   Base state for thermodynamic properties:
 *
 *   The thermodynamic base state for water is set to the NIST basis here
 *   by specifying constants EW_Offset and SW_Offset. These offsets are
 *   specified so that the following properties hold:
 *
 *   Delta_Hfo_gas(298.15) = -241.826 kJ/gmol
 *   So_gas(298.15, 1bar)  = 188.835 J/gmolK
 *
 *   (http://webbook.nist.gov)
 *
 *   The "o" here refers to a hypothetical ideal gas state. The way
 *   we achieve this in practice is to evaluate at a very low pressure
 *   and then use the theoretical ideal gas results to scale up to
 *   higher pressures:
 *
 *   Ho(1bar) = H(P0)
 *
 *   So(1bar) = S(P0) + RT ln(1bar/P0)
 *
 *   The offsets used in the steam tables are different than NIST's.
 *   They assume u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the
 *   triple point conditions.
 *
 * @ingroup pdssthermo
 */
class PDSS_Water : public PDSS
{
public:
    //! @name  Constructors
    //! @{

    //! Bare constructor
    /*!
     *  eliminate?
     */
    PDSS_Water();

    //! Constructor that initializes the object by examining the XML entries
    //! from the ThermoPhase object
    /*!
     *  This function calls the constructPDSS member function.
     *
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS_Water(VPStandardStateTP* tp, int spindex);

    //! Copy Constructor
    /*!
     * @param b object to be copied
     */
    PDSS_Water(const PDSS_Water& b);

    //! Assignment operator
    /*!
     * @param b Object to be copied
     */
    PDSS_Water& operator=(const PDSS_Water& b);

    //! Constructor that initializes the object by examining the input file
    //! of the variable pressure ThermoPhase object
    /*!
     *  This function calls the constructPDSSFile member function.
     *
     *  @param tp        Pointer to the variable pressure ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param inputFile String name of the input file
     *  @param id        String name of the phase in the input file. The default
     *                   is the empty string, in which case the first phase in the
     *                   file is used.
     */
    PDSS_Water(VPStandardStateTP* tp, int spindex,
               const std::string& inputFile, const std::string& id = "");

    //! Constructor that initializes the object by examining the input file
    //! of the variable pressure ThermoPhase object
    /*!
     *  This function calls the constructPDSSXML member function.
     *
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param speciesNode Reference to the species XML tree.
     *  @param phaseRef  Reference to the XML tree containing the phase information.
     *  @param spInstalled Is the species already installed.
     */
    PDSS_Water(VPStandardStateTP* tp, int spindex, const XML_Node& speciesNode,
               const XML_Node& phaseRef, bool spInstalled);

    //! Destructor
    virtual ~PDSS_Water();

    //! Duplication routine for objects which inherit from %PDSS
    /*!
     *  This virtual routine can be used to duplicate %PDSS  objects
     *  inherited from %PDSS even if the application only has
     *  a pointer to %PDSS to work with.
     *
     * @return returns a pointer to the base %PDSS object type
     */
    virtual PDSS* duplMyselfAsPDSS() const;

    //! @}
    //! @name  Molar Thermodynamic Properties of the Species Standard State in the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_mole() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;
    virtual doublereal molarVolume() const;
    virtual doublereal density() const;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    //! Returns a reference pressure value that can be safely calculated by the
    //! underlying real equation of state for water
    /*!
     *  Note, this function is needed because trying to calculate a one atm
     *  value around the critical point will cause a crash
     *
     * @param temp  Temperature (Kelvin)
     */
    doublereal pref_safe(doublereal temp) const;

    virtual doublereal gibbs_RT_ref() const;
    virtual doublereal enthalpy_RT_ref() const;
    virtual doublereal entropy_R_ref() const;
    virtual doublereal cp_R_ref() const;
    virtual doublereal molarVolume_ref() const;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal pres);
    virtual void setTemperature(doublereal temp);
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! Set the density of the water phase
    /*!
     *  This is a non-virtual function because it specific
     *  to this object.
     *
     * @param dens Density of the water (kg/m3)
     */
    void setDensity(doublereal dens);

    virtual doublereal thermalExpansionCoeff() const;

    //! Return the derivative of the volumetric thermal expansion coefficient. Units: 1/K2.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal dthermalExpansionCoeffdT() const;

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  or
     * \f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    //! @}
    //! @name Miscellaneous properties of the standard state
    //! @{

    virtual doublereal critTemperature() const;
    virtual doublereal critPressure() const;
    virtual doublereal critDensity() const;
    virtual doublereal satPressure(doublereal t);

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return m_waterProps;
    }

    //! @}
    //! @name Initialization of the Object
    //! @{

    //! Internal routine that initializes the underlying water model
    void constructSet();

    //! Initialization of a PDSS object using an
    //! input XML file.
    /*!
     * This routine is a precursor to constructPDSSXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param vptp_ptr    Pointer to the Variable pressure %ThermoPhase object
     *                    This object must have already been malloced.
     *
     * @param spindex     Species index within the phase
     *
     * @param inputFile   XML file containing the description of the phase
     *
     * @param id          Optional parameter identifying the name of the
     *                    phase. If none is given, the first XML
     *                    phase element will be used.
     */
    void constructPDSSFile(VPStandardStateTP* vptp_ptr, int spindex,
                           const std::string& inputFile, const std::string& id);

    //!Initialization of a PDSS object using an xml tree
    /*!
     * This routine is a driver for the initialization of the
     * object.
     *
     *   basic logic:
     *     - initThermo()                 (cascade)
     *     - getStuff from species Part of XML file
     *     - initThermoXML(phaseNode)      (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure ThermoPhase object
     *                   This object must have already been malloced.
     *
     * @param spindex    Species index within the phase
     *
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     *
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, int spindex,
                          const XML_Node& phaseNode, const std::string& id);

    virtual void initThermo();
    virtual void initThermoXML(const XML_Node& phaseNode, const std::string& id);
    //@}

private:
    //! Pointer to the WaterPropsIAPWS object, which does the actual calculations
    //! for the real equation of state
    /*!
     * This object owns m_sub
     */
    mutable WaterPropsIAPWS* m_sub;

    //! Pointer to the WaterProps object
    /*!
     *   This class is used to house several approximation
     *   routines for properties of water.
     *
     * This object owns m_waterProps, and the WaterPropsIAPWS object used by
     * WaterProps is m_sub, which is defined above.
     */
    WaterProps* m_waterProps;

    //! State of the system - density
    /*!
     * Density is the independent variable here, but it's hidden behind the
     * object's interface.
     */
    doublereal m_dens;

    //! state of the fluid
    /*!
     *  @code
     *    0  WATER_GAS
     *    1  WATER_LIQUID
     *    2  WATER_SUPERCRIT
     *    3  WATER_UNSTABLELIQUID
     *    4  WATER_UNSTABLEGAS
     *  @endcode
     */
    int m_iState;

    /**
     *  Offset constants used to obtain consistency with the NIST database.
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    doublereal EW_Offset;

    /**
     *  Offset constant used to obtain consistency with NIST convention.
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    doublereal SW_Offset;

    //! Verbose flag - used?
    bool m_verbose;

public:
    /**
     *  Since this phase represents a liquid phase, it's an error to
     *  return a gas-phase answer. However, if the below is true, then
     *  a gas-phase answer is allowed. This is used to check the thermodynamic
     *  consistency with ideal-gas thermo functions for example.
     */
    bool m_allowGasPhase;
};

}

#endif
