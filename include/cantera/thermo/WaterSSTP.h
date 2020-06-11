/**
 *  @file WaterSSTP.h
 * Declares a ThermoPhase class consisting of pure water (see \ref thermoprops
 * and class \link Cantera::WaterSSTP WaterSSTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WATERSSTP_H
#define CT_WATERSSTP_H

#include "SingleSpeciesTP.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/thermo/WaterProps.h"

namespace Cantera
{

class WaterPropsIAPWS;
class WaterProps;
//! Class for single-component water. This is designed to cover just the liquid
//! part of water.
/*!
 * The reference is W. Wagner, A. Pruss, "The IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use," J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * ## Specification of Species Standard State Properties
 *
 * The offsets used in the steam tables are different than NIST's. They assume
 * u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the triple point conditions:
 *
 *      -  u(273.16, rho)    = 0.0
 *      -  s(273.16, rho)    = 0.0
 *      -  psat(273.16)      = 611.655 Pascal
 *      -  rho(273.16, psat) = 999.793 kg m-3
 *
 * These "steam table" assumptions are used by the WaterPropsIAPWS class.
 * Therefore, offsets must be calculated to make the thermodynamic properties
 * calculated within this class to be consistent with thermo properties within
 * Cantera.
 *
 * The thermodynamic base state for water is set to the NIST basis here by
 * specifying constants, #EW_Offset and #SW_Offset, one for energy quantities
 * and one for entropy quantities. The offsets are specified so that the
 * following properties hold:
 *
 * - Delta_Hfo_idealgas(298.15) = -241.826 kJ/gmol
 * - So_idealgas(298.15, 1bar)  =  188.835 J/gmolK
 *
 * (From http://webbook.nist.gov)
 *
 * The "o" here refers to a hypothetical ideal gas state. The way we achieve
 * this in practice is to evaluate at a very low pressure and then use the
 * theoretical ideal gas results to scale up to higher pressures:
 *
 * Ho(1bar) = H(P0)
 *
 * So(1bar) = S(P0) + RT ln(1bar/P0)
 *
 * ## %Application within Kinetics Managers
 *
 * This is unimplemented.
 *
 * ## Instantiation of the Class
 *
 * A new WaterSSTP object may be created by the following code snippets,
 * combined with an XML file given in the XML example section.
 *
 * @code
 *      ThermoPhase* w = newPhase("waterSSTPphase.xml");
 * @endcode
 *
 * or
 *
 * @code
 *      WaterSSTP *w = new WaterSSTP("waterSSTPphase.xml","");
 * @endcode
 *
 * or
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "waterSSTPphase.xml#water", 0);
 *    WaterSSTP *w = new WaterSSTP(*xm);
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "waterSSTPphase.xml#water", 0);
 *    WaterSSTP water;
 *    importPhase(*xm, &water);
 * @endcode
 *
 * ## XML Example
 *
 * An example of an XML Element named phase setting up a WaterSSTP object with
 * id "water" is given below.
 *
 * @code
 * <!-- phase water     -->
 * <phase dim="3" id="water">
 *   <elementArray datasrc="elements.xml">O  H </elementArray>
 *   <speciesArray datasrc="#species_data">H2O</speciesArray>
 *   <state>
 *     <temperature units="K">300.0</temperature>
 *     <pressure units="Pa">101325.0</pressure>
 *   </state>
 *   <thermo model="PureLiquidWater"/>
 *   <kinetics model="none"/>
 * </phase>
 * @endcode
 *
 *  Note the model "PureLiquidWater" indicates the usage of the WaterSSTP object.
 *
 * @ingroup thermoprops
 */
class WaterSSTP : public SingleSpeciesTP
{
public:
    //! Base constructor
    WaterSSTP();

    //! Full constructor for a water phase
    /*!
     * @param inputFile String name of the input file
     * @param id        string id of the phase name
     */
    explicit WaterSSTP(const std::string& inputFile, const std::string& id = "");

    //! Full constructor for a water phase
    /*!
     * @param phaseRef  XML node referencing the water phase.
     * @param id        string id of the phase name
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    explicit WaterSSTP(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "Water";
    }

    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    virtual doublereal cv_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties
    //@{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal p);
    virtual doublereal isothermalCompressibility() const;
    virtual doublereal thermalExpansionCoeff() const;

    //! Return the derivative of the volumetric thermal expansion coefficient.
    //! Units: 1/K2.
    virtual doublereal dthermalExpansionCoeffdT() const;

    //! @}
    //! @name Properties of the Standard State of the Species in the Solution
    //! @{

    virtual void getStandardChemPotentials(doublereal* gss) const;
    virtual void getGibbs_RT(doublereal* grt) const;
    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getCp_R(doublereal* cpr) const;
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //@}
    //! @name Thermodynamic Values for the Species Reference State
    /*!
     *  All functions in this group need to be overridden, because the
     *  m_spthermo MultiSpeciesThermo function is not adequate for the real
     *  equation of state.
     */
    //@{

    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const;
    virtual void getGibbs_ref(doublereal* g) const;
    virtual void getEntropy_R_ref(doublereal* er) const;
    virtual void getCp_R_ref(doublereal* cprt) const;
    virtual void getStandardVolumes_ref(doublereal* vol) const;
    //! @}

    virtual doublereal critTemperature() const;
    virtual doublereal critPressure() const;
    virtual doublereal critDensity() const;

    virtual doublereal satPressure(doublereal t);

    virtual bool compatibleWithMultiPhase() const {
        return false;
    }

    //! Return the fraction of vapor at the current conditions
    /*!
     * Below Tcrit, this routine will always return 0, by definition of the
     * functionality of the routine. Above Tcrit, we query the density to toggle
     * between 0 and 1.
     */
    virtual doublereal vaporFraction() const;

    //! Set the temperature of the phase
    /*!
     * The density and composition of the phase is constant during this
     * operator.
     *
     * @param temp Temperature (Kelvin)
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the density of the phase
    /*!
     * The temperature and composition of the phase is constant during this
     * operator.
     *
     * @param dens value of the density in kg m-3
     */
    virtual void setDensity(const doublereal dens);

    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return &m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return m_waterProps.get();
    }

protected:
    /**
     * @internal This internal routine must be overridden because it is not
     *        applicable.
     */
    void _updateThermo() const;

private:
    //! WaterPropsIAPWS that calculates the real properties of water.
    mutable WaterPropsIAPWS m_sub;

    //! Pointer to the WaterProps object
    /*!
     * This class is used to house several approximation routines for properties
     * of water. This object owns m_waterProps, and the WaterPropsIAPWS object
     * used by WaterProps is m_sub, which is defined above.
     */
    std::unique_ptr<WaterProps> m_waterProps;

    //! Molecular weight of Water -> Cantera assumption
    doublereal m_mw;

    //! Offset constants used to obtain consistency with the NIST database.
    /*!
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    doublereal EW_Offset;

    //! Offset constant used to obtain consistency with NIST convention.
    /*!
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    doublereal SW_Offset;

    //! Boolean is true if object has been properly initialized for calculation
    bool m_ready;

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
