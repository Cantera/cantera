/**
 *  @file LTPspecies.h
 *  Header file defining class LTPspecies and its child classes
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_LTPSPECIES_H
#define CT_LTPSPECIES_H

#include "TransportBase.h"

namespace Cantera
{

//! Enumeration of the types of transport properties that can be
//! handled by the variables in the various Transport classes.
/*!
 * Not all of these are handled by each class and each class
 * should handle exceptions where the transport property is not handled.
 *
 * Transport properties currently on the list
 *
 *    0  - viscosity
 *    1  - Ionic conductivity
 *    2  - Mobility Ratio
 *    3  - Self Diffusion coefficient
 *    4  - Thermal conductivity
 *    5  - species diffusivity
 *    6  - hydrodynamic radius
 *    7  - electrical conductivity
 */
enum TransportPropertyType {
    TP_UNKNOWN = -1,
    TP_VISCOSITY = 0,
    TP_IONCONDUCTIVITY,
    TP_MOBILITYRATIO,
    TP_SELFDIFFUSION,
    TP_THERMALCOND,
    TP_DIFFUSIVITY,
    TP_HYDRORADIUS,
    TP_ELECTCOND,
    TP_DEFECTCONC,
    TP_DEFECTDIFF
};

//! Temperature dependence type for standard state species properties
/*!
 *  Types of temperature dependencies:
 *     0  - Independent of temperature
 *     1  - extended arrhenius form
 *     2  - polynomial in temperature form
 *     3 -  exponential temperature polynomial
 */
enum LTPTemperatureDependenceType {
    LTP_TD_NOTSET=-1,
    LTP_TD_CONSTANT,
    LTP_TD_ARRHENIUS,
    LTP_TD_POLY,
    LTP_TD_EXPT
};

//! Class LTPspecies holds transport parameterizations for a specific liquid-
//! phase species.
/*!
 * @deprecated To be removed after Cantera 2.4
 *
 * Subclasses handle different means of specifying transport properties like
 * constant, %Arrhenius or polynomial temperature fits.  In its current state,
 * it is primarily suitable for specifying temperature dependence, but the
 * adjustCoeffsForComposition() method can be implemented to adjust for the
 * composition dependence.
 *
 * Mixing rules for computing mixture transport properties are handled
 * separately in the LiquidTranInteraction subclasses.
 */
class LTPspecies
{
public:
    LTPspecies();
    virtual ~LTPspecies() {}

    //! Duplication routine
    /*!
     *   (virtual from LTPspecies)
     *
     *  @returns a copy of this routine as a pointer to LTPspecies
     */
    virtual LTPspecies* duplMyselfAsLTPspecies() const;

    //! Set the ThermoPhase object, which is used to find the temperature
    void setThermo(thermo_t* thermo) {
        m_thermo = thermo;
    }

    //! Set the species name
    void setName(const std::string& name) {
        m_speciesName = name;
    }

    //! TransportPropertyType containing the property id that this object is
    //! creating a parameterization for (e.g., viscosity)
    void setTransportPropertyType(TransportPropertyType tp_ind) {
        m_property = tp_ind;
    }

    //! Set up the species transport property from the XML node, `<propNode>`
    //! that is a child of the `<transport>` node in the species block and
    //! specifies a type of transport property (like viscosity)
    virtual void setupFromXML(const XML_Node& propNode) {}

    //! Returns the vector of standard state species transport property
    /*!
     * The standard state species transport property is returned.  Any
     * temperature and composition dependence will be adjusted internally
     * according to the information provided by the subclass object.
     *
     * @returns a single double containing the property evaluation at the
     *           current ThermoPhase temperature.
     */
    virtual doublereal getSpeciesTransProp();

    //! Check to see if the property evaluation will be positive
    virtual bool checkPositive() const;

    //! Return the weight mixture
    doublereal getMixWeight() const;



private:
    //! Internal model to adjust species-specific properties for the composition.
    /*!
     * Currently just a place holder, but this method could take the composition
     * from the thermo object and adjust coefficients according to some yet
     * unspecified model.
     */
    virtual void adjustCoeffsForComposition();

protected:
    //! Species Name for the property that is being described
    std::string m_speciesName;

    //! Model type for the temperature dependence
    LTPTemperatureDependenceType m_model;

    //! enum indicating which property this is (i.e viscosity)
    TransportPropertyType m_property;

    //! Model temperature-dependence ceofficients
    vector_fp m_coeffs;

    //! Pointer to a const thermo object to get current temperature
    const thermo_t* m_thermo;

    //! Weighting used for mixing.
    /*!
     * This weighting can be employed to allow salt transport properties to be
     * represented by specific ions. For example, to have Li+ and Ca+ represent
     * the mixing transport properties of LiCl and CaCl2, the weightings for Li+
     * would be 2.0, for K+ would be 3.0 and for Cl- would be 0.0. The transport
     * properties for Li+ would be those for LiCl and the transport properties
     * for Ca+ would be those for CaCl2. The transport properties for Cl- should
     * be something innoccuous like 1.0--note that 0.0 is not innocuous if there
     * are logarithms involved.
     */
    doublereal m_mixWeight;
};

//! Class LTPspecies_Const holds transport parameters for a specific liquid-
//! phase species (LTPspecies) when the transport property is just a constant
//! value.
/*!
 * As an example of the input required for LTPspecies_Const consider the
 * following XML fragment
 *
 * \verbatim
 *    <species>
 *      <!-- thermodynamic properties -->
 *      <transport>
 *        <hydrodynamicRadius model="Constant" units="A">
 *            1.000
 *        </hydrodynamicRadius>
 *        <!-- other transport properties -->
 *      </transport>
 *    </species>
 * \endverbatim
 */
class LTPspecies_Const : public LTPspecies
{
public:
    LTPspecies_Const();

    virtual LTPspecies* duplMyselfAsLTPspecies() const;
    virtual void setupFromXML(const XML_Node& propNode);

    // Set the (constant) property value
    void setCoeff(double C);

    doublereal getSpeciesTransProp();
};

//! Class LTPspecies_Arrhenius holds transport parameters for a specific liquid-
//! phase species (LTPspecies) when the transport property is expressed in
//! Arrhenius form.
/*!
 * Used for standard state species properties with equations of the form
 * \f[
 *      x = A T^b \exp( - E / RT )
 * \f]
 * where A, b, and E are passed in the XML input file.
 *
 * As an example of the input required for LTPspecies_Arrhenius consider the
 * following XML fragment
 *
 * \verbatim
 *    <species>
 *      <!-- thermodynamic properties -->
 *      <transport>
 *        <viscosity model="Arrhenius">
 *           <!-- Janz, JPCRD, 17, supplement 2, 1988 -->
 *           <A>6.578e-5</A>
 *           <b>0.0</b>
 *           <E units="J/kmol">23788.e3</E>
 *        </viscosity>
 *        <!-- other transport properties -->
 *      </transport>
 *    </species>
 * \endverbatim
 */
class LTPspecies_Arrhenius : public LTPspecies
{
public:
    LTPspecies_Arrhenius();
    virtual LTPspecies* duplMyselfAsLTPspecies() const;
    virtual void setupFromXML(const XML_Node& propNode);

    //! Return the standard state species value for this transport property
    //! evaluated from the Arrhenius expression
    /*!
     * In general the Arrhenius expression is
     *
     * \f[
     *      \mu = A T^n \exp( - E / R T ).
     * \f]
     *
     * Note that for viscosity, the convention is such that a positive
     * activation energy corresponds to the typical case of a positive
     * argument to the exponential so that the Arrhenius expression is
     *
     * \f[
     *      \mu = A T^n \exp( + E / R T ).
     * \f]
     *
     * Any temperature and composition dependence will be adjusted internally
     * according to the information provided.
     */
    doublereal getSpeciesTransProp();

    // Set the coefficients in the Arrhenius expression
    void setCoeffs(double A, double n, double Tact);

protected:
    //! temperature from thermo object
    doublereal m_temp;

    //! logarithm of current temperature
    doublereal m_logt;

    //! most recent evaluation of transport property
    doublereal m_prop;

    //! logarithm of most recent evaluation of transport property
    doublereal m_logProp;
};

//! Class LTPspecies_Poly holds transport parameters for a specific liquid-phase
//! species (LTPspecies) when the transport property is expressed as a
//! polynomial in temperature.
/*!
 * Used for standard state species properties with equations of the form
 * \f[
 *    x = f[0] + f[1] T + ... + f[N] T^N
 * \f]
 * where f[i] are elements of the float array passed in.
 *
 * As an example of the input required for LTPspecies_Poly consider the
 * following XML fragment
 *
 * \verbatim
 *    <species>
 *      <!-- thermodynamic properties -->
 *      <transport>
 *        <thermalConductivity model="coeffs">
 *           <floatArray size="2">  0.6, -15.0e-5 </floatArray>
 *        </thermalConductivity>
 *        <!-- other transport properties -->
 *      </transport>
 *    </species>
 * \endverbatim
 */
class LTPspecies_Poly : public LTPspecies
{
public:
    LTPspecies_Poly();

    virtual LTPspecies* duplMyselfAsLTPspecies() const;
    virtual void setupFromXML(const XML_Node& propNode);
    doublereal getSpeciesTransProp();

    // Set the coefficients in the polynomial expression
    void setCoeffs(size_t N, const double* coeffs);

protected:
    //! temperature from thermo object
    doublereal m_temp;

    //! most recent evaluation of transport property
    doublereal m_prop;
};

//! Class LTPspecies_ExpT holds transport parameters for a specific liquid-
//! phase species (LTPspecies) when the transport property is expressed as an
//! exponential in temperature.
/*!
 * Used for standard state species properties with equations of the form
 *
 * \f[
 *    x = f[0] \exp( f[1] T + ... + f[N] T^{N} )
 * \f]
 *
 * where f[i] are elements of the float array passed in.
 *
 * As an example of the input required for LTPspecies_ExpT consider the
 * following XML fragment
 *
 * \verbatim
 *    <species>
 *      <!-- thermodynamic properties -->
 *      <transport>
 *        <thermalConductivity model="expT">
 *           <floatArray size="2">  0.6, -15.0e-5 </floatArray>
 *        </thermalConductivity>
 *        <!-- other transport properties -->
 *      </transport>
 *    </species>
 * \endverbatim
 */
class LTPspecies_ExpT : public LTPspecies
{
public:
    LTPspecies_ExpT();

    virtual LTPspecies* duplMyselfAsLTPspecies() const;
    virtual void setupFromXML(const XML_Node& propNode);

    // Set the coefficients in the polynomial expression
    void setCoeffs(size_t N, const double* coeffs);

    doublereal getSpeciesTransProp();

protected:
    //! temperature from thermo object
    doublereal m_temp;

    //! most recent evaluation of the transport property
    doublereal m_prop;
};

}
#endif
