/**
 *  @file PDSS_IdealGas.h
 *   Declarations for the class PDSS_IdealGas (pressure dependent standard state)
 *    which handles calculations for a single ideal gas species in a phase
 *    (see \ref pdssthermo and class \link Cantera::PDSS_IdealGas PDSS_IdealGas\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_IDEALGAS_H
#define CT_PDSS_IDEALGAS_H

#include "PDSS.h"

namespace Cantera
{
//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * This class is for a single Ideal Gas species.
 *
 * @ingroup pdssthermo
 */
class PDSS_IdealGas : public PDSS_Nondimensional
{
public:
    //! @name  Constructors
    //! @{

    //! Constructor
    /*!
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS_IdealGas(VPStandardStateTP* tp, int spindex);

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSXML member function.
     *
     *  @param vptp_ptr    Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex     Species index of the species in the phase
     *  @param speciesNode Reference to the species XML tree.
     *  @param phaseRef    Reference to the XML tree containing the phase information.
     *  @param spInstalled Boolean indicating whether the species is installed yet
     *                     or not.
     */
    PDSS_IdealGas(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
                  const XML_Node& phaseRef, bool spInstalled);

    //! @}
    //! @name Molar Thermodynamic Properties of the Species Standard State in the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_RT() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_R() const;
    virtual doublereal gibbs_RT() const;
    virtual doublereal cp_R() const;
    virtual doublereal cv_mole() const;
    virtual doublereal molarVolume() const;
    virtual doublereal density() const;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

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
    virtual doublereal temperature() const;
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //! @name Initialization of the Object
    //! @{

    //!Initialization of a PDSS object using an XML tree
    /*!
     * This routine is a driver for the initialization of the object.
     *
     * basic logic:
     *   - initThermo()                 (cascade)
     *   - get stuff from species part of XML file
     *   - initThermoXML(phaseNode)     (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure ThermoPhase object
     * @param spindex    Species index within the phase
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& phaseNode, const std::string& id);

    virtual void initThermo();
    //@}
};
}

#endif
