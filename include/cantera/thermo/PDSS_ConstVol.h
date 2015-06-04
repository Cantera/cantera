/**
 *  @file PDSS_ConstVol.h
 *    Declarations for the class PDSS_ConstVol (pressure dependent standard state)
 *    which handles calculations for a single species with a constant molar volume in a phase
 *    (see class \ref pdssthermo and \link Cantera::PDSS_ConstVol PDSS_ConstVol\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_PDSS_CONSTVOL_H
#define CT_PDSS_CONSTVOL_H

#include "PDSS.h"

namespace Cantera
{
//! Class for pressure dependent standard states that use a constant volume model
/*!
 * @ingroup pdssthermo
 */
class PDSS_ConstVol : public PDSS_Nondimensional
{
public:
    //! @name  Constructors
    //! @{

    //! Constructor
    /*!
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     */
    PDSS_ConstVol(VPStandardStateTP* tp, size_t spindex);

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     *  This function calls the constructPDSSFile member function.
     *
     *  @param tp        Pointer to the ThermoPhase object pertaining to the phase
     *  @param spindex   Species index of the species in the phase
     *  @param inputFile String name of the input file
     *  @param id        String name of the phase in the input file. The default
     *                   is the empty string, in which case the first phase in the
     *                   file is used.
     */
    PDSS_ConstVol(VPStandardStateTP* tp, size_t spindex,
                  const std::string& inputFile, const std::string& id = "");

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
    PDSS_ConstVol(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
                  const XML_Node& phaseRef, bool spInstalled);

    //! Copy Constructor
    /*!
     * @param b Object to be copied
     */
    PDSS_ConstVol(const PDSS_ConstVol& b);

    //! Assignment operator
    /*!
     * @param b Object to be copied
     */
    PDSS_ConstVol& operator=(const PDSS_ConstVol& b);

    virtual PDSS* duplMyselfAsPDSS() const;

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

    virtual void setPressure(doublereal pres);
    virtual void setTemperature(doublereal temp);
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //!  @name  Miscellaneous properties of the standard state
    //! @{

    virtual doublereal satPressure(doublereal t);

    //! @}
    //! @name Initialization of the Object
    //! @{

    virtual void initThermo();

    //! Initialization of a PDSS object using an
    //! input XML file.
    /*!
     * This routine is a precursor to constructPDSSXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param vptp_ptr    Pointer to the Variable pressure ThermoPhase object
     *                    This object must have already been malloced.
     * @param spindex     Species index within the phase
     * @param inputFile   XML file containing the description of the phase
     * @param id          Optional parameter identifying the name of the
     *                    phase. If none is given, the first XML
     *                    phase element will be used.
     */
    void constructPDSSFile(VPStandardStateTP* vptp_ptr, size_t spindex,
                           const std::string& inputFile, const std::string& id);

    //!  Initialization of a PDSS object using an XML tree
    /*!
     * This routine is a driver for the initialization of the object.
     *
     *   basic logic:
     *     - initThermo()                 (cascade)
     *     - getStuff from species Part of XML file
     *     - initThermoXML(phaseNode)      (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure ThermoPhase object
     *                   This object must have already been malloced.
     * @param spindex    Species index within the phase
     * @param speciesNode XML Node containing the species information
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     * @param spInstalled  Boolean indicating whether the species is
     *                     already installed.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& speciesNode,
                          const XML_Node& phaseNode, bool spInstalled);

    virtual void initThermoXML(const XML_Node& phaseNode, const std::string& id);

    //@}

private:
    //! Value of the constant molar volume for the species
    /*!
     *    m3 / kmol
     */
    doublereal m_constMolarVolume;
};

}

#endif
