/**
 *  @file PDSS_IonsFromNeutral.h
 *   Declarations for the class PDSS_IonsFromNeutral (
 *    which handles calculations for a single ion in a fluid, whose properties
 *    are calculated from another neutral molecule.
 *    (see \ref pdssthermo and class \link Cantera::PDSS_IonsFromNeutral PDSS_IonsFromNeutral\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_IONSFROMNEUTRAL_H
#define CT_PDSS_IONSFROMNEUTRAL_H

#include "PDSS.h"

namespace Cantera
{
class ThermoPhase;

//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * This class is for a single Ideal Gas species.
 *
 * @ingroup pdssthermo
 */
class PDSS_IonsFromNeutral : public PDSS_Nondimensional
{
public:
    //! @name  Constructors
    //! @{

    //! Constructor
    /*!
     * @param tp       Pointer to the ThermoPhase object pertaining to the phase
     * @param spindex  Species index of the species in the phase
     */
    PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex);

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
     * @deprecated To be removed after Cantera 2.3.
     */
    PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex,
                         const std::string& inputFile, const std::string& id = "");

    //! Constructor that initializes the object by examining the input file
    //! of the ThermoPhase object
    /*!
     * This function calls the constructPDSSXML member function.
     *
     * @param vptp_ptr    Pointer to the ThermoPhase object pertaining to the phase
     * @param spindex     Species index of the species in the phase
     * @param speciesNode Reference to the species XML tree.
     * @param phaseRef    Reference to the XML tree containing the phase information.
     * @param spInstalled Boolean indicating whether the species is installed
     *                    yet or not.
     */
    PDSS_IonsFromNeutral(VPStandardStateTP* vptp_ptr, size_t spindex, const XML_Node& speciesNode,
                         const XML_Node& phaseRef, bool spInstalled);

    PDSS_IonsFromNeutral(const PDSS_IonsFromNeutral& b);
    PDSS_IonsFromNeutral& operator=(const PDSS_IonsFromNeutral& b);
    virtual PDSS* duplMyselfAsPDSS() const;
    virtual void initAllPtrs(VPStandardStateTP* vptp_ptr, VPSSMgr* vpssmgr_ptr,
                             MultiSpeciesThermo* spthermo_ptr);

    //! @}
    //! @name  Molar Thermodynamic Properties of the Species Standard State in the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_RT() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_R() const;

    /*!
     * @copydoc PDSS::gibbs_RT()
     *
     * \f[
     *   \frac{\mu^o_k}{RT} = \sum_{m}{ \alpha_{m , k} \frac{\mu^o_{m}}{RT}} + ( 1 - \delta_{k,sp}) 2.0 \ln{2.0}
     * \f]
     *
     * *m* is the neutral molecule species index. \f$ \alpha_{m , k} \f$ is the
     * stoiciometric coefficient for the neutral molecule, *m*, that creates the
     * thermodynamics for the ionic species  *k*. A factor  \f$ 2.0 \ln{2.0} \f$
     * is added to all ions except for the species ionic species, which in this
     * case is the single anion species, with species index *sp*.
     */
    virtual doublereal gibbs_RT() const;
    virtual doublereal cp_R() const;
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

    virtual void setTemperature(doublereal temp);
    virtual doublereal temperature() const;
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //! @name Initialization of the Object
    //! @{

    //! Initialization of a PDSS object using an input XML file.
    /*!
     * This routine is a precursor to constructPDSSXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param vptp_ptr    Pointer to the Variable pressure ThermoPhase object
     * @param spindex     Species index within the phase
     * @param inputFile   XML file containing the description of the phase
     * @param id          Optional parameter identifying the name of the
     *                    phase. If none is given, the first XML
     *                    phase element will be used.
     * @deprecated To be removed after Cantera 2.3.
     */
    void constructPDSSFile(VPStandardStateTP* vptp_ptr, size_t spindex,
                           const std::string& inputFile, const std::string& id);

    //! Initialization of a PDSS object using an XML tree
    /*!
     * This routine is a driver for the initialization of the object.
     *
     *   basic logic:
     *     - initThermo()                 (cascade)
     *     - getStuff from species Part of XML file
     *     - initThermoXML(phaseNode)      (cascade)
     *
     * @param vptp_ptr   Pointer to the Variable pressure ThermoPhase object
     * @param spindex    Species index within the phase
     * @param speciesNode  Reference to the phase Information for the species
     *                     that this standard state refers to
     * @param phaseNode  Reference to the phase Information for the phase
     *                   that owns this species.
     * @param id         Optional parameter identifying the name of the
     *                   phase. If none is given, the first XML
     *                   phase element will be used.
     */
    void constructPDSSXML(VPStandardStateTP* vptp_ptr, size_t spindex,
                          const XML_Node& speciesNode,
                          const XML_Node& phaseNode, const std::string& id);

    virtual void initThermo();
    //@}

protected:
    //! Pointer to the Neutral Molecule ThermoPhase object
    /*!
     *  This is a shallow pointer.
     */
    const ThermoPhase* neutralMoleculePhase_;

public:
    //! Number of neutral molecule species that make up the stoichiometric
    //! vector for this species, in terms of calculating thermodynamic functions
    size_t numMult_;

    //! Vector of species indices in the neutral molecule ThermoPhase
    std::vector<size_t> idNeutralMoleculeVec;

    //! Stoichiometric coefficient for this species using the Neutral Molecule
    //! Species in the vector idNeutralMoleculeVec
    vector_fp factorVec;

    //! Add 2RTln2 to the entropy and Gibbs free energies for this species
    /*!
     *  This is true if this species is not the special species
     */
    bool add2RTln2_;

    //! Vector of length equal to the number of species in the neutral molecule phase
    mutable vector_fp tmpNM;

    //! True if this species is the special species
    int specialSpecies_;
};
}

#endif
