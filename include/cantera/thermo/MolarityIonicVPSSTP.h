/**
 *  @file MolarityIonicVPSSTP.h (see \ref thermoprops and class \link
 *      Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles variable pressure
 * standard state methods for calculating thermodynamic properties that are
 * further based upon activities based on the molarity scale.  In this class, we
 * expect that there are ions, but they are treated on the molarity scale.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLARITYIONICVPSSTP_H
#define CT_MOLARITYIONICVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

/*!
 * MolarityIonicVPSSTP is a derived class of GibbsExcessVPSSTP that handles
 * variable pressure standard state methods for calculating thermodynamic
 * properties that are further based on expressing the Excess Gibbs free energy
 * as a function of the mole fractions (or pseudo mole fractions) of the
 * constituents. This category is the workhorse for describing ionic systems
 * which are not on the molality scale.
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * @deprecated To be removed after Cantera 2.4
 *
 * This class adds additional functions onto the ThermoPhase interface that
 * handles the calculation of the excess Gibbs free energy. The ThermoPhase
 * class includes a member function, ThermoPhase::activityConvention() that
 * indicates which convention the activities are based on. The default is to
 * assume activities are based on the molar convention. That default is used
 * here.
 *
 * All of the Excess Gibbs free energy formulations in this area employ
 * symmetrical formulations.
 *
 * This layer will massage the mole fraction vector to implement cation and
 * anion based mole numbers in an optional manner, such that it is expected that
 * there exists a charge balance at all times. One of the ions must be a
 * "special ion" in the sense that its' thermodynamic functions are set to zero,
 * and the thermo functions of all other ions are based on a valuation relative
 * to that special ion.
 */
class MolarityIonicVPSSTP : public GibbsExcessVPSSTP
{
public:
    MolarityIonicVPSSTP();

    //! Construct and initialize a MolarityIonicVPSSTP ThermoPhase object
    //! directly from an XML input file
    /*!
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    MolarityIonicVPSSTP(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a MolarityIonicVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    MolarityIonicVPSSTP(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "MolarityIonic";
    }

    /**
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and pressure.
     * @{
     */

    virtual void getLnActivityCoefficients(doublereal* lnac) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    virtual void getChemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the molality-based
     * activity coefficient wrt temperature
     *
     * \f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     * \f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     * \f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     * \f]
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Returns an array of partial molar entropies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     * \f[
     *   ???????????????
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *   ???????????????
     * \f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}

    //! Calculate pseudo binary mole fractions
    virtual void calcPseudoBinaryMoleFractions() const;

    /// @name Initialization
    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see importPhase().
    /// @{

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);
    //! @}

    virtual std::string report(bool show_thermo=true,
                               doublereal threshold=1e-14) const;

private:
    //! Initialize lengths of local variables after all species have been
    //! identified.
    void initLengths();

    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe the
     * Redlich-Kister model for a particular binary interaction. This function
     * reads the XML file and writes the coefficients it finds to an internal
     * data structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named
     *     "binaryNeutralSpeciesParameters" containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node& xmlBinarySpecies);

    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally stored natural
     * logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients wrt T
    /*!
     * This function will be called to update the internally stored derivative
     * of the natural logarithm of the activity coefficients wrt temperature.
     */
    void s_update_dlnActCoeff_dT() const;

    //! Internal routine that calculates the derivative of the activity
    //! coefficients wrt the mole fractions.
    /*!
     * This routine calculates the the derivative of the activity coefficients
     * wrt to mole fraction with all other mole fractions held constant. This is
     * strictly not permitted. However, if the resulting matrix is multiplied by
     * a permissible deltaX vector then everything is ok.
     *
     * This is the natural way to handle concentration derivatives in this
     * routine.
     */
    void s_update_dlnActCoeff_dX_() const;

protected:
    // Pseudobinary type
    /*!
     *  - `PBTYPE_PASSTHROUGH` - All species are passthrough species
     *  - `PBTYPE_SINGLEANION` - there is only one anion in the mixture
     *  - `PBTYPE_SINGLECATION` - there is only one cation in the mixture
     *  - `PBTYPE_MULTICATIONANION` - Complex mixture
     */
    int PBType_;

    //! Number of pseudo binary species
    size_t numPBSpecies_;

    mutable vector_fp PBMoleFractions_;

    //! Vector of cation indices in the mixture
    std::vector<size_t> cationList_;

    std::vector<size_t> anionList_;

    std::vector<size_t> passThroughList_;
    size_t neutralPBindexStart;

    mutable vector_fp moleFractionsTmp_;
};

#define PBTYPE_PASSTHROUGH 0
#define PBTYPE_SINGLEANION 1
#define PBTYPE_SINGLECATION 2
#define PBTYPE_MULTICATIONANION 3

}

#endif
