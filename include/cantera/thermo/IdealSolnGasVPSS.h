/**
 *  @file IdealSolnGasVPSS.h
 * Definition file for a derived class of ThermoPhase that assumes
 * an ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALSOLNGASVPSS_H
#define CT_IDEALSOLNGASVPSS_H

#include "VPStandardStateTP.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 * An ideal solution approximation of a phase. Uses variable
 * pressure standard state methods for calculating thermodynamic properties.
 */
class IdealSolnGasVPSS : public VPStandardStateTP
{
public:
    /// Create an object from an input file
    explicit IdealSolnGasVPSS(const std::string& infile="", std::string id="");

    //! @name  Utilities (IdealSolnGasVPSS)
    //! @{

    virtual std::string type() const {
        return "ideal-solution-VPSS";
    }

    virtual bool isIdeal() const {
        return true;
    }

    //! Set the standard concentration model
    /*
     * Must be one of 'unity', 'species-molar-volume', or 'solvent-molar-volume'.
     */
    void setStandardConcentrationModel(const std::string& model);

    //! @}
    //! @name Molar Thermodynamic Properties
    //! @{

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Properties
    //! @{

    void setPressure(doublereal p);

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input. The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the species standard state molar
     * volumes. The species molar volumes may be functions of temperature and
     * pressure.
     */
    virtual void calcDensity();
    //! @}

public:
    virtual Units standardConcentrationUnits() const;
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     *  For ideal gases, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;


    //! @name  Partial Molar Properties of the Solution
    //! @{

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;
    virtual void getPartialMolarCp(doublereal* cpbar) const;
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    //! @}

public:
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing the phase
     * and setting its parameters from a specification in an input file. They
     * are not normally used in application programs. To see how they are used,
     * see importPhase().
     */
    //! @{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void initThermo();
    virtual void getParameters(AnyMap& phaseNode) const;
    virtual void setToEquilState(const doublereal* lambda_RT);
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! @}

protected:
    //! form of the generalized concentrations
    /*!
     *    - 0 unity (default)
     *    - 1 1/V_k
     *    - 2 1/V_0
     */
    int m_formGC;

    //! Temporary storage - length = m_kk.
    vector_fp m_pp;
};
}

#endif
