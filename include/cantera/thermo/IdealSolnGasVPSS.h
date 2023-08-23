/**
 *  @file IdealSolnGasVPSS.h
 * Definition file for a derived class of ThermoPhase that assumes
 * an ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see @ref thermoprops and
 * class @link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS@endlink).
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
    //! Create an object from an input file
    explicit IdealSolnGasVPSS(const string& infile="", string id="");

    //! @name  Utilities (IdealSolnGasVPSS)
    //! @{

    string type() const override {
        return "ideal-solution-VPSS";
    }

    bool isIdeal() const override {
        return true;
    }

    //! Set the standard concentration model
    /**
     * Must be one of 'unity', 'species-molar-volume', or 'solvent-molar-volume'.
     */
    void setStandardConcentrationModel(const string& model);

    //! @}
    //! @name Molar Thermodynamic Properties
    //! @{

    double enthalpy_mole() const override;
    double entropy_mole() const override;
    double cp_mole() const override;
    double cv_mole() const override;

    //! @}
    //! @name Mechanical Properties
    //! @{

    void setPressure(double p) override;

protected:
    void calcDensity() override;
    //! @}

public:
    Units standardConcentrationUnits() const override;
    void getActivityConcentrations(double* c) const override;

    //! Returns the standard concentration @f$ C^0_k @f$, which is used to
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
    double standardConcentration(size_t k=0) const override;

    void getActivityCoefficients(double* ac) const override;


    //! @name  Partial Molar Properties of the Solution
    //! @{

    void getChemPotentials(double* mu) const override;
    void getPartialMolarEnthalpies(double* hbar) const override;
    void getPartialMolarEntropies(double* sbar) const override;
    void getPartialMolarIntEnergies(double* ubar) const override;
    void getPartialMolarCp(double* cpbar) const override;
    void getPartialMolarVolumes(double* vbar) const override;
    //! @}

public:
    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing the phase
    //! and setting its parameters from a specification in an input file. They
    //! are not normally used in application programs. To see how they are used,
    //! see importPhase().
    //! @{

    bool addSpecies(shared_ptr<Species> spec) override;
    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;
    void setToEquilState(const double* lambda_RT) override;

    //! @}

protected:
    //! form of the generalized concentrations
    /*!
     *    - 0 unity (default)
     *    - 1 1/V_k
     *    - 2 1/V_0
     */
    int m_formGC = 0;

    //! Temporary storage - length = m_kk.
    vector<double> m_pp;
};
}

#endif
