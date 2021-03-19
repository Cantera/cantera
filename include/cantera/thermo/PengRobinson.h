//! @file PengRobinson.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PENGROBINSON_H
#define CT_PENGROBINSON_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Peng-Robinson equation of state
 *
 * @ingroup thermoprops
 */
class PengRobinson : public MixtureFugacityTP
{
public:
    //! @name Constructors and Duplicators
    //! @{

    //! Base constructor.
    PengRobinson();

    //! Construct and initialize a PengRobinson object directly from an
    //! ASCII input file
    /*!
     * @param infile    Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the empty
     *     string.
     */
    PengRobinson(const std::string& infile, const std::string& id="");

    //! Construct and initialize a PengRobinson object directly from an
    //! XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id       id attribute containing the name of the phase.  (default
     *      is the empty string)
     */
    PengRobinson(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "PengRobinson";
    }

    //! @name Molar Thermodynamic properties
    //! @{

    virtual double cp_mole() const;
    virtual double cv_mole() const;

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  Since the mass density, temperature, and mass fractions are stored,
     *  this method uses these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots, Y_K) \f$.
     *
     * \f[
     *    P = \frac{RT}{v-b_{mix}} - \frac{\left(\alpha a\right)_{mix}}{v^2 + 2b_{mix}v - b_{mix}^2 }
     * \f]
     *
     *  where:
     *
     * \f[
     *    \alpha = \left[ 1 + \kappa \left(1-T_r^{0.5}\right)\right]^2
     * \f]
     *
     *  and
     *
     * \f[
     *    \kappa = \left(0.37464 + 1.54226\omega - 0.26992\omega^2\right)                              for omega <= 0.491
     *    \kappa = \left(0.379642 + 1.487503\omega - 0.164423\omega^2 +  0.016667\omega^3 \right)      for omega > 0.491
     * \f]
     *
     *Coefficients a_mix, b_mix and (a \alpha)_{mix} are caclulated as
     *
     *\f[
     *a_{mix} = \sum_i \sum_j X_i X_j a_{i, j} = \sum_i \sum_j X_i X_j sqrt{a_i a_j}
     *\f]
     *
     *\f[
     *   b_{mix} = \sum_i X_i b_i
     *\f]
     *
     *\f[
     *   {a \alpha}_{mix} = \sum_i \sum_j X_i X_j {a \alpha}_{i, j} = \sum_i \sum_j X_i X_j sqrt{a_i a_j} sqrt{\alpha_i \alpha_j}
     *\f]
     *
     *
     */
    
    virtual double pressure() const;

    // @}

protected:

    virtual void setTemperature(const double temp);
    virtual void compositionChanged();

public:

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase. 
     * The ideal gas mixture is considered as the standard or reference state here. 
     * Since the activity for an ideal gas mixture is simply the mole fraction, for an ideal gas
     * \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default is to
     *          assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual double standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(double* ac) const;

    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the array of non-dimensional species chemical potentials.
    //! These are partial molar Gibbs free energies.
    /*!
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function here calling getChemPotentials() and
     * then dividing by RT. No need for child classes to handle.
     *
     * @param mu    Output vector of non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    virtual void getChemPotentials_RT(double* mu) const;

    virtual void getChemPotentials(double* mu) const;
    virtual void getPartialMolarEnthalpies(double* hbar) const;
    virtual void getPartialMolarEntropies(double* sbar) const;
    virtual void getPartialMolarIntEnergies(double* ubar) const;
    virtual void getPartialMolarCp(double* cpbar) const;
    virtual void getPartialMolarVolumes(double* vbar) const;

    //! Calculate the temperature dependent interaction parameter alpha needed for P-R EoS
    /*
    *  The temperature dependent parameter in P-R EoS is calculated as
    *       \alpha = [1 + \kappa(1 - sqrt{T/T_crit}]^2
    *  kappa is a function calulated based on the accentric factor.
    *  Units: unitless
    */
    virtual void calculateAlpha(const std::string& species, double a, double b, double w);
    //@}
    /// @name Critical State Properties.
    //@{

    virtual double critTemperature() const;
    virtual double critPressure() const;
    virtual double critVolume() const;
    virtual double critCompressibility() const;
    virtual double critDensity() const;
    virtual double speciesCritTemperature(double a, double b) const;

public:
    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void setParametersFromXML(const XML_Node& thermoNode);
    virtual void setToEquilState(const double* lambda_RT);
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Retrieve a and b coefficients by looking up tabulated critical parameters
    /*!
    *  If pureFluidParameters are not provided for any species in the phase,
    *  consult the critical properties tabulated in build/data/thermo/critProperties.xml.
    *  If the species is found there, calculate pure fluid parameters a_k and b_k as:
    *  \f[ a_k = 0.4278*R**2*T_c^2/P_c \f]
    *
    *  and:
    *  \f[ b_k = 0.08664*R*T_c/P_c \f]
    *
    *  @param iName    Name of the species
    */
    virtual std::vector<double> getCoeff(const std::string& iName);

    //! Set the pure fluid interaction parameters for a species
    /*!
     *  The "a" parameter for species *i* in the Peng-Robinson model is assumed
     *  to be a linear function of temperature:
     *  \f[ a = a_0 + a_1 T \f]
     *
     *  @param species   Name of the species
     *  @param a0        constant term in the expression for the "a" parameter
     *      of the specified species [Pa-m^6/kmol^2]
     *  @param a1        temperature-proportional term in the expression for the
     *      "a" parameter of the specified species [Pa-m^6/kmol^2/K]
     *  @param b         "b" parameter in the Peng-Robinson model [m^3/kmol]
     *  @param alpha     dimensionless function of T_r and \omega
     *  @param omega     acentric factor
     */
    void setSpeciesCoeffs(const std::string& species, double a, double b,
                              double w);

    //! Set values for the interaction parameter between two species
    /*!
     *  The "a" parameter for interactions between species *i* and *j* is
     *  assumed by default to be computed as:
     *  \f[ a_{ij} = \sqrt(a_{i, 0} a_{j, 0}) + \sqrt(a_{i, 1} a_{j, 1}) T \f]
     *
     *  This function overrides the defaults with the specified parameters:
     *  \f[ a_{ij} = a_{ij, 0} + a_{ij, 1} T \f]
     *
     *  @param species_i   Name of one species
     *  @param species_j   Name of the other species
     *  @param a0          constant term in the "a" expression [Pa-m^6/kmol^2]
     *  @param a1          temperature-proportional term in the "a" expression
     *      [Pa-m^6/kmol^2/K]
     */
    void setBinaryCoeffs(const std::string& species_i,
                         const std::string& species_j, double a0, double a1);

private:
    //! Read the pure species PengRobinson input parameters
    /*!
     *  @param pureFluidParam   XML_Node for the pure fluid parameters
     */
    void readXMLPureFluid(XML_Node& pureFluidParam);

    //! Read the cross species PengRobinson input parameters
    /*!
     *  @param crossFluidParam   XML_Node for the cross fluid parameters
     */
    void readXMLCrossFluid(XML_Node& crossFluidParam);

    // @}

protected:
    // Special functions inherited from MixtureFugacityTP
    virtual double sresid() const;
    virtual double hresid() const;

public:
    virtual double liquidVolEst(double TKelvin, double& pres) const;
    virtual double densityCalc(double TKelvin, double pressure, int phase, double rhoguess);

    virtual double densSpinodalLiquid() const;
    virtual double densSpinodalGas() const;
    virtual double pressureCalc(double TKelvin, double molarVol) const;
    virtual double dpdVCalc(double TKelvin, double molarVol, double& presCalc) const;

    //! Calculate dpdV and dpdT at the current conditions
    /*!
     *  These are stored internally.
     */
    void pressureDerivatives() const;

    virtual void updateMixingExpressions();

    //! Update the a and b parameters
    /*!
     *  The a and the b parameters depend on the mole fraction and the
     *  temperature. This function updates the internal numbers based on the
     *  state of the object.
     */
    void updateAB();

    //! Calculate the a and the b parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param temp  Temperature (TKelvin)
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     */
    void calculateAB(double temp, double& aCalc, double& bCalc, double& aAlpha) const;

    // Special functions not inherited from MixtureFugacityTP

    double daAlpha_dT() const;
    double d2aAlpha_dT2() const;

    void calcCriticalConditions(double a, double b,double& pc, double& tc, double& vc) const;

    //! Prepare variables and call the function to solve the cubic equation of state
    int NicholsCall(double T, double pres, double a, double b, double aAlpha,
                     double Vroot[3]) const;
protected:
    //! Form of the temperature parameterization
    /*!
     *  0 = There is no temperature parameterization of a or b
     *  1 = The a_ij parameter is a linear function of the temperature
     */
    int m_formTempParam;

    //! Value of b in the equation of state
    /*!
     *  m_b_current is a function of the temperature and the mole fractions.
     */
    double m_b_current;

    //! Value of a and alpha in the equation of state
    /*!
     *  m_aAlpha_current is a function of the temperature and the mole fractions. m_a_current depends only on the mole fractions.
     */
    double m_a_current;
    double m_aAlpha_current;

    // Vectors required to store a_coeff, b_coeff, alpha, kappa and other values for every species. Length = m_kk
    vector_fp a_vec_Curr_;
    vector_fp b_vec_Curr_;
    vector_fp aAlpha_vec_Curr_;
    vector_fp alpha_vec_Curr_;
    vector_fp kappa_vec_;
    mutable vector_fp dalphadT_vec_Curr_;
    mutable vector_fp d2alphadT2_;

    Array2D a_coeff_vec;
    Array2D aAlpha_coeff_vec;
    
    int NSolns_;

    double Vroot_[3];

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_pp;

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_tmpV;

    // Partial molar volumes of the species
    mutable vector_fp m_partialMolarVolumes;

    //! The derivative of the pressure with respect to the volume
    /*!
     * Calculated at the current conditions. temperature and mole number kept
     * constant
     */
    mutable double dpdV_;

    //! The derivative of the pressure with respect to the temperature
    /*!
     *  Calculated at the current conditions. Total volume and mole number kept
     *  constant
     */
    mutable double dpdT_;

    //! Vector of derivatives of pressure with respect to mole number
    /*!
     *  Calculated at the current conditions. Total volume, temperature and
     *  other mole number kept constant
     */
    mutable vector_fp dpdni_;

public:
    //! Omega constants: a0 (= omega_a) and b0 (= omega_b) values used in Peng-Robinson equation of state
    /*!
     *  These values are calculated by solving P-R cubic equation at the critical point.
     */
    static const double omega_a;

    //! Omega constant for b
    static const double omega_b;

    //! Omega constant for the critical molar volume
    static const double omega_vc;
};
}

#endif
