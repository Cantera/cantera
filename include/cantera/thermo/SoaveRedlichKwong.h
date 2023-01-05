//! @file SoaveRedlichKwong.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOAVEREDLICHKWONG_H
#define CT_SOAVEREDLICHKWONG_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Soave-Redlich-Kwong equation of state
 *
 * @ingroup thermoprops
 */
class SoaveRedlichKwong : public MixtureFugacityTP
{
public:
    explicit SoaveRedlichKwong(const std::string& infile="",
                               const std::string& id="");

    virtual std::string type() const {
        return "Soave-Redlich-Kwong";
    }


    //! Calculate species-specific critical temperature
    /*!
     *  The temperature dependent parameter in SRK EoS is calculated as
     *       \f[ T_{crit} = (0.0778 a)/(0.4572 b R) \f] TODO
     *  Units: Kelvin
     *
     * @param a    species-specific coefficients used in P-R EoS
     * @param b    species-specific coefficients used in P-R EoS
     */
    double speciesCritTemperature(double a, double b) const;

    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void initThermo();
    virtual void getSpeciesParameters(const std::string& name,
                                      AnyMap& speciesNode) const;

protected:

    virtual double dpdVCalc(double T, double molarVol, double& presCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    //! Calculate temperature derivative \f$d(a \alpha)/dT\f$
    /*!
     *  These are stored internally.
     */
    double daAlpha_dT() const;

    //! Calculate second derivative \f$d^2(a \alpha)/dT^2\f$
    /*!
     *  These are stored internally.
     */
    double d2aAlpha_dT2() const;

public:
    //! Calculate \f$dp/dV\f$ and \f$dp/dT\f$ at the current conditions
    /*!
     *  These are stored internally.
     */
    void calculatePressureDerivatives() const;

    //! Prepare variables and call the function to solve the cubic equation of state
    int solveCubic(double T, double pres, double a, double b, double aAlpha,
                   double Vroot[3]) const;

protected:
    //! Value of \f$b\f$ in the equation of state
    /*!
     *  `m_b` is a function of the mole fractions and species-specific b values.
     */
    double m_b;

    //! Value of \f$a\f$ in the equation of state
    /*!
     *  `m_a` depends only on the mole fractions.
     */
    double m_a;

    //! Value of \f$\alpha\f$ in the equation of state
    /*!
     *  `m_aAlpha_mix` is a function of the temperature and the mole fractions.
     */
    double m_aAlpha_mix;

    //! The derivative of the pressure with respect to the volume
    /*!
     * Calculated at the current conditions. temperature and mole number kept
     * constant
     */
    mutable double m_dpdV;

    //! The derivative of the pressure with respect to the temperature
    /*!
     *  Calculated at the current conditions. Total volume and mole number kept
     *  constant
     */
    mutable double m_dpdT;


    enum class CoeffSource { EoS, CritProps, Database };
    //! For each species, specifies the source of the a, b, and omega coefficients
    std::vector<CoeffSource> m_coeffSource;

private:
    static const double omega_a;
    static const double omega_b;
    static const double omega_vc;

};
}

#endif
