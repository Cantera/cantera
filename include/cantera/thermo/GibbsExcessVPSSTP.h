/**
 *  @file GibbsExcessVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ Gibbs excess free energy based formulations
 *  (see \ref thermoprops
 * and class \link Cantera::GibbsExcessVPSSTP GibbsExcessVPSSTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GIBBSEXCESSVPSSTP_H
#define CT_GIBBSEXCESSVPSSTP_H

#include "VPStandardStateTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{

/*!
 * GibbsExcessVPSSTP is a derived class of ThermoPhase that handles variable
 * pressure standard state methods for calculating thermodynamic properties that
 * are further based on expressing the Excess Gibbs free energy as a function of
 * the mole fractions (or pseudo mole fractions) of constituents. This category
 * is the workhorse for describing molten salts, solid-phase mixtures of
 * semiconductors, and mixtures of miscible and semi-miscible compounds.
 *
 * It includes
 *   - regular solutions
 *   - Margules expansions
 *   - NTRL equation
 *   - Wilson's equation
 *   - UNIQUAC equation of state.
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
 * Chemical potentials of species k, \f$ \mu_o \f$, has the following general
 * format:
 *
 * \f[
 *    \mu_k = \mu^o_k(T,P) + R T ln( \gamma_k X_k )
 * \f]
 *
 * where \f$ \gamma_k^{\triangle} \f$ is a molar based activity coefficient for
 * species \f$k\f$.
 *
 * GibbsExcessVPSSTP contains an internal vector with the current mole fraction
 * vector. That's one of its primary usages. In order to keep the mole fraction
 * vector constant, all of the setState functions are redesigned at this layer.
 *
 * ### Activity Concentrations: Relationship of ThermoPhase to %Kinetics Expressions
 *
 * As explained in a similar discussion in the ThermoPhase class, the actual
 * units used in kinetics expressions must be specified in the ThermoPhase class
 * for the corresponding species. These units vary with the field of study.
 * %Cantera uses the concept of activity concentrations to represent this.
 * Activity concentrations are used directly in the expressions for kinetics.
 * Standard concentrations are used as the multiplicative constant that takes
 * the activity of a species and turns it into an activity concentration.
 * Standard concentrations must not depend on the concentration of the species
 * in the phase.
 *
 * Here we set a standard for the specification of the standard concentrations
 * for this class and all child classes underneath it. We specify here that the
 * standard concentration is equal to 1 for all species. Therefore, the
 * activities appear directly in kinetics expressions involving species in
 * underlying GibbsExcessVPSSTP phases.
 *
 * ### SetState Strategy
 *
 * All setState functions that set the internal state of the ThermoPhase object
 * are overloaded at this level, so that a current mole fraction vector is
 * maintained within the object.
 */
class GibbsExcessVPSSTP : public VPStandardStateTP
{
public:
    //! @name Constructors
    //! @{

    GibbsExcessVPSSTP() {}

    //! @}

    //! @}
    //! @name Mechanical Properties
    //! @{

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     *
     * NOTE: This is a non-virtual function, which is not a member of the
     *       ThermoPhase base class.
     */
    void calcDensity();

public:
    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T) + \hat R T \log a_k. \f] The
     * quantity \f$\mu_k^0(T,P)\f$ is the chemical potential at unit activity,
     * which depends only on temperature and pressure.
     * @{
     */

    virtual Units standardConcentrationUnits() const;
    virtual void getActivityConcentrations(doublereal* c) const;

    /**
     * The standard concentration \f$ C^0_k \f$ used to normalize the
     * generalized concentration. In many cases, this quantity will be the same
     * for all species in a phase - for example, for an ideal gas
     * \f$ C^0_k = P/\hat R T \f$. For this reason, this method returns a single
     * value, instead of an array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of different
     * sizes), this method may be called with an optional parameter indicating
     * the species.
     *
     * The standard concentration for defaulted to 1. In other words
     * the activity concentration is assumed to be 1.
     *
     * @param k species index. Defaults to zero.
     */
    virtual doublereal standardConcentration(size_t k=0) const;
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Get the array of non-dimensional activities (molality based for this
    //! class and classes that derive from it) at the current solution
    //! temperature, pressure, and solution concentration.
    /*!
     * \f[
     *  a_i^\triangle = \gamma_k^{\triangle} \frac{m_k}{m^\triangle}
     * \f]
     *
     * This function must be implemented in derived classes.
     *
     * @param ac     Output vector of molality-based activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const;

    virtual void getActivityCoefficients(doublereal* ac) const;

    //! Get the array of temperature derivatives of the log activity coefficients
    /*!
     * This function is virtual, and first appears in GibbsExcessVPSSTP.
     *
     * units = 1/Kelvin
     *
     * @param dlnActCoeffdT    Output vector of temperature derivatives of the
     *                         log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdT(doublereal* dlnActCoeffdT) const {
        throw NotImplementedError("GibbsExcessVPSSTP::getdlnActCoeffdT");
    }

    virtual void getdlnActCoeffdlnN(const size_t ld, doublereal* const dlnActCoeffdlnN) {
        throw NotImplementedError("GibbsExcessVPSSTP::getdlnActCoeffdlnN: "
                                  "nonzero and nonimplemented");
    }

    //! Get the array of log concentration-like derivatives of the log activity
    //! coefficients
    /*!
     * This function is a virtual method.  For ideal mixtures (unity activity
     * coefficients), this can return zero. Implementations should take the
     * derivative of the logarithm of the activity coefficient with respect to
     * the logarithm of the concentration-like variable (i.e. number of moles in
     * in a unit volume. ) that represents the standard state. This quantity is
     * to be used in conjunction with derivatives of that concentration-like
     * variable when the derivative of the chemical potential is taken.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnX    Output vector of derivatives of the
     *                         log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdlnX(doublereal* dlnActCoeffdlnX) const {
        throw NotImplementedError("GibbsExcessVPSSTP::getdlnActCoeffdlnX");
    }

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  Frequently, for this class of thermodynamics representations,
     *  the excess Volume due to mixing is zero. Here, we set it as
     *  a default. It may be overridden in derived classes.
     *
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    virtual const vector_fp& getPartialMolarVolumesVector() const;

    virtual bool addSpecies(shared_ptr<Species> spec);

protected:
    virtual void compositionChanged();

    //! utility routine to check mole fraction sum
    /*!
     * @param x   vector of mole fractions.
     */
    double checkMFSum(const doublereal* const x) const;

    //! Storage for the current values of the mole fractions of the species
    /*!
     * This vector is kept up-to-date when the setState functions are called.
     * Therefore, it may be considered to be an independent variable.
     *
     * Note in order to do this, the setState functions are redefined to always
     * keep this vector current.
     */
    mutable vector_fp moleFractions_;

    //! Storage for the current values of the activity coefficients of the
    //! species
    mutable vector_fp lnActCoeff_Scaled_;

    //! Storage for the current derivative values of the gradients with respect
    //! to temperature of the log of the activity coefficients of the species
    mutable vector_fp dlnActCoeffdT_Scaled_;

    //! Storage for the current derivative values of the gradients with respect
    //! to temperature of the log of the activity coefficients of the species
    mutable vector_fp d2lnActCoeffdT2_Scaled_;

    //! Storage for the current derivative values of the gradients with respect
    //! to logarithm of the mole fraction of the log of the activity
    //! coefficients of the species
    mutable vector_fp dlnActCoeffdlnN_diag_;

    //! Storage for the current derivative values of the gradients with respect
    //! to logarithm of the mole fraction of the log of the activity
    //! coefficients of the species
    mutable vector_fp dlnActCoeffdlnX_diag_;

    //! Storage for the current derivative values of the gradients with respect
    //! to logarithm of the species mole number of the log of the activity
    //! coefficients of the species
    /*!
     *  dlnActCoeffdlnN_(k, m)  is the derivative of ln(gamma_k) wrt ln mole
     *  number of species m
     */
    mutable Array2D dlnActCoeffdlnN_;
};

}

#endif
