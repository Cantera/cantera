/**
 *  @file SurfPhase.h
 *  Header for a simple thermodynamics model of a surface phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model
 *  (see @ref thermoprops and class @link Cantera::SurfPhase SurfPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "ThermoPhase.h"

namespace Cantera
{

//! A simple thermodynamic model for a surface phase, assuming an ideal solution
//! model.
/*!
 * The surface consists of a grid of equivalent sites. Surface species may be
 * defined to occupy one or more sites. The surface species are assumed to be
 * independent, and thus the species form an ideal solution.
 *
 * The density of surface sites is given by the variable @f$ n_0 @f$,
 * which has SI units of kmol m-2.
 *
 * ## Specification of Species Standard State Properties
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
 * state thermodynamics is resolved at this level.
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * Therefore, The standard state internal energy for species *k* is equal to the
 * enthalpy for species *k*.
 *
 * @f[
 *      u^o_k = h^o_k
 * @f]
 *
 * Also, the standard state chemical potentials, entropy, and heat capacities
 * are independent of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The activity of species defined in the phase is given by
 * @f[
 *      a_k = \theta_k
 * @f]
 *
 * The chemical potential for species *k* is equal to
 * @f[
 *      \mu_k(T,P) = \mu^o_k(T) + R T \ln \theta_k
 * @f]
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * The internal energy for species k is equal to the enthalpy for species *k*
 * @f[
 *      u_k = h_k
 * @f]
 *
 * The entropy for the phase is given by the following relation, which is
 * independent of the pressure:
 *
 * @f[
 *      s_k(T,P) = s^o_k(T) - R \ln \theta_k
 * @f]
 *
 * ## Application within Kinetics Managers
 *
 * The activity concentration,@f$  C^a_k @f$, used by the kinetics manager, is equal to
 * the actual concentration, @f$ C^s_k @f$, and is given by the following
 * expression.
 * @f[
 *      C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}
 * @f]
 *
 * The standard concentration for species *k* is:
 * @f[
 *      C^0_k = \frac{n_0}{s_k}
 * @f]
 *
 * An example phase definition is given in the
 * <a href="../../sphinx/html/yaml/phases.html#ideal-surface"> YAML API Reference</a>.
 *
 * @ingroup thermoprops
 */
class SurfPhase : public ThermoPhase
{
public:
    //! Construct and initialize a SurfPhase ThermoPhase object directly from an
    //! input file
    /*!
     * @param infile name of the input file. If blank, an empty phase will be created.
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    explicit SurfPhase(const string& infile="", const string& id="");

    string type() const override {
        return "ideal-surface";
    }

    bool isCompressible() const override {
        return false;
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     * @f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T),
     * @f]
     * and is a function only of temperature. The standard-state pure-species
     * Enthalpies @f$ \hat h^0_k(T) @f$ are computed by the species
     * thermodynamic property manager.
     *
     * \see MultiSpeciesThermo
     */
    double enthalpy_mole() const override;

    //! Return the Molar Internal Energy. Units: J/kmol
    /**
     * For a surface phase, the pressure is not a relevant thermodynamic
     * variable, and so the Enthalpy is equal to the Internal Energy.
     */
    double intEnergy_mole() const override;

    //! Return the Molar Entropy. Units: J/kmol-K
    /**
     * @f[
     *  \hat s(T,P) = \sum_k X_k (\hat s^0_k(T) - R \ln \theta_k)
     * @f]
     */
    double entropy_mole() const override;

    double cp_mole() const override;
    double cv_mole() const override;

    void getChemPotentials(double* mu) const override;
    void getPartialMolarEnthalpies(double* hbar) const override;
    void getPartialMolarEntropies(double* sbar) const override;
    void getPartialMolarCp(double* cpbar) const override;
    void getPartialMolarVolumes(double* vbar) const override;
    void getStandardChemPotentials(double* mu0) const override;

    //! Return a vector of activity concentrations for each species
    /*!
     * For this phase the activity concentrations,@f$ C^a_k @f$, are defined to
     * be equal to the actual concentrations, @f$ C^s_k @f$. Activity
     * concentrations are
     *
     * @f[
     *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}
     * @f]
     *
     * where @f$ \theta_k @f$ is the surface site fraction for species k,
     * @f$ n_0 @f$ is the surface site density for the phase, and
     * @f$ s_k @f$ is the surface size of species k.
     *
     * @f$ C^a_k @f$ that are defined such that @f$ a_k = C^a_k / C^0_k, @f$
     * where @f$ C^0_k @f$ is a standard concentration defined below and @f$ a_k
     * @f$ are activities used in the thermodynamic functions.  These activity
     * concentrations are used by kinetics manager classes to compute the
     * forward and reverse rates of elementary reactions. Note that they may or
     * may not have units of concentration --- they might be partial pressures,
     * mole fractions, or surface coverages,
     *
     * @param c vector of activity concentration (kmol m-2).
     */
    void getActivityConcentrations(double* c) const override;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration @f$ C^0_k @f$ used to normalize the activity
     * (that is, generalized) concentration. For this phase, the standard
     * concentration is species- specific
     *
     * @f[
     *            C^0_k = \frac{n_0}{s_k}
     * @f]
     *
     * This definition implies that the activity is equal to @f$ \theta_k @f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return the standard concentration in units of kmol/m^2 for surface phases or
     *     kmol/m for edge phases.
     */
    double standardConcentration(size_t k=0) const override;
    double logStandardConc(size_t k=0) const override;

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;

    bool addSpecies(shared_ptr<Species> spec) override;

    //! Since interface phases have no volume, this returns 0.0.
    double molarVolume() const override {
        return 0.0;
    }

    //! Returns the site density
    /*!
     * Site density kmol m-2
     */
    double siteDensity() const {
        return m_n0;
    }

    //! Returns the number of sites occupied by one molecule of species *k*.
    double size(size_t k) const {
        return m_speciesSize[k];
    }

    //! Set the site density of the surface phase (kmol m-2)
    /*!
     *  @param n0 Site density of the surface phase (kmol m-2)
     */
    void setSiteDensity(double n0);

    void getGibbs_RT(double* grt) const override;
    void getEnthalpy_RT(double* hrt) const override;
    void getEntropy_R(double* sr) const override;
    void getCp_R(double* cpr) const override;
    void getStandardVolumes(double* vol) const override;

    //! Return the thermodynamic pressure (Pa).
    double pressure() const override {
        return m_press;
    }

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     *  @param p input Pressure (Pa)
     */
    void setPressure(double p) override {
        m_press = p;
    }

    void getPureGibbs(double* g) const override;
    void getGibbs_RT_ref(double* grt) const override;
    void getEnthalpy_RT_ref(double* hrt) const override;
    void getEntropy_R_ref(double* er) const override;
    void getCp_R_ref(double* cprt) const override;

    //! Set the surface site fractions to a specified state.
    /*!
     * This routine converts to concentrations in kmol/m2, using m_n0, the
     * surface site density, and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule. It then calls
     * Phase::setConcentrations to set the internal concentration in the object.
     *
     * @param theta    This is the surface site fraction for the kth species in
     *                 the surface phase. This is a dimensionless quantity.
     *
     * This routine normalizes the theta's to 1, before application
     */
    void setCoverages(const double* theta);

    //! Set the surface site fractions to a specified state.
    /*!
     * This routine converts to concentrations in kmol/m2, using m_n0, the
     * surface site density, and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule. It then calls
     * Phase::setConcentrations to set the internal concentration in the object.
     *
     * @param theta    This is the surface site fraction for the kth species in
     *                 the surface phase. This is a dimensionless quantity.
     */
    void setCoveragesNoNorm(const double* theta);

    //! Set the coverages from a string of colon-separated name:value pairs.
    /*!
     *  @param cov  String containing colon-separated name:value pairs
     */
    void setCoveragesByName(const string& cov);

    //! Set the coverages from a map of name:value pairs
    void setCoveragesByName(const Composition& cov);

    //! Return a vector of surface coverages
    /*!
     * Get the coverages.
     *
     * @param theta Array theta must be at least as long as the number of
     *              species.
     */
    void getCoverages(double* theta) const;

    //! @copydoc ThermoPhase::setState
    /*!
     * Additionally uses the key `coverages` to set the fractional coverages.
     */
    void setState(const AnyMap& state) override;

protected:
    void compositionChanged() override;

    //! Surface site density (kmol m-2)
    double m_n0 = 1.0;

    //! Vector of species sizes (number of sites occupied). length m_kk.
    vector<double> m_speciesSize;

    //! log of the surface site density
    double m_logn0;

    //! Current value of the pressure (Pa)
    double m_press = OneAtm;

    //! Temporary storage for the reference state enthalpies
    mutable vector<double> m_h0;

    //! Temporary storage for the reference state entropies
    mutable vector<double> m_s0;

    //! Temporary storage for the reference state heat capacities
    mutable vector<double> m_cp0;

    //! Temporary storage for the reference state Gibbs energies
    mutable vector<double> m_mu0;

    //! Temporary work array
    mutable vector<double> m_work;

    //! vector storing the log of the size of each species.
    /*!
     * The size of each species is defined as the number of surface sites each
     * species occupies.
     */
    mutable vector<double> m_logsize;

    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only reevaluated if
     * the temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevaluation of the thermo
     *               polynomials. default = false.
     */
    void _updateThermo(bool force=false) const;
};

}

#endif
