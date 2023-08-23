/**
 *  @file LatticeSolidPhase.h Header for a simple thermodynamics model of a bulk
 *      solid phase derived from ThermoPhase, assuming an ideal solution model
 *      based on a lattice of solid atoms (see @ref thermoprops and class @link
 *      Cantera::LatticeSolidPhase LatticeSolidPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LATTICESOLID_H
#define CT_LATTICESOLID_H

#include "ThermoPhase.h"

namespace Cantera
{

//! A phase that is comprised of a fixed additive combination of other lattice
//! phases
/*!
 * This is the main way %Cantera describes semiconductors and other solid
 * phases. This ThermoPhase object calculates its properties as a sum over other
 * LatticePhase objects. Each of the LatticePhase objects is a ThermoPhase
 * object by itself.
 *
 * The results from this LatticeSolidPhase model reduces to the LatticePhase
 * model when there is one lattice phase and the molar densities of the
 * sublattice and the molar density within the LatticeSolidPhase have the same
 * values.
 *
 * The mole fraction vector is redefined within the LatticeSolidPhase object.
 * Each of the mole fractions sum to one on each of the sublattices.  The
 * routine getMoleFraction() and setMoleFraction() have been redefined to use
 * this convention.
 *
 * ## Specification of Species Standard State Properties
 *
 * The standard state properties are calculated in the normal way for each of
 * the sublattices. The normal way here means that a thermodynamic polynomial in
 * temperature is developed. Also, a constant volume approximation for the
 * pressure dependence is assumed.  All of these properties are on a Joules per
 * kmol of sublattice constituent basis.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The sum over the LatticePhase objects is carried out by weighting each
 * LatticePhase object value with the molar density (kmol m-3) of its
 * LatticePhase. Then the resulting quantity is divided by the molar density of
 * the total compound. The LatticeSolidPhase object therefore only contains a
 * listing of the number of LatticePhase object that comprises the solid, and it
 * contains a value for the molar density of the entire mixture. This is the
 * same thing as saying that
 *
 * @f[
 *     L_i = L^{solid}   \theta_i
 * @f]
 *
 * @f$ L_i @f$ is the molar volume of the ith lattice. @f$ L^{solid} @f$ is the
 * molar volume of the entire solid. @f$  \theta_i @f$ is a fixed weighting
 * factor for the ith lattice representing the lattice stoichiometric
 * coefficient. For this object the @f$ \theta_i @f$ values are fixed.
 *
 * Let's take FeS2 as an example, which may be thought of as a combination of
 * two lattices: Fe and S lattice. The Fe sublattice has a molar density of 1
 * gmol cm-3. The S sublattice has a molar density of 2 gmol cm-3. We then
 * define the LatticeSolidPhase object as having a nominal composition of FeS2,
 * and having a molar density of 1 gmol cm-3.  All quantities pertaining to the
 * FeS2 compound will be have weights associated with the sublattices. The Fe
 * sublattice will have a weight of 1.0 associated with it. The S sublattice
 * will have a weight of 2.0 associated with it.
 *
 * ### Specification of Solution Density Properties
 *
 * Currently, molar density is not a constant within the object, even though the
 * species molar volumes are a constant.  The basic idea is that a swelling of
 * one of the sublattices will result in a swelling of of all of the lattices.
 * Therefore, the molar volumes of the individual lattices are not independent
 * of one another.
 *
 * The molar volume of the Lattice solid is calculated from the following
 * formula
 *
 *  @f[
 *         V = \sum_i{ \theta_i V_i^{lattice}}
 *  @f]
 *
 * where @f$ V_i^{lattice} @f$ is the molar volume of the ith sublattice. This
 * is calculated from the following standard formula.
 *
 * @f[
 *     V_i = \sum_k{ X_k V_k}
 * @f]
 *
 * where k is a species in the ith sublattice.
 *
 * The mole fraction vector is redefined within the the LatticeSolidPhase object.
 * Each of the mole fractions sum to one on each of the sublattices.  The
 * routine getMoleFraction() and setMoleFraction() have been redefined to use
 * this convention.
 *
 * (This object is still under construction)
 */
class LatticeSolidPhase : public ThermoPhase
{
public:
    //! Base empty constructor
    LatticeSolidPhase() = default;

    string type() const override {
        return "compound-lattice";
    }

    //! String indicating the mechanical phase of the matter in this Phase.
    /*!
     * `LatticeSolid` phases only represent solids.
     */
    string phaseOfMatter() const override {
        return "solid";
    }

    bool isCompressible() const override {
        return false;
    }

    map<string, size_t> nativeState() const override {
        return { {"T", 0}, {"P", 1}, {"X", 2} };
    }

    double minTemp(size_t k=npos) const override;
    double maxTemp(size_t k=npos) const override;
    double refPressure() const override;

    //! This method returns the convention used in specification of the standard
    //! state, of which there are currently two, temperature based, and variable
    //! pressure based.
    /*!
     *  All of the thermo is determined by slave ThermoPhase routines.
     */
    int standardStateConvention() const override {
        return cSS_CONVENTION_SLAVE;
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where @f$
     * \theta_n @f$ is the lattice stoichiometric coefficient of the nth lattice
     *
     * @f[
     *   \tilde h(T,P) = {\sum_n \theta_n  \tilde h_n(T,P) }
     * @f]
     *
     * @f$ \tilde h_n(T,P) @f$ is the enthalpy of the nth lattice.
     *
     *  units J/kmol
     */
    double enthalpy_mole() const override;

    //! Return the Molar Internal Energy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where @f$
     * \theta_n @f$ is the lattice stoichiometric coefficient of the nth lattice
     *
     * @f[
     *   \tilde u(T,P) = {\sum_n \theta_n \tilde u_n(T,P) }
     * @f]
     *
     * @f$ \tilde u_n(T,P) @f$ is the internal energy of the nth lattice.
     *
     *  units J/kmol
     */
    double intEnergy_mole() const override;

    //! Return the Molar Entropy. Units: J/kmol/K.
    /*!
     * The molar enthalpy is determined by the following formula, where @f$
     * \theta_n @f$ is the lattice stoichiometric coefficient of the nth lattice
     *
     * @f[
     *   \tilde s(T,P) = \sum_n \theta_n \tilde s_n(T,P)
     * @f]
     *
     * @f$ \tilde s_n(T,P) @f$ is the molar entropy of the nth lattice.
     *
     *  units J/kmol/K
     */
    double entropy_mole() const override;

    //! Return the Molar Gibbs energy. Units: J/kmol.
    /*!
     * The molar Gibbs free energy is determined by the following formula, where
     * @f$ \theta_n @f$ is the lattice stoichiometric coefficient of the nth
     * lattice
     *
     * @f[
     *   \tilde h(T,P) = {\sum_n \theta_n \tilde h_n(T,P) }
     * @f]
     *
     * @f$ \tilde h_n(T,P) @f$ is the enthalpy of the nth lattice.
     *
     *  units J/kmol
     */
    double gibbs_mole() const override;

    //! Return the constant pressure heat capacity. Units: J/kmol/K
    /*!
     * The molar constant pressure heat capacity is determined by the following
     * formula, where @f$ C_n @f$ is the lattice molar density of the nth
     * lattice, and @f$ C_T @f$ is the molar density of the solid compound.
     *
     * @f[
     *   \tilde c_{p,n}(T,P) = \frac{\sum_n C_n \tilde c_{p,n}(T,P) }{C_T},
     * @f]
     *
     * @f$ \tilde c_{p,n}(T,P) @f$ is the heat capacity of the nth lattice.
     *
     *  units J/kmol/K
     */
    double cp_mole() const override;

    //! Return the constant volume heat capacity. Units: J/kmol/K
    /*!
     * The molar constant volume heat capacity is determined by the following
     * formula, where @f$ C_n @f$ is the lattice molar density of the nth
     * lattice, and @f$ C_T @f$ is the molar density of the solid compound.
     *
     * @f[
     *   \tilde c_{v,n}(T,P) = \frac{\sum_n C_n \tilde c_{v,n}(T,P) }{C_T},
     * @f]
     *
     * @f$ \tilde c_{v,n}(T,P) @f$ is the heat capacity of the nth lattice.
     *
     *  units J/kmol/K
     */
    double cv_mole() const override {
        return cp_mole();
    }

    //! Report the Pressure. Units: Pa.
    /*!
     *  This method simply returns the stored pressure value.
     */
    double pressure() const override {
        return m_press;
    }

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * @param p Pressure (units - Pa)
     */
    void setPressure(double p) override;

    //! Calculate the density of the solid mixture
    /*!
     * The formula for this is
     *
     * @f[
     *      \rho = \sum_n{ \rho_n \theta_n }
     * @f]
     *
     * where @f$ \rho_n @f$  is the density of the nth sublattice
     */
    double calcDensity();

    //! Set the mole fractions to the specified values, and then normalize them
    //! so that they sum to 1.0 for each of the subphases
    /*!
     * On input, the mole fraction vector is assumed to sum to one for each of
     * the sublattices. The sublattices are updated with this mole fraction
     * vector. The mole fractions are also stored within this object, after they
     * are normalized to one by dividing by the number of sublattices.
     *
     * @param x  Input vector of mole fractions. There is no restriction on the
     *           sum of the mole fraction vector. Internally, this object will
     *           pass portions of this vector to the sublattices which assume
     *           that the portions individually sum to one. Length is m_kk.
     */
    void setMoleFractions(const double* const x) override;

    //! Get the species mole fraction vector.
    /*!
     * On output the mole fraction vector will sum to one for each of the
     * subphases which make up this phase.
     *
     * @param x On return, x contains the mole fractions. Must have a length
     *          greater than or equal to the number of species.
     *
     * @todo This method shadows but does not override Phase::getMoleFractions(). The
     *       LatticeSolidPhase class should be revised to avoid needing to override the
     *       methods for getting composition. See
     *       https://github.com/Cantera/cantera/issues/1310 for additional information.
     */
    void getMoleFractions(double* const x) const;

    void setMassFractions(const double* const y) override {
        throw NotImplementedError("LatticeSolidPhase::setMassFractions");
    }

    void setMassFractions_NoNorm(const double* const y) override {
        throw NotImplementedError("LatticeSolidPhase::setMassFractions_NoNorm");
    }

    void getConcentrations(double* const c) const override {
        throw NotImplementedError("LatticeSolidPhase::getConcentrations");
    }

    double concentration(size_t k) const override {
        throw NotImplementedError("LatticeSolidPhase::concentration");
    }

    void setConcentrations(const double* const conc) override {
        throw NotImplementedError("LatticeSolidPhase::setConcentrations");
    }

    Units standardConcentrationUnits() const override;

    void getActivityConcentrations(double* c) const override;

    void getActivityCoefficients(double* ac) const override;

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solution at the current temperature, pressure and mole fraction of the
     * solution.
     *
     * This returns the underlying lattice chemical potentials, as the units are
     * kmol-1 of the sublattice species.
     *
     * @param mu  Output vector of species chemical potentials. Length: m_kk.
     *            Units: J/kmol
     */
    void getChemPotentials(double* mu) const override;

    //! Returns an array of partial molar enthalpies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol). For this phase, the partial molar enthalpies are equal to
     * the pure species enthalpies
     *  @f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * @f]
     * The reference-state pure-species enthalpies, @f$ \hat h^{ref}_k(T) @f$,
     * at the reference pressure,@f$ P_{ref} @f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    void getPartialMolarEnthalpies(double* hbar) const override;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K. For this phase, the partial molar entropies
     * are equal to the pure species entropies plus the ideal solution
     * contribution.
     *  @f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R \ln(X_k)
     * @f]
     * The reference-state pure-species entropies,@f$ \hat s^{ref}_k(T) @f$, at
     * the reference pressure, @f$ P_{ref} @f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    void getPartialMolarEntropies(double* sbar) const override;

    /**
     * Returns an array of partial molar Heat Capacities at constant pressure of
     * the species in the solution. Units: J/kmol/K. For this phase, the partial
     * molar heat capacities are equal to the standard state heat capacities.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    void getPartialMolarCp(double* cpbar) const override;

    /**
     * returns an array of partial molar volumes of the species in the solution.
     * Units: m^3 kmol-1.
     *
     * For this solution, the partial molar volumes are equal to the constant
     * species molar volumes.
     *
     * @param vbar  Output vector of partial molar volumes. Length: m_kk.
     */
    void getPartialMolarVolumes(double* vbar) const override;

    //! Get the array of standard state chemical potentials at unit activity for
    //! the species at their standard states at the current *T* and *P* of the
    //! solution.
    /*!
     * These are the standard state chemical potentials @f$ \mu^0_k(T,P) @f$.
     * The values are evaluated at the current temperature and pressure of the
     * solution.
     *
     * This returns the underlying lattice standard chemical potentials, as the
     * units are kmol-1 of the sublattice species.
     *
     * @param mu0    Output vector of chemical potentials.
     *               Length: m_kk. Units: J/kmol
     */
    void getStandardChemPotentials(double* mu0) const override;

    double standardConcentration(size_t k=0) const override;
    double logStandardConc(size_t k=0) const override;

    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    void getGibbs_RT_ref(double* grt) const override;
    void getGibbs_ref(double* g) const override;
    //! @}

    bool addSpecies(shared_ptr<Species> spec) override;

    //! Add a lattice to this phase
    void addLattice(shared_ptr<ThermoPhase> lattice);

    //! Set the lattice stoichiometric coefficients, @f$ \theta_i @f$
    void setLatticeStoichiometry(const Composition& comp);

    void setParameters(const AnyMap& phaseNode,
                       const AnyMap& rootNode=AnyMap()) override;
    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;

    //! Set the Lattice mole fractions using a string
    /*!
     * @param n  Integer value of the lattice whose mole fractions are being set
     * @param x  string containing Name:value pairs that will specify the mole
     *           fractions of species on a particular lattice
     */
    void setLatticeMoleFractionsByName(int n, const string& x);

    void modifyOneHf298SS(const size_t k, const double Hf298New) override;
    void resetHf298(const size_t k=npos) override;

protected:
    //! Current value of the pressure
    double m_press = -1.0;

    //! Current value of the molar density
    double m_molar_density = 0.0;

    //! Vector of sublattice ThermoPhase objects
    vector<shared_ptr<ThermoPhase>> m_lattice;

    //! Vector of mole fractions
    /*!
     *  Note these mole fractions sum to one when summed over all phases.
     *  However, this is not what's passed down to the lower m_lattice objects.
     */
    mutable vector<double> m_x;

    //! Lattice stoichiometric coefficients
    vector<double> theta_;

    //! Temporary vector
    mutable vector<double> tmpV_;

    vector<size_t> lkstart_;

    //! Root node of the AnyMap which contains this phase definition.
    //! Used to look up the phase definitions for the constituent phases.
    AnyMap m_rootNode;

private:
    //! Update the reference thermodynamic functions
    void _updateThermo() const;
};
}

#endif
