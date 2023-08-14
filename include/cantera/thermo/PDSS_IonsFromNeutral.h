/**
 *  @file PDSS_IonsFromNeutral.h
 *   Declarations for the class PDSS_IonsFromNeutral (
 *    which handles calculations for a single ion in a fluid, whose properties
 *    are calculated from another neutral molecule.
 *    (see @ref pdssthermo and class @link Cantera::PDSS_IonsFromNeutral PDSS_IonsFromNeutral@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_IONSFROMNEUTRAL_H
#define CT_PDSS_IONSFROMNEUTRAL_H

#include "PDSS.h"

namespace Cantera
{
class ThermoPhase;

//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to %Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of  %Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * This class is for a single Ideal Gas species.
 *
 * @deprecated To be removed after %Cantera 3.0
 *
 * @ingroup pdssthermo
 */
class PDSS_IonsFromNeutral : public PDSS_Nondimensional
{
public:
    //! Default constructor
    PDSS_IonsFromNeutral();

    //! @name  Molar Thermodynamic Properties of the Species Standard State
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    double enthalpy_RT() const override;
    double intEnergy_mole() const override;
    double entropy_R() const override;

    /**
     * @copydoc PDSS::gibbs_RT()
     *
     * @f[
     *   \frac{\mu^o_k}{RT} = \sum_{m}{ \alpha_{m , k} \frac{\mu^o_{m}}{RT}} + ( 1 - \delta_{k,sp}) 2.0 \ln{2.0}
     * @f]
     *
     * *m* is the neutral molecule species index. @f$ \alpha_{m , k} @f$ is the
     * stoichiometric coefficient for the neutral molecule, *m*, that creates the
     * thermodynamics for the ionic species  *k*. A factor  @f$ 2.0 \ln{2.0} @f$
     * is added to all ions except for the species ionic species, which in this
     * case is the single anion species, with species index *sp*.
     */
    double gibbs_RT() const override;
    double cp_R() const override;
    double molarVolume() const override;
    double density() const override;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    double gibbs_RT_ref() const override;
    double enthalpy_RT_ref() const override;
    double entropy_R_ref() const override;
    double cp_R_ref() const override;
    double molarVolume_ref() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    void setState_TP(double temp, double pres) override;

    //! @}
    //! @name Initialization of the Object
    //! @{

    void setParent(VPStandardStateTP* phase, size_t k) override;

    void setNeutralSpeciesMultiplier(const string& species, double mult);
    void setSpecialSpecies(bool special=true);
    void getParameters(AnyMap& eosNode) const override;
    void initThermo() override;
    //! @}

protected:
    //! Pointer to the Neutral Molecule ThermoPhase object
    shared_ptr<ThermoPhase> neutralMoleculePhase_;

    map<string, double> neutralSpeciesMultipliers_;

    //! Number of neutral molecule species that make up the stoichiometric
    //! vector for this species, in terms of calculating thermodynamic functions
    size_t numMult_ = 0;

    //! Vector of species indices in the neutral molecule ThermoPhase
    vector<size_t> idNeutralMoleculeVec;

    //! Stoichiometric coefficient for this species using the Neutral Molecule
    //! Species in the vector idNeutralMoleculeVec
    vector<double> factorVec;

    //! Add 2RTln2 to the entropy and Gibbs free energies for this species
    /*!
     *  This is true if this species is not the special species
     */
    bool add2RTln2_ = true;

    //! Vector of length equal to the number of species in the neutral molecule phase
    mutable vector<double> tmpNM;
};
}

#endif
