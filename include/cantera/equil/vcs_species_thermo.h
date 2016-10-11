//! @file vcs_species_thermo.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef VCS_SPECIES_THERMO_H
#define VCS_SPECIES_THERMO_H

#include <cstdlib>

namespace Cantera
{

class vcs_VolPhase;

// Models for the species standard state Naught temperature dependence
#define VCS_SS0_NOTHANDLED -1
#define VCS_SS0_CONSTANT 0
//#define VCS_SS0_NASA_POLY 1
#define VCS_SS0_CONSTANT_CP 2

// Models for the species standard state extra pressure dependence
#define VCS_SSSTAR_NOTHANDLED -1
#define VCS_SSSTAR_CONSTANT 0
#define VCS_SSSTAR_IDEAL_GAS 1

/*!
 * Identifies the thermo model for the species. This structure is shared by
 * volumetric and surface species. However, each will have its own types of
 * thermodynamic models. These quantities all have appropriate units.
 */
class VCS_SPECIES_THERMO
{
    /*
     * All objects are public for ease of development
     */
public:
    //! Index of the phase that this species belongs to.
    size_t IndexPhase;

    //! Index of this species in the current phase.
    size_t IndexSpeciesPhase;

    //! Pointer to the owning phase object.
    vcs_VolPhase* OwningPhase;

    //! Integer representing the models for the species standard state Naught
    //! temperature dependence. They are listed above and start with VCS_SS0_...
    int SS0_Model;

    //! Internal storage of the last calculation of the reference naught Gibbs
    //! free energy at SS0_TSave. (always in units of Kelvin)
    double SS0_feSave;

    //! Internal storage of the last temperature used in the calculation of the
    //! reference naught Gibbs free energy. units = kelvin
    double SS0_TSave;

    //! Base temperature used in the VCS_SS0_CONSTANT_CP model
    double SS0_T0;

    //! Base enthalpy used in the VCS_SS0_CONSTANT_CP model
    double SS0_H0;

    //! Base entropy used in the VCS_SS0_CONSTANT_CP model
    double SS0_S0;

    //! Base heat capacity used in the VCS_SS0_CONSTANT_CP model
    double SS0_Cp0;

    //! Value of the pressure for the reference state.
    //! defaults to 1.01325E5 = 1 atm
    double SS0_Pref;

    //! Integer value representing the star state model.
    int SSStar_Model;

    //! Models for the standard state volume of each species
    int SSStar_Vol_Model;

    //! parameter that is used in the VCS_SSVOL_CONSTANT model.
    double SSStar_Vol0;

    VCS_SPECIES_THERMO(size_t indexPhase, size_t indexSpeciesPhase);
    virtual ~VCS_SPECIES_THERMO() {}

    //! Duplication function for derived classes.
    virtual VCS_SPECIES_THERMO* duplMyselfAsVCS_SPECIES_THERMO();

    /**
     * This function calculates the standard state Gibbs free energy
     * for species, kspec, at the temperature TKelvin and pressure, Pres.
     *
     * @param kspec     species global index
     * @param TKelvin   Temperature in Kelvin
     * @param pres      pressure in Pa
     * @return standard state free energy in units of Kelvin.
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    virtual double GStar_R_calc(size_t kspec, double TKelvin, double pres);

    /**
     * This function calculates the standard state Gibbs free energy for
     * species, kspec, at the temperature TKelvin
     *
     * @param kglob    species global index.
     * @param TKelvin  Temperature in Kelvin
     * @return standard state free energy in Kelvin.
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    virtual double G0_R_calc(size_t kglob, double TKelvin);

    /**
     * This function calculates the standard state molar volume for species,
     * kspec, at the temperature TKelvin and pressure, Pres,
     *
     * @return standard state volume in m**3 / kmol
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    virtual double VolStar_calc(size_t kglob, double TKelvin, double Pres);

    /**
     * This function evaluates the activity coefficient for species, kspec
     *
     * Note, T, P and mole fractions are obtained from the single private
     * instance of VCS_SOLVE
     *
     * @param kspec index of the species in the global species list within
     *     VCS_SOLVE. Phase and local species id can be looked up within object.
     * @return activity coefficient for species kspec
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    virtual double eval_ac(size_t kspec);
};

}

#endif
