/**
 * @file vcs_prob.h
 *  Header for the Interface class for the vcs thermo equilibrium solver package,
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef _VCS_PROB_H
#define _VCS_PROB_H

#include "cantera/base/Array.h"

namespace Cantera
{

class vcs_VolPhase;
class VCS_SPECIES_THERMO;

//! Interface class for the vcs thermo equilibrium solver package,
//! which generally describes the problem to be solved.
class VCS_PROB
{
public:
    //! Problem type. I.e., the identity of what is held constant. Currently, T
    //! and P are held constant, and this input is ignored
    int prob_type;

    //! Total number of species in the problems
    size_t nspecies;

    //! Species number used to size data structures
    size_t NSPECIES0;

    //! Number of element constraints in the equilibrium problem
    size_t ne;

    //! Number of element constraints used to size data structures
    //! involving elements
    size_t NE0;

    //! Number of phases in the problem
    size_t NPhase;

    //! Number of phases used to size data structures
    size_t NPHASE0;

    //! Vector of chemical potentials of the species. This is a calculated
    //! output quantity. length = number of species.
    vector_fp m_gibbsSpecies;

    //! Total number of moles of the kth species.
    /*!
     * This is both an input and an output variable. On input, this is an
     * estimate of the mole numbers. The actual element abundance vector
     * contains the problem specification.
     *
     * On output, this contains the solution for the total number of moles of
     * the kth species.
     */
    vector_fp w;

    //! Mole fraction vector. This is a calculated vector, calculated from w[].
    //! length number of species.
    vector_fp mf;

    //! Element abundances for jth element
    /*!
     * This is input from the input file and is considered a constant from
     * thereon within the vcs_solve_TP().
     */
    vector_fp gai;

    //! Formula Matrix for the problem
    /*!
     * FormulaMatrix(kspec,j) = Number of elements, j, in the kspec species
     */
    Array2D FormulaMatrix;

    //! Specifies the species unknown type
    /*!
     * There are two types. One is the straightforward species, with the mole
     * number w[k], as the unknown. The second is the an interfacial voltage
     * where w[k] refers to the interfacial voltage in volts.
     *
     * These species types correspond to metallic electrons corresponding to
     * electrodes. The voltage and other interfacial conditions sets up an
     * interfacial current, which is set to zero in this initial treatment.
     * Later we may have non-zero interfacial currents.
     */
    vector_int SpeciesUnknownType;

    //! Temperature (Kelvin)
    /*!
     * Specification of the temperature for the equilibrium problem
     */
    double T;

    //! Pressure
    double PresPA;

    //! Volume of the entire system
    /*!
     * Note, this is an output variable atm
     */
    double Vol;

    //! Partial Molar Volumes of species
    /*!
     * This is a calculated vector, calculated from w[].
     * length number of species.
     */
    vector_fp VolPM;

    //! Specification of the initial estimate method
    /*!
     *  * 0: user estimate
     *  * 1: user estimate if satisifies elements
     *  * -1: machine estimate
     */
    int iest;

    //! Tolerance requirement for major species
    double tolmaj;

    //! Tolerance requirement for minor species
    double tolmin;

    //! Mapping between the species and the phases
    std::vector<size_t> PhaseID;

    //! Vector of strings containing the species names
    std::vector<std::string> SpName;

    //! vector of strings containing the element names
    std::vector<std::string> ElName;

    //! vector of Element types
    vector_int m_elType;

    //! Specifies whether an element constraint is active
    /*!
     * The default is true
     * Length = nelements
     */
    vector_int ElActive;

    //! Molecular weight of species
    /*!
     * WtSpecies[k] = molecular weight of species in gm/mol
     */
    vector_fp WtSpecies;

    //! Charge of each species
    vector_fp Charge;

    //! Array of phase structures
    std::vector<vcs_VolPhase*> VPhaseList;

    // String containing the title of the run
    std::string Title;

    //! Vector of pointers to thermo structures which identify the model and
    //! parameters for evaluating the thermodynamic functions for that
    //! particular species
    std::vector<VCS_SPECIES_THERMO*> SpeciesThermo;

    //! Number of iterations. This is an output variable
    int m_Iterations;

    //! Number of basis optimizations used. This is an output variable.
    int m_NumBasisOptimizations;

    //! Print level for print routines
    int m_printLvl;

    //! Debug print lvl
    int vcs_debug_print_lvl;

    //! Constructor
    /*!
     * This constructor initializes the sizes within the object to parameter
     * values.
     *
     * @param nsp number of species
     * @param nel number of elements
     * @param nph number of phases
     */
    VCS_PROB(size_t nsp, size_t nel, size_t nph);

    ~VCS_PROB();

    //! Resizes all of the element lists within the structure
    /*!
     * Note, this doesn't change the number of element constraints in the
     * problem. It will change #NE0 if `nel` is greater than #NE0.
     *
     * @param nel     size to dimension all the elements lists
     * @param force   If true, this will dimension the size to be equal to `nel`
     *                even if `nel` is less than the current value of #NE0
     */
    void resizeElements(size_t nel, int force);

    //! Calculate the element abundance vector from the mole numbers
    void set_gai();

    //! Print out the  problem specification in all generality as it currently
    //! exists in the VCS_PROB object
    /*!
     *  @param print_lvl Parameter lvl for printing
     *      * 0 - no printing
     *      * 1 - all printing
     */
    void prob_report(int print_lvl);

    //! Add elements to the local element list
    /*!
     * This routine sorts through the elements defined in the vcs_VolPhase
     * object. It then adds the new elements to the VCS_PROB object, and creates
     * a global map, which is stored in the vcs_VolPhase object. Id and matching
     * of elements is done strictly via the element name, with case not
     * mattering.
     *
     * The routine also fills in the position of the element in the vcs_VolPhase
     * object's ElGlobalIndex field.
     *
     * @param volPhase  Object containing the phase to be added. The elements in
     *     this phase are parsed for addition to the global element list
     */
    void addPhaseElements(vcs_VolPhase* volPhase);

    //! This routine resizes the number of elements in the VCS_PROB object by
    //! adding a new element to the end of the element list
    /*!
     * The element name is added. Formula vector entries ang element abundances
     * for the new element are set to zero.
     *
     * @param elNameNew New name of the element
     * @param elType    Type of the element
     * @param elactive  boolean indicating whether the element is active
     * @returns the index number of the new element
     */
    size_t addElement(const char* elNameNew, int elType, int elactive);

    //! This routines adds entries for the formula matrix for one species
    /*!
     * This routines adds entries for the formula matrix for this object for one
     * species
     *
     * This object also fills in the index filed, IndSpecies, within the
     * volPhase object.
     *
     * @param volPhase object containing the species
     * @param k        Species number within the volPhase k
     * @param kT       global Species number within this object
     *
     */
    size_t addOnePhaseSpecies(vcs_VolPhase* volPhase, size_t k, size_t kT);

    void reportCSV(const std::string& reportFile);

    //! Set the debug level
    /*!
     * @param vcs_debug_print_lvl input debug level
     */
    void setDebugPrintLvl(int vcs_debug_print_lvl);
};

}

#endif
