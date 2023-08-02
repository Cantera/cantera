//! @file MultiPhaseEquil.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIPHASE_EQUIL
#define CT_MULTIPHASE_EQUIL

#include "MultiPhase.h"

namespace Cantera
{

/**
 * Multiphase chemical equilibrium solver. Class MultiPhaseEquil is designed
 * to be used to set a mixture containing one or more phases to a state of
 * chemical equilibrium. It implements the VCS algorithm, described in Smith
 * and Missen @cite smith1982.
 *
 * This class only handles chemical equilibrium at a specified temperature and
 * pressure. To compute equilibrium holding other properties fixed, it is
 * necessary to iterate on T and P in an "outer" loop, until the specified
 * properties have the desired values. This is done, for example, in method
 * equilibrate of class MultiPhase.
 *
 * This class is primarily meant to be used internally by the equilibrate
 * method of class MultiPhase, although there is no reason it cannot be used
 * directly in application programs if desired.
 *
 * @ingroup equilGroup
 */
class MultiPhaseEquil
{
public:
    //! Construct a multiphase equilibrium manager for a multiphase mixture.
    //! @param mix Pointer to a multiphase mixture object.
    //! @param start If true, the initial composition will be determined by a
    //!     linear Gibbs minimization, otherwise the initial mixture
    //!     composition will be used.
    //! @param loglevel Desired level of debug printing. loglevel = 0 suppresses
    //!     printing. Higher values request more verbose logging output.
    MultiPhaseEquil(MultiPhase* mix, bool start=true, int loglevel = 0);

    virtual ~MultiPhaseEquil() {}

    size_t constituent(size_t m) {
        if (m < m_nel) {
            return m_order[m];
        } else {
            return npos;
        }
    }

    void getStoichVector(size_t rxn, vector<double>& nu) {
        nu.resize(m_nsp, 0.0);
        if (rxn > nFree()) {
            return;
        }
        for (size_t k = 0; k < m_nsp; k++) {
            nu[m_order[k]] = m_N(k, rxn);
        }
    }

    int iterations() {
        return m_iter;
    }

    double equilibrate(int XY, double err = 1.0e-9,
                       int maxsteps = 1000, int loglevel=-99);
    double error();

    string reactionString(size_t j) {
        return "";
    }
    void setInitialMixMoles(int loglevel = 0) {
        setInitialMoles(loglevel);
        finish();
    }

    size_t componentIndex(size_t n) {
        return m_species[m_order[n]];
    }

    void reportCSV(const string& reportFile);

    double phaseMoles(size_t iph) const;

protected:
    //! This method finds a set of component species and a complete set of
    //! formation reactions for the non-components in terms of the components.
    //! In most cases, many different component sets are possible, and
    //! therefore neither the components returned by this method nor the
    //! formation reactions are unique. The algorithm used here is described
    //! in Smith and Missen @cite smith1982.
    //!
    //! The component species are taken to be the first M species in array
    //! 'species' that have linearly-independent compositions.
    //!
    //! @param order On entry, vector @e order should contain species index
    //!     numbers in the order of decreasing desirability as a component.
    //!     For example, if it is desired to choose the components from among
    //!     the major species, this array might list species index numbers in
    //!     decreasing order of mole fraction. If array 'species' does not
    //!     have length = nSpecies(), then the species will be considered as
    //!     candidates to be components in declaration order, beginning with
    //!     the first phase added.
    void getComponents(const vector<size_t>& order);

    //! Estimate the initial mole numbers. This is done by running each
    //! reaction as far forward or backward as possible, subject to the
    //! constraint that all mole numbers remain non-negative. Reactions for
    //! which @f$ \Delta \mu^0 @f$ are positive are run in reverse, and ones
    //! for which it is negative are run in the forward direction. The end
    //! result is equivalent to solving the linear programming problem of
    //! minimizing the linear Gibbs function subject to the element and non-
    //! negativity constraints.
    int setInitialMoles(int loglevel = 0);

    void computeN();

    //! Take one step in composition, given the gradient of G at the starting
    //! point, and a vector of reaction steps dxi.
    double stepComposition(int loglevel = 0);

    //! Re-arrange a vector of species properties in sorted form
    //! (components first) into unsorted, sequential form.
    void unsort(vector<double>& x);

    void step(double omega, vector<double>& deltaN, int loglevel = 0);

    //! Compute the change in extent of reaction for each reaction.
    double computeReactionSteps(vector<double>& dxi);

    void updateMixMoles();

    //! Clean up the composition. The solution algorithm can leave some
    //! species in stoichiometric condensed phases with very small negative
    //! mole numbers. This method simply sets these to zero.
    void finish();

    // moles of the species with sorted index ns
    double moles(size_t ns) const {
        return m_moles[m_order[ns]];
    }
    double& moles(size_t ns) {
        return m_moles[m_order[ns]];
    }
    int solutionSpecies(size_t n) const {
        return m_dsoln[m_order[n]];
    }
    bool isStoichPhase(size_t n) const {
        return (m_dsoln[m_order[n]] == 0);
    }
    double mu(size_t n) const {
        return m_mu[m_species[m_order[n]]];
    }
    string speciesName(size_t n) const {
        return
            m_mix->speciesName(m_species[m_order[n]]);
    }

    //! Number of degrees of freedom
    size_t nFree() const {
        return (m_nsp > m_nel) ? m_nsp - m_nel : 0;
    }

    size_t m_nel_mix, m_nsp_mix;
    size_t m_nel = 0;
    size_t m_nsp = 0;
    size_t m_eloc = 1000;
    int m_iter;
    MultiPhase* m_mix;
    double m_press, m_temp;
    vector<size_t> m_order;
    DenseMatrix m_N, m_A;
    vector<double> m_work, m_work2, m_work3;
    vector<double> m_moles, m_lastmoles, m_dxi;
    vector<double> m_deltaG_RT, m_mu;
    vector<bool> m_majorsp;
    vector<size_t> m_sortindex;
    vector<int> m_lastsort;
    vector<int> m_dsoln;
    vector<int> m_incl_element, m_incl_species;

    // Vector of indices for species that are included in the calculation.
    // This is used to exclude pure-phase species with invalid thermo data
    vector<size_t> m_species;
    vector<size_t> m_element;
    vector<bool> m_solnrxn;
    bool m_force = true;
};

}

#endif
