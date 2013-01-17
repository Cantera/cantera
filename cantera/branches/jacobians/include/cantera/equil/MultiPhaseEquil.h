#ifndef CT_MULTIPHASE_EQUIL
#define CT_MULTIPHASE_EQUIL

#include "cantera/base/ct_defs.h"
#include "MultiPhase.h"

namespace Cantera
{

/**
 * Multiphase chemical equilibrium solver.  Class MultiPhaseEquil
 * is designed to be used to set a mixture containing one or more
 * phases to a state of chemical equilibrium.  It implements the
 * VCS algorithm, described in Smith and Missen, "Chemical
 * Reaction Equilibrium."
 *
 * This class only handles chemical equilibrium at a specified
 * temperature and pressure. To compute equilibrium holding other
 * properties fixed, it is necessary to iterate on T and P in an
 * "outer" loop, until the specified properties have the desired
 * values. This is done, for example, in method equilibrate of
 * class MultiPhase.
 *
 * This class is primarily meant to be used internally by the
 * equilibrate method of class MultiPhase, although there is no
 * reason it cannot be used directly in application programs if
 * desired.
 *
 * @ingroup equil
 */

class MultiPhaseEquil
{

public:
    MultiPhaseEquil(MultiPhase* mix, bool start=true, int loglevel = 0);

    virtual ~MultiPhaseEquil() {}

    size_t constituent(size_t m) {
        if (m < m_nel) {
            return m_order[m];
        } else {
            return npos;
        }
    }

    void getStoichVector(size_t rxn, vector_fp& nu) {
        size_t k;
        nu.resize(m_nsp, 0.0);
        if (rxn > nFree()) {
            return;
        }
        for (k = 0; k < m_nsp; k++) {
            nu[m_order[k]] = m_N(k, rxn);
        }
    }

    int iterations() {
        return m_iter;
    }

    doublereal equilibrate(int XY, doublereal err = 1.0e-9,
                           int maxsteps = 1000, int loglevel=-99);
    doublereal error();

#if defined(WITH_HTML_LOGS)
    std::string reactionString(size_t j);
    void printInfo(int loglevel);
#else
    inline std::string reactionString(size_t j) {
        return std::string("");
    }
    inline void printInfo(int loglevel) {}
#endif

    void setInitialMixMoles(int loglevel = 0) {
        setInitialMoles(loglevel);
        finish();
    }

    size_t componentIndex(size_t n) {
        return m_species[m_order[n]];
    }

    void reportCSV(const std::string& reportFile);

    double phaseMoles(size_t iph) const;

protected:

    void getComponents(const std::vector<size_t>& order);
    int setInitialMoles(int loglevel = 0);
    void computeN();
    doublereal stepComposition(int loglevel = 0);
    //void sort(vector_fp& x);
    void unsort(vector_fp& x);
    void step(doublereal omega, vector_fp& deltaN, int loglevel = 0);
    doublereal computeReactionSteps(vector_fp& dxi);
    void updateMixMoles();
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
    doublereal mu(size_t n) const {
        return m_mu[m_species[m_order[n]]];
    }
    std::string speciesName(size_t n) const {
        return
            m_mix->speciesName(m_species[m_order[n]]);
    }

    //! Number of degrees of freedom
    size_t nFree() const {
        return (m_nsp > m_nel) ? m_nsp - m_nel : 0;
    }

    size_t m_nel_mix, m_nsp_mix, m_np;
    size_t m_nel, m_nsp;
    size_t m_eloc;
    int m_iter;
    MultiPhase* m_mix;
    doublereal m_press, m_temp;
    std::vector<size_t> m_order;
    DenseMatrix m_N, m_A;
    vector_fp m_work, m_work2, m_work3;
    vector_fp m_moles, m_lastmoles, m_dxi;
    vector_fp m_deltaG_RT, m_mu;
    std::vector<bool> m_majorsp;
    std::vector<size_t> m_sortindex;
    vector_int m_lastsort;
    vector_int m_dsoln;
    vector_int m_incl_element, m_incl_species;

    // Vector of indices for species that are included in the
    // calculation.  This is used to exclude pure-phase species
    // with invalid thermo data
    std::vector<size_t> m_species;
    std::vector<size_t> m_element;
    std::vector<bool> m_solnrxn;
    bool m_force;
};

}


#endif
