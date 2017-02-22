/**
 * @file GasQSSKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.
// Created by Hao Wu (wuhao@stanford.edu), Feb. 2017

#ifndef CT_GASQSSKINETICS_H
#define CT_GASQSSKINETICS_H

#include "GasKinetics.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera {

/**
 * Kinetics manager for elementary gas-phase chemistry with QSSA species.
 * @ingroup kinetics
 */
class GasQSSKinetics : public GasKinetics
{
public:

    GasQSSKinetics(doublereal *_buff_full = nullptr, thermo_t *thermo = 0);

    virtual ~GasQSSKinetics();

    virtual std::string kineticsType() const {
        return "GasQSSA";
    }

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(doublereal* kfwd);

    virtual void getDeltaGibbs(doublereal* deltaG);

    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    virtual void init();
    virtual bool addReaction(shared_ptr<Reaction> r);

    /**
     * Species creation rates [kmol/m^3/s or kmol/m^2/s]. Return the species
     * creation rates in array cdot, which must be dimensioned at least as
     * large as the total number of species in all phases. @see nTotalSpecies.
     *
     * @param cdot   Output vector of creation rates. Length: m_kk.
     */
    virtual void getCreationRates(doublereal* cdot);

    /**
     * Species destruction rates [kmol/m^3/s or kmol/m^2/s]. Return the species
     * destruction rates in array ddot, which must be dimensioned at least as
     * large as the total number of species. @see nTotalSpecies.
     *
     * @param ddot   Output vector of destruction rates. Length: m_kk.
     */
    virtual void getDestructionRates(doublereal* ddot);

    /**
     * Species net production rates [kmol/m^3/s or kmol/m^2/s]. Return the
     * species net production rates (creation - destruction) in array wdot,
     * which must be dimensioned at least as large as the total number of
     * species. @see nTotalSpecies.
     *
     * @param wdot   Output vector of net production rates. Length: m_kk.
     */
    virtual void getNetProductionRates(doublereal* wdot);

    void updateROP();

    virtual void update_rates_T();

    virtual void update_rates_C();

    //! Update the equilibrium constants in molar units.
    void updateKc();

    //! Calculate concentration for QSS species
    virtual void calc_conc_QSS(doublereal* conc_qss);

    //! Update rate for QSS
    virtual void update_ROP_QSS(const doublereal* conc_qss);

    //! Initialize rate info for QSS
    virtual void init_QSS();
protected:
    virtual bool addReactionQSS(shared_ptr<Reaction> r);
    //! can be removed after development
    virtual void printQSS();

    doublereal * const m_buff_full;
    double m_rel_density_qss;
    bool m_QSS_init;
    bool m_QSS_ok;
    size_t m_nSpeciesQSS;
    // index of destruction rxn for QSS speciies in forward direction
    // m_rodf_qss[i] = vector of reaction, whose reactants are QSS species
    std::vector<std::vector<size_t>> m_rodf_qss;
    // index of destruction rxn for QSS speciies in reverse direction
    // m_rodr_qss[i] = vector of reaction, whose products are QSS species
    std::vector<std::vector<size_t>> m_rodr_qss;
    // index of destruction rxn with multiple QSS speciies in forward direction
    std::vector<size_t> m_rodf_mqss;
    // index of destruction rxn with multiple QSS speciies in reverse direction
    std::vector<size_t> m_rodr_mqss;
    // index of production rxn for QSS speciies in forward direction
    // from non-qss species
    std::vector<std::vector<size_t>> m_ropf_noqss;
    // index of production rxn for QSS speciies in reverse direction
    // from non-qss species
    std::vector<std::vector<size_t>> m_ropr_noqss;
    // index of production rxn for QSS speciies in forward direction
    // from qss species [from k, to i]
    std::vector<std::vector<std::vector<size_t>>> m_ropf_qss_tmp;
    std::vector<std::vector<size_t>> m_ropf_qss;
    // index of production rxn for QSS speciies in reverse direction
    // from qss species [from k, to i]
    std::vector<std::vector<std::vector<size_t>>> m_ropr_qss_tmp;
    std::vector<std::vector<size_t>> m_ropr_qss;
    //
    static const char IROPF = 1;
    static const char IROPR = 1 << 1;
    std::vector<char> m_ifr_qss;
    //
    Eigen::VectorXd m_rod_qss;
    Eigen::SparseMatrix<double> m_rop_qss;
    Eigen::VectorXd m_rop_noqss;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        m_solver_qss;

};
}

#endif
