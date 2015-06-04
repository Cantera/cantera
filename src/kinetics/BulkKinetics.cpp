#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

BulkKinetics::BulkKinetics(thermo_t* thermo) :
    m_ROP_ok(false),
    m_temp(0.0),
    m_finalized(false)
{
    if (thermo) {
        addPhase(*thermo);
    }
}


Kinetics* BulkKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    BulkKinetics* kin = new BulkKinetics(*this);
    kin->assignShallowPointers(tpVector);
    return kin;
}

bool BulkKinetics::isReversible(size_t i) {
    if (std::find(m_revindex.begin(), m_revindex.end(), i)
            < m_revindex.end()) {
        return true;
    } else {
        return false;
    }
}


void BulkKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the ideal gas solution.
    thermo().getChemPotentials(&m_grt[0]);
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(&m_grt[0], deltaG);
}

void BulkKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species in the ideal gas.
    thermo().getPartialMolarEnthalpies(&m_grt[0]);
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(&m_grt[0], deltaH);
}

void BulkKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in the solid solution.
    thermo().getPartialMolarEntropies(&m_grt[0]);
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(&m_grt[0], deltaS);
}

void BulkKinetics::getDeltaSSGibbs(doublereal* deltaG)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity. We define these here as
    // the chemical potentials of the pure species at the temperature and
    // pressure of the solution.
    thermo().getStandardChemPotentials(&m_grt[0]);
    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(&m_grt[0], deltaG);
}

void BulkKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    // Get the standard state enthalpies of the species.
    thermo().getEnthalpy_RT(&m_grt[0]);
    doublereal RT = thermo().temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= RT;
    }
    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(&m_grt[0], deltaH);
}

void BulkKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    // Get the standard state entropy of the species. We define these here as
    // the entropies of the pure species at the temperature and pressure of the
    // solution.
    thermo().getEntropy_R(&m_grt[0]);
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= R;
    }
    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(&m_grt[0], deltaS);
}

void BulkKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
{
    /*
     * go get the forward rate constants. -> note, we don't
     * really care about speed or redundancy in these
     * informational routines.
     */
    getFwdRateConstants(krev);

    if (doIrreversible) {
        getEquilibriumConstants(&m_ropnet[0]);
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] /=  m_ropnet[i];
        }
    } else {
        // m_rkcn[] is zero for irreversible reactions
        for (size_t i = 0; i < m_ii; i++) {
            krev[i] *= m_rkcn[i];
        }
    }
}

void BulkKinetics::addReaction(ReactionData& r)
{
    Kinetics::addReaction(r);
    m_dn.push_back(accumulate(r.pstoich.begin(), r.pstoich.end(), 0.0) -
                   accumulate(r.rstoich.begin(), r.rstoich.end(), 0.0));

    if (r.reversible) {
        m_revindex.push_back(nReactions());
    } else {
        m_irrev.push_back(nReactions());
    }
}

bool BulkKinetics::addReaction(shared_ptr<Reaction> r)
{
    bool added = Kinetics::addReaction(r);
    if (!added) {
        return false;
    }
    double dn = 0.0;
    for (Composition::const_iterator iter = r->products.begin();
         iter != r->products.end();
         ++iter) {
        dn += iter->second;
    }
    for (Composition::const_iterator iter = r->reactants.begin();
         iter != r->reactants.end();
         ++iter) {
        dn -= iter->second;
    }

    m_dn.push_back(dn);

    if (r->reversible) {
        m_revindex.push_back(nReactions()-1);
    } else {
        m_irrev.push_back(nReactions()-1);
    }
    return true;
}

void BulkKinetics::addElementaryReaction(ReactionData& r)
{
    m_rates.install(nReactions(), r);
}

void BulkKinetics::addElementaryReaction(ElementaryReaction& r)
{
    m_rates.install(nReactions()-1, r.rate);
}

void BulkKinetics::modifyElementaryReaction(size_t i, ElementaryReaction& rNew)
{
    m_rates.replace(i, rNew.rate);
}

void BulkKinetics::init()
{
    m_kk = thermo().nSpecies();
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_conc.resize(m_kk);
    m_grt.resize(m_kk);
}

void BulkKinetics::finalize()
{
    m_finalized = true;

    // Guarantee that these arrays can be converted to double* even in the
    // special case where there are no reactions defined.
    if (!m_ii) {
        m_perturb.resize(1, 1.0);
        m_ropf.resize(1, 0.0);
        m_ropr.resize(1, 0.0);
        m_ropnet.resize(1, 0.0);
        m_rkcn.resize(1, 0.0);
    }
}

bool BulkKinetics::ready() const
{
    return m_finalized;
}

void BulkKinetics::setMultiplier(size_t i, double f) {
    Kinetics::setMultiplier(i, f);
    m_ROP_ok = false;
}

}
