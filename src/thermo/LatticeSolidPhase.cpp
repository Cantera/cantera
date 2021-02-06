/**
 *  @file LatticeSolidPhase.cpp
 *  Definitions for a simple thermodynamics model of a bulk solid phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
LatticeSolidPhase::LatticeSolidPhase() :
    m_press(-1.0),
    m_molar_density(0.0)
{
}

doublereal LatticeSolidPhase::minTemp(size_t k) const
{
    if (k != npos) {
        for (size_t n = 0; n < m_lattice.size(); n++) {
            if (lkstart_[n+1] < k) {
                return m_lattice[n]->minTemp(k-lkstart_[n]);
            }
        }
    }
    doublereal mm = 1.0E300;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        double ml = m_lattice[n]->minTemp();
        mm = std::min(mm, ml);
    }
    return mm;
}

doublereal LatticeSolidPhase::maxTemp(size_t k) const
{
    if (k != npos) {
        for (size_t n = 0; n < m_lattice.size(); n++) {
            if (lkstart_[n+1] < k) {
                return (m_lattice[n])->maxTemp(k - lkstart_[n]);
            }
        }
    }
    doublereal mm = -1.0E300;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        double ml = m_lattice[n]->maxTemp();
        mm = std::max(mm, ml);
    }
    return mm;
}

doublereal LatticeSolidPhase::refPressure() const
{
    return m_lattice[0]->refPressure();
}

doublereal LatticeSolidPhase::enthalpy_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->enthalpy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::intEnergy_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->intEnergy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::entropy_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->entropy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::gibbs_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->gibbs_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::cp_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->cp_mole();
    }
    return sum;
}

Units LatticeSolidPhase::standardConcentrationUnits() const
{
    return Units(1.0);
}

void LatticeSolidPhase::getActivityConcentrations(doublereal* c) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        m_lattice[n]->getMoleFractions(c+strt);
        strt += m_lattice[n]->nSpecies();
    }
}

void LatticeSolidPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

doublereal LatticeSolidPhase::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal LatticeSolidPhase::logStandardConc(size_t k) const
{
    return 0.0;
}

void LatticeSolidPhase::setPressure(doublereal p)
{
    m_press = p;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        m_lattice[n]->setPressure(m_press);
    }
    calcDensity();
}

doublereal LatticeSolidPhase::calcDensity()
{
    double sum = 0.0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        sum += theta_[n] * m_lattice[n]->density();
    }
    Phase::assignDensity(sum);
    return sum;
}

void LatticeSolidPhase::setMoleFractions(const doublereal* const x)
{
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nsp = m_lattice[n]->nSpecies();
        m_lattice[n]->setMoleFractions(x + strt);
        strt += nsp;
    }
    for (size_t k = 0; k < strt; k++) {
        m_x[k] = x[k] / m_lattice.size();
    }
    Phase::setMoleFractions(m_x.data());
    calcDensity();
}

void LatticeSolidPhase::getMoleFractions(doublereal* const x) const
{
    size_t strt = 0;
    // the ifdef block should be the way we calculate this.!!!!!
    Phase::getMoleFractions(x);
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nsp = m_lattice[n]->nSpecies();
        double sum = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            sum += (x + strt)[k];
        }
        for (size_t k = 0; k < nsp; k++) {
            (x + strt)[k] /= sum;
        }

        // At this point we can check against the mole fraction vector of the
        // underlying LatticePhase objects and get the same answer.
        m_lattice[n]->getMoleFractions(&m_x[strt]);
        for (size_t k = 0; k < nsp; k++) {
            if (fabs((x + strt)[k] - m_x[strt+k]) > 1.0E-14) {
                throw CanteraError("LatticeSolidPhase::getMoleFractions",
                                   "internal error");
            }
        }
        strt += nsp;
    }
}

void LatticeSolidPhase::getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nlsp = m_lattice[n]->nSpecies();
        m_lattice[n]->getChemPotentials(mu+strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nlsp = m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarEnthalpies(hbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nlsp = m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarEntropies(sbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarCp(doublereal* cpbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nlsp = m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarCp(cpbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nlsp = m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarVolumes(vbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getStandardChemPotentials(doublereal* mu0) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        m_lattice[n]->getStandardChemPotentials(mu0+strt);
        strt += m_lattice[n]->nSpecies();
    }
}

void LatticeSolidPhase::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    for (size_t n = 0; n < m_lattice.size(); n++) {
        m_lattice[n]->getGibbs_RT_ref(grt + lkstart_[n]);
    }
}

void LatticeSolidPhase::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void LatticeSolidPhase::setParameters(const AnyMap& phaseNode,
                                      const AnyMap& rootNode)
{
    ThermoPhase::setParameters(phaseNode, rootNode);
    m_rootNode = rootNode;
}

void LatticeSolidPhase::initThermo()
{
    if (m_input.hasKey("composition")) {
        compositionMap composition = m_input["composition"].asMap<double>();
        for (auto& item : composition) {
            AnyMap& node = m_rootNode["phases"].getMapWhere("name", item.first);
            addLattice(newPhase(node, m_rootNode));
        }
        setLatticeStoichiometry(composition);
    }

    setMoleFractions(m_x.data());
    ThermoPhase::initThermo();
}

bool LatticeSolidPhase::addSpecies(shared_ptr<Species> spec)
{
    // Species are added from component phases in addLattice()
    return false;
}

void LatticeSolidPhase::addLattice(shared_ptr<ThermoPhase> lattice)
{
    m_lattice.push_back(lattice);
    if (lkstart_.empty()) {
        lkstart_.push_back(0);
    }
    lkstart_.push_back(lkstart_.back() + lattice->nSpecies());

    if (theta_.size() == 0) {
        theta_.push_back(1.0);
    } else {
        theta_.push_back(0.0);
    }

    for (size_t k = 0; k < lattice->nSpecies(); k++) {
        ThermoPhase::addSpecies(lattice->species(k));
        vector_fp constArr(lattice->nElements());
        const vector_fp& aws = lattice->atomicWeights();
        for (size_t es = 0; es < lattice->nElements(); es++) {
            addElement(lattice->elementName(es), aws[es], lattice->atomicNumber(es),
                       lattice->entropyElement298(es), lattice->elementType(es));
        }
        m_x.push_back(lattice->moleFraction(k));
        tmpV_.push_back(0.0);
    }
}

void LatticeSolidPhase::setLatticeStoichiometry(const compositionMap& comp)
{
    for (size_t i = 0; i < m_lattice.size(); i++) {
        theta_[i] = getValue(comp, m_lattice[i]->name(), 0.0);
    }
    // Add in the lattice stoichiometry constraint
    for (size_t i = 1; i < m_lattice.size(); i++) {
        string econ = fmt::format("LC_{}_{}", i, name());
        size_t m = addElement(econ, 0.0, 0, 0.0, CT_ELEM_TYPE_LATTICERATIO);
        size_t mm = nElements();
        for (size_t k = 0; k < m_lattice[0]->nSpecies(); k++) {
            m_speciesComp[k * mm + m] = -theta_[0];
        }
        for (size_t k = 0; k < m_lattice[i]->nSpecies(); k++) {
            size_t ks = lkstart_[i] + k;
            m_speciesComp[ks * mm + m] = theta_[i];
        }
    }
}

void LatticeSolidPhase::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        getMoleFractions(m_x.data());
        size_t strt = 0;
        for (size_t n = 0; n < m_lattice.size(); n++) {
            m_lattice[n]->setTemperature(tnow);
            m_lattice[n]->setMoleFractions(&m_x[strt]);
            m_lattice[n]->setPressure(m_press);
            strt += m_lattice[n]->nSpecies();
        }
        m_tlast = tnow;
    }
}

void LatticeSolidPhase::setLatticeMoleFractionsByName(int nn, const std::string& x)
{
    m_lattice[nn]->setMoleFractionsByName(x);
    size_t loc = 0;
    for (size_t n = 0; n < m_lattice.size(); n++) {
        size_t nsp = m_lattice[n]->nSpecies();
        double ndens = m_lattice[n]->molarDensity();
        for (size_t k = 0; k < nsp; k++) {
            m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
            loc++;
        }
    }
    setMoleFractions(m_x.data());
}

void LatticeSolidPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","LatticeSolid");
    XML_Node& la = eosdata.child("LatticeArray");
    std::vector<XML_Node*> lattices = la.getChildren("phase");
    for (auto lattice : lattices) {
        addLattice(shared_ptr<ThermoPhase>(newPhase(*lattice)));
    }
    setLatticeStoichiometry(parseCompString(eosdata.child("LatticeStoichiometry").value()));
}

void LatticeSolidPhase::modifyOneHf298SS(const size_t k, const doublereal Hf298New)
{
    for (size_t n = 0; n < m_lattice.size(); n++) {
        if (lkstart_[n+1] < k) {
            size_t kk = k-lkstart_[n];
            MultiSpeciesThermo& l_spthermo = m_lattice[n]->speciesThermo();
            l_spthermo.modifyOneHf298(kk, Hf298New);
        }
    }
    invalidateCache();
    _updateThermo();
}

void LatticeSolidPhase::resetHf298(const size_t k)
{
    if (k != npos) {
        for (size_t n = 0; n < m_lattice.size(); n++) {
            if (lkstart_[n+1] < k) {
                size_t kk = k-lkstart_[n];
                m_lattice[n]->speciesThermo().resetHf298(kk);
            }
        }
    } else {
        for (size_t n = 0; n < m_lattice.size(); n++) {
            m_lattice[n]->speciesThermo().resetHf298(npos);
        }
    }
    invalidateCache();
    _updateThermo();
}

} // End namespace Cantera
