/**
 *  @file LatticeSolidPhase.cpp
 *  Definitions for a simple thermodynamics model of a bulk solid phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).
 */

#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
LatticeSolidPhase::LatticeSolidPhase() :
    m_press(-1.0),
    m_molar_density(0.0),
    m_nlattice(0)
{
}

LatticeSolidPhase::LatticeSolidPhase(const LatticeSolidPhase& right) :
    m_press(-1.0),
    m_molar_density(0.0),
    m_nlattice(0)
{
    *this = right;
}

LatticeSolidPhase&
LatticeSolidPhase::operator=(const LatticeSolidPhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_tlast = right.m_tlast;
        m_press = right.m_press;
        m_molar_density = right.m_molar_density;
        m_nlattice = right.m_nlattice;
        deepStdVectorPointerCopy<LatticePhase>(right.m_lattice, m_lattice);
        m_x = right.m_x;
        theta_ = right.theta_;
        tmpV_ = right.tmpV_;
    }
    return *this;
}

LatticeSolidPhase::~LatticeSolidPhase()
{
    // We own the sublattices. So we have to delete the sublattices
    for (size_t n = 0; n < m_nlattice; n++) {
        delete m_lattice[n];
        m_lattice[n] = 0;
    }
}

ThermoPhase* LatticeSolidPhase::duplMyselfAsThermoPhase() const
{
    return new LatticeSolidPhase(*this);
}

doublereal LatticeSolidPhase::minTemp(size_t k) const
{
    if (k != npos) {
        for (size_t n = 0; n < m_nlattice; n++) {
            if (lkstart_[n+1] < k) {
                return (m_lattice[n])->minTemp(k-lkstart_[n]);
            }
        }
    }
    doublereal mm = 1.0E300;
    for (size_t n = 0; n < m_nlattice; n++) {
        double ml = (m_lattice[n])->minTemp();
        mm = std::min(mm, ml);
    }
    return mm;
}

doublereal LatticeSolidPhase::maxTemp(size_t k) const
{
    if (k != npos) {
        for (size_t n = 0; n < m_nlattice; n++) {
            if (lkstart_[n+1] < k) {
                return (m_lattice[n])->maxTemp(k - lkstart_[n]);
            }
        }
    }
    doublereal mm = -1.0E300;
    for (size_t n = 0; n < m_nlattice; n++) {
        double ml = (m_lattice[n])->maxTemp();
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
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->enthalpy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::intEnergy_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->intEnergy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::entropy_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->entropy_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::gibbs_mole() const
{
    _updateThermo();
    doublereal  sum = 0.0;
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->gibbs_mole();
    }
    return sum;
}

doublereal LatticeSolidPhase::cp_mole() const
{
    _updateThermo();
    doublereal sum = 0.0;
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->cp_mole();
    }
    return sum;
}

void LatticeSolidPhase::getActivityConcentrations(doublereal* c) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
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

void  LatticeSolidPhase::setPressure(doublereal p)
{
    m_press = p;
    for (size_t n = 0; n < m_nlattice; n++) {
        m_lattice[n]->setPressure(m_press);
    }
    calcDensity();
}

doublereal  LatticeSolidPhase::calcDensity()
{
    double sum = 0.0;
    for (size_t n = 0; n < m_nlattice; n++) {
        sum += theta_[n] * m_lattice[n]->density();
    }
    Phase::setDensity(sum);
    return sum;
}

void LatticeSolidPhase::setMoleFractions(const doublereal* const x)
{
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->setMoleFractions(x + strt);
        strt += nsp;
    }
    for (size_t k = 0; k < strt; k++) {
        m_x[k] = x[k] / m_nlattice;
    }
    Phase::setMoleFractions(DATA_PTR(m_x));
    calcDensity();
}

void LatticeSolidPhase::getMoleFractions(doublereal* const x) const
{
    size_t strt = 0;
    // the ifdef block should be the way we calculate this.!!!!!
    Phase::getMoleFractions(x);
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nsp =  m_lattice[n]->nSpecies();
        double sum = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            sum += (x + strt)[k];
        }
        for (size_t k = 0; k < nsp; k++) {
            (x + strt)[k] /= sum;
        }
        /*
         * At this point we can check against the mole fraction vector of the underlying LatticePhase objects and
         * get the same answer.
         */
        if (DEBUG_MODE_ENABLED) {
            m_lattice[n]->getMoleFractions(&(m_x[strt]));
            for (size_t k = 0; k < nsp; k++) {
                if (fabs((x + strt)[k] - m_x[strt+k]) > 1.0E-14) {
                    throw CanteraError("LatticeSolidPhase::getMoleFractions()",
                                       "internal error");
                }
            }
        }
        strt += nsp;
    }
}

void LatticeSolidPhase::getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nlsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->getChemPotentials(mu+strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nlsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarEnthalpies(hbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nlsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarEntropies(sbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarCp(doublereal* cpbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nlsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarCp(cpbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nlsp =  m_lattice[n]->nSpecies();
        m_lattice[n]->getPartialMolarVolumes(vbar + strt);
        strt += nlsp;
    }
}

void LatticeSolidPhase::getStandardChemPotentials(doublereal* mu0) const
{
    _updateThermo();
    size_t strt = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        m_lattice[n]->getStandardChemPotentials(mu0+strt);
        strt += m_lattice[n]->nSpecies();
    }
}

void LatticeSolidPhase::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    for (size_t n = 0; n < m_nlattice; n++) {
        m_lattice[n]->getGibbs_RT_ref(grt + lkstart_[n]);
    }
}

void LatticeSolidPhase::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= GasConstant * temperature();
    }
}

void LatticeSolidPhase::installSlavePhases(Cantera::XML_Node* phaseNode)
{
    size_t kk = 0;
    size_t kstart = 0;
    m_speciesData.clear();

    XML_Node& la = phaseNode->child("thermo").child("LatticeArray");
    std::vector<XML_Node*> lattices = la.getChildren("phase");
    for (size_t n = 0; n < m_nlattice; n++) {
        LatticePhase* lp = m_lattice[n];
        vector<doublereal> constArr(lp->nElements());
        const vector_fp& aws = lp->atomicWeights();
        for (size_t es = 0; es < lp->nElements(); es++) {
            addElement(lp->elementName(es), aws[es], lp->atomicNumber(es),
                       lp->entropyElement298(es), lp->elementType(es));
        }
        kstart = kk;

        for (size_t k = 0; k < lp->nSpecies(); k++) {
            addSpecies(lp->species(k));
            kk++;
        }
        /*
         *  Add in the lattice stoichiometry constraint
         */
        if (n > 0) {
            string econ = "LC_" + int2str(n) + "_" + id();
            size_t m = addElement(econ, 0.0, 0, 0.0, CT_ELEM_TYPE_LATTICERATIO);
            size_t mm = nElements();
            size_t nsp0 =  m_lattice[0]->nSpecies();
            for (size_t k = 0; k < nsp0; k++) {
                m_speciesComp[k * mm + m] = -theta_[0];
            }
            for (size_t k = 0; k < lp->nSpecies(); k++) {
                size_t ks = kstart + k;
                m_speciesComp[ks * mm + m] = theta_[n];
            }
        }
    }
}

void LatticeSolidPhase::initThermo()
{
    initLengths();
    size_t loc = 0;
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nsp = m_lattice[n]->nSpecies();
        lkstart_[n] = loc;
        for (size_t k = 0; k < nsp; k++) {
            m_x[loc] =m_lattice[n]->moleFraction(k) / (double) m_nlattice;
            loc++;
        }
        lkstart_[n+1] = loc;
    }
    setMoleFractions(DATA_PTR(m_x));
    ThermoPhase::initThermo();
}

void LatticeSolidPhase::initLengths()
{
    theta_.resize(m_nlattice,0);
    lkstart_.resize(m_nlattice+1);
    m_x.resize(m_kk, 0.0);
    tmpV_.resize(m_kk, 0.0);
}

void LatticeSolidPhase::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        getMoleFractions(DATA_PTR(m_x));
        size_t strt = 0;
        for (size_t n = 0; n < m_nlattice; n++) {
            m_lattice[n]->setTemperature(tnow);
            m_lattice[n]->setMoleFractions(DATA_PTR(m_x) + strt);
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
    for (size_t n = 0; n < m_nlattice; n++) {
        size_t nsp = m_lattice[n]->nSpecies();
        double ndens = m_lattice[n]->molarDensity();
        for (size_t k = 0; k < nsp; k++) {
            m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
            loc++;
        }
    }
    setMoleFractions(DATA_PTR(m_x));
}

void LatticeSolidPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","LatticeSolid");
    XML_Node& la = eosdata.child("LatticeArray");
    std::vector<XML_Node*> lattices = la.getChildren("phase");
    m_nlattice = lattices.size();
    for (size_t n = 0; n < m_nlattice; n++) {
        m_lattice.push_back((LatticePhase*)newPhase(*lattices[n]));
    }
    std::vector<string> pnam;
    std::vector<string> pval;
    int np = getPairs(eosdata.child("LatticeStoichiometry"), pnam, pval);
    theta_.resize(m_nlattice);
    for (int i = 0; i < np; i++) {
        double val = fpValueCheck(pval[i]);
        bool found = false;
        for (size_t j = 0; j < m_nlattice; j++) {
            ThermoPhase& tp = *(m_lattice[j]);
            string idj = tp.id();
            if (idj == pnam[i]) {
                theta_[j] = val;
                found = true;
                break;
            }
        }
        if (!found) {
            throw CanteraError("", "not found");
        }
    }

}

void LatticeSolidPhase::modifyOneHf298SS(const size_t k, const doublereal Hf298New)
{
    for (size_t n = 0; n < m_nlattice; n++) {
        if (lkstart_[n+1] < k) {
            size_t kk = k-lkstart_[n];
            SpeciesThermo& l_spthermo =  m_lattice[n]->speciesThermo();
            l_spthermo.modifyOneHf298(kk, Hf298New);
        }
    }
    m_tlast += 0.0001234;
    _updateThermo();
}

} // End namespace Cantera
