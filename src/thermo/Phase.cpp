/**
 *  @file Phase.cpp
 *   Definition file for class Phase.
 */

// Copyright 2001  California Institute of Technology

#include "cantera/thermo/Phase.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

Phase::Phase() :
    m_kk(0),
    m_ndim(3),
    m_xml(new XML_Node("phase")),
    m_id("<phase>"),
    m_name(""),
    m_temp(0.0),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1),
    m_speciesFrozen(false) ,
    m_Elements(new Elements())
{
    m_Elements->subscribe();
}

Phase::Phase(const Phase& right) :
    m_kk(0),
    m_ndim(3),
    m_xml(0),
    m_id("<phase>"),
    m_name(""),
    m_temp(0.0),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1),
    m_speciesFrozen(false) ,
    m_Elements(0)
{
    // Use the assignment operator to do the actual copying
    *this = operator=(right);
}

Phase& Phase::operator=(const Phase& right)
{
    // Check for self assignment.
    if (this == &right) {
        return *this;
    }

    // Handle our own data
    m_kk = right.m_kk;
    m_ndim = right.m_ndim;
    m_temp = right.m_temp;
    m_dens = right.m_dens;
    m_mmw = right.m_mmw;
    m_ym = right.m_ym;
    m_y = right.m_y;
    m_molwts = right.m_molwts;
    m_rmolwts = right.m_rmolwts;
    m_stateNum = -1;
    m_speciesFrozen = right.m_speciesFrozen;
    if (m_Elements) {
        int nleft = m_Elements->unsubscribe();
        if (nleft <= 0) {
            vector<Elements*>::iterator it;
            for (it  = Elements::Global_Elements_List.begin();
                    it != Elements::Global_Elements_List.end(); ++it) {
                if (*it == m_Elements) {
                    Elements::Global_Elements_List.erase(it);
                    break;
                }
            }
            delete m_Elements;
        }
    }
    m_Elements = right.m_Elements;
    if (m_Elements) {
        m_Elements->subscribe();
    }
    m_speciesNames = right.m_speciesNames;
    m_speciesComp = right.m_speciesComp;
    m_speciesCharge = right.m_speciesCharge;
    m_speciesSize = right.m_speciesSize;

    /*
     * This is a little complicated. -> Because we delete m_xml
     * in the destructor, we own m_xml completely, and we need
     * to have our own individual copies of the XML data tree
     * in each object
     */
    if (m_xml) {
        delete m_xml;
        m_xml = 0;
    }
    if (right.m_xml) {
        m_xml = new XML_Node();
        (right.m_xml)->copy(m_xml);
    }
    m_id    = right.m_id;
    m_name  = right.m_name;

    return *this;
}

// Destructor.
Phase::~Phase()
{
    if (m_xml) {
        delete m_xml;
        m_xml = 0;
    }

    int ileft = m_Elements->unsubscribe();
    /*
     * Here we may delete Elements Objects or not. Right now, we
     * will delete them. We also delete the global pointer entry
     * to keep everything consistent.
     */
    if (ileft <= 0) {
        vector<Elements*>::iterator it;
        for (it  = Elements::Global_Elements_List.begin();
                it != Elements::Global_Elements_List.end(); ++it) {
            if (*it == m_Elements) {
                Elements::Global_Elements_List.erase(it);
                break;
            }
        }
        delete m_Elements;
    }
}

inline void Phase::stateMFChangeCalc(bool forcerChange)
{
    // Right now we assume that the mole fractions have changed every time
    // the function is called
    m_stateNum++;
    if (m_stateNum > 1000000) {
        m_stateNum = -10000000;
    }
}

XML_Node& Phase::xml()
{
    return *m_xml;
}

std::string Phase::id() const
{
    return m_id;
}

void Phase::setID(std::string id)
{
    m_id = id;
}

std::string Phase::name() const
{
    return m_name;
}

void Phase::setName(std::string nm)
{
    m_name = nm;
}

size_t Phase::nElements() const
{
    return m_Elements->nElements();
}

string Phase::elementName(size_t m) const
{
    return m_Elements->elementName(m);
}

size_t Phase::elementIndex(std::string name) const
{
    return m_Elements->elementIndex(name);
}

const vector<string>& Phase::elementNames() const
{
    return m_Elements->elementNames();
}

doublereal Phase::atomicWeight(size_t m) const
{
    return m_Elements->atomicWeight(m);
}

doublereal Phase::entropyElement298(size_t m) const
{
    return m_Elements->entropyElement298(m);
}

const vector_fp& Phase::atomicWeights() const
{
    return m_Elements->atomicWeights();
}

int Phase::atomicNumber(size_t m) const
{
    return m_Elements->atomicNumber(m);
}

int Phase::elementType(size_t m) const
{
    return m_Elements->elementType(m);
}

doublereal Phase::nAtoms(size_t k, size_t m) const
{
    const size_t m_mm = m_Elements->nElements();
    if (m >= m_mm) {
        throw IndexError("Phase::nAtoms", "", m, nElements());
    }
    if (k >= nSpecies()) {
        throw IndexError("Phase::nAtoms", "", k, nSpecies());
    }
    return m_speciesComp[m_mm * k + m];
}

void Phase::getAtoms(size_t k, double* atomArray) const
{
    const size_t m_mm = m_Elements->nElements();
    for (size_t m = 0; m < m_mm; m++) {
        atomArray[m] = (double) m_speciesComp[m_mm * k + m];
    }
}

size_t Phase::speciesIndex(std::string nameStr) const
{
    std::string pn;
    std::string sn = parseSpeciesName(nameStr, pn);
    if (pn == "" || pn == m_name || pn == m_id) {
        vector<string>::const_iterator it = m_speciesNames.begin();
        for (size_t k = 0; k < m_kk; k++) {
            if (*it == sn) {
                return k;
            }
            ++it;
        }
        return npos;
    }
    return npos;
}

string Phase::speciesName(size_t k) const
{
    if (k >= nSpecies())
        throw IndexError("Phase::speciesName", "m_speciesNames", k, nSpecies());
    return m_speciesNames[k];
}

const vector<string>& Phase::speciesNames() const
{
    return m_speciesNames;
}

std::string Phase::speciesSPName(int k) const
{
    std::string sn = speciesName(k);
    return(m_name + ":" + sn);
}

void Phase::saveState(vector_fp& state) const
{
    state.resize(nSpecies() + 2);
    saveState(state.size(),&(state[0]));
}
void Phase::saveState(size_t lenstate, doublereal* state) const
{
    state[0] = temperature();
    state[1] = density();
    getMassFractions(state + 2);
}

void Phase::restoreState(const vector_fp& state)
{
    restoreState(state.size(),&state[0]);
}

void Phase::restoreState(size_t lenstate, const doublereal* state)
{
    if (lenstate >= nSpecies() + 2) {
        setMassFractions_NoNorm(state + 2);
        setTemperature(state[0]);
        setDensity(state[1]);
    } else {
        throw ArraySizeError("Phase::restoreState",
                             lenstate,nSpecies()+2);
    }
}

void Phase::setMoleFractions(const doublereal* const x)
{
    doublereal sum = dot(x, x + m_kk, m_molwts.begin());
    doublereal rsum = 1.0/sum;
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(rsum));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(),
              m_y.begin(), multiplies<double>());
    doublereal norm = accumulate(x, x + m_kk, 0.0);
    m_mmw = sum/norm;

    // Call a routine to determine whether state has changed.
    stateMFChangeCalc();
}

void Phase::setMoleFractions_NoNorm(const doublereal* const x)
{
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    doublereal rmmw = 1.0/m_mmw;
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(rmmw));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(),
              m_y.begin(), multiplies<double>());

    // Call a routine to determine whether state has changed.
    stateMFChangeCalc();
}

void Phase::setMoleFractionsByName(compositionMap& xMap)
{
    size_t kk = nSpecies();
    doublereal x;
    vector_fp mf(kk, 0.0);
    for (size_t k = 0; k < kk; k++) {
        x = xMap[speciesName(k)];
        if (x > 0.0) {
            mf[k] = x;
        }
    }
    setMoleFractions(&mf[0]);
}

void Phase::setMoleFractionsByName(const std::string& x)
{
    size_t kk = nSpecies();
    compositionMap xx;
    for (size_t k = 0; k < kk; k++) {
        xx[speciesName(k)] = -1.0;
    }
    parseCompString(x, xx);
    setMoleFractionsByName(xx);
}

void Phase::setMassFractions(const doublereal* const y)
{
    doublereal norm = 0.0, sum = 0.0;
    norm = accumulate(y, y + m_kk, 0.0);
    copy(y, y + m_kk, m_y.begin());
    scale(y, y + m_kk, m_y.begin(), 1.0/norm);

    transform(m_y.begin(), m_y.begin() + m_kk, m_rmolwts.begin(),
              m_ym.begin(), multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.begin() + m_kk, 0.0);
    m_mmw = 1.0/sum;

    // Call a routine to determine whether state has changed.
    stateMFChangeCalc();
}

void Phase::setMassFractions_NoNorm(const doublereal* const y)
{
    doublereal sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
              multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    m_mmw = 1.0/sum;

    // Call a routine to determine whether state has changed.
    stateMFChangeCalc();
}

void Phase::setMassFractionsByName(compositionMap& yMap)
{
    size_t kk = nSpecies();
    doublereal y;
    vector_fp mf(kk, 0.0);
    for (size_t k = 0; k < kk; k++) {
        y = yMap[speciesName(k)];
        if (y > 0.0) {
            mf[k] = y;
        }
    }
    setMassFractions(&mf[0]);
}

void Phase::setMassFractionsByName(const std::string& y)
{
    size_t kk = nSpecies();
    compositionMap yy;
    for (size_t k = 0; k < kk; k++) {
        yy[speciesName(k)] = -1.0;
    }
    parseCompString(y, yy);
    setMassFractionsByName(yy);
}

void Phase::setState_TRX(doublereal t, doublereal dens, const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setDensity(dens);
}

void Phase::setState_TNX(doublereal t, doublereal n, const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setMolarDensity(n);
}

void Phase::setState_TRX(doublereal t, doublereal dens, compositionMap& x)
{
    setMoleFractionsByName(x);
    setTemperature(t);
    setDensity(dens);
}

void Phase::setState_TRY(doublereal t, doublereal dens, const doublereal* y)
{
    setMassFractions(y);
    setTemperature(t);
    setDensity(dens);
}

void Phase::setState_TRY(doublereal t, doublereal dens, compositionMap& y)
{
    setMassFractionsByName(y);
    setTemperature(t);
    setDensity(dens);
}

void Phase::setState_TR(doublereal t, doublereal rho)
{
    setTemperature(t);
    setDensity(rho);
}

void Phase::setState_TX(doublereal t, doublereal* x)
{
    setTemperature(t);
    setMoleFractions(x);
}

void Phase::setState_TY(doublereal t, doublereal* y)
{
    setTemperature(t);
    setMassFractions(y);
}

void Phase::setState_RX(doublereal rho, doublereal* x)
{
    setMoleFractions(x);
    setDensity(rho);
}

void Phase::setState_RY(doublereal rho, doublereal* y)
{
    setMassFractions(y);
    setDensity(rho);
}

doublereal Phase::molecularWeight(size_t k) const
{
    if (k >= nSpecies()) {
        throw IndexError("Phase::molecularWeight", "m_weight", k, nSpecies());
    }
    return m_molwts[k];
}

void Phase::getMolecularWeights(vector_fp& weights) const
{
    const vector_fp& mw = molecularWeights();
    if (weights.size() < mw.size()) {
        weights.resize(mw.size());
    }
    copy(mw.begin(), mw.end(), weights.begin());
}

void Phase::getMolecularWeights(int iwt, doublereal* weights) const
{
    const vector_fp& mw = molecularWeights();
    copy(mw.begin(), mw.end(), weights);
}

void Phase::getMolecularWeights(doublereal* weights) const
{
    const vector_fp& mw = molecularWeights();
    copy(mw.begin(), mw.end(), weights);
}

const vector_fp& Phase::molecularWeights() const
{
    return m_molwts;
}

void Phase::getMoleFractionsByName(compositionMap& x) const
{
    x.clear();
    size_t kk = nSpecies();
    for (size_t k = 0; k < kk; k++) {
        x[speciesName(k)] = Phase::moleFraction(k);
    }
}

void Phase::getMoleFractions(doublereal* const x) const
{
    scale(m_ym.begin(), m_ym.end(), x, m_mmw);
}

doublereal Phase::moleFraction(size_t k) const
{
    if (k < m_kk) {
        return m_ym[k] * m_mmw;
    } else {
        throw CanteraError("Phase::moleFraction",
                           "illegal species index number");
    }
    return 0.0;
}

doublereal Phase::moleFraction(std::string nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return moleFraction(iloc);
    } else {
        return 0.0;
    }
}

const doublereal* Phase::moleFractdivMMW() const
{
    return &m_ym[0];
}

doublereal Phase::massFraction(size_t k) const
{
    if (k < m_kk) {
        return m_y[k];
    }
    throw CanteraError("State:massFraction", "illegal species index number");
    return 0.0;
}

doublereal Phase::massFraction(std::string nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return massFractions()[iloc];
    } else {
        return 0.0;
    }
}

void Phase::getMassFractions(doublereal* const y) const
{
    copy(m_y.begin(), m_y.end(), y);
}

doublereal Phase::concentration(const size_t k) const
{
    if (k < m_kk) {
        return m_y[k] * m_dens * m_rmolwts[k] ;
    }
    throw CanteraError("State:massFraction", "illegal species index number");
    return 0.0;
}

void Phase::getConcentrations(doublereal* const c) const
{
    scale(m_ym.begin(), m_ym.end(), c, m_dens);
}

void Phase::setConcentrations(const doublereal* const conc)
{
    doublereal sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        sum += conc[k]*m_molwts[k];
        norm += conc[k];
    }
    m_mmw = sum/norm;
    setDensity(sum);
    doublereal rsum = 1.0/sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = conc[k] * rsum;
        m_y[k] =  m_ym[k] * m_molwts[k];
    }

    // Call a routine to determine whether state has changed.
    stateMFChangeCalc();
}

doublereal Phase::molarDensity() const
{
    return density()/meanMolecularWeight();
}

void Phase::setMolarDensity(const doublereal molarDensity)
{
    m_dens = molarDensity*meanMolecularWeight();
}

doublereal Phase::molarVolume() const
{
    return 1.0/molarDensity();
}

doublereal Phase::charge(size_t k) const
{
    return m_speciesCharge[k];
}

doublereal Phase::chargeDensity() const
{
    size_t kk = nSpecies();
    doublereal cdens = 0.0;
    for (size_t k = 0; k < kk; k++) {
        cdens += charge(k)*moleFraction(k);
    }
    cdens *= Faraday;
    return cdens;
}

doublereal Phase::mean_X(const doublereal* const Q) const
{
    return m_mmw*std::inner_product(m_ym.begin(), m_ym.end(), Q, 0.0);
}

doublereal Phase::mean_Y(const doublereal* const Q) const
{
    return dot(m_y.begin(), m_y.end(), Q);
}

doublereal Phase::sum_xlogx() const
{
    return m_mmw* Cantera::sum_xlogx(m_ym.begin(), m_ym.end()) + log(m_mmw);
}

doublereal Phase::sum_xlogQ(doublereal* Q) const
{
    return m_mmw * Cantera::sum_xlogQ(m_ym.begin(), m_ym.end(), Q);
}

void Phase::addElement(const std::string& symbol, doublereal weight)
{
    m_Elements->addElement(symbol, weight);
}

void Phase::addElement(const XML_Node& e)
{
    m_Elements->addElement(e);
}

void Phase::addUniqueElement(const std::string& symbol, doublereal weight,
                             int atomicNumber, doublereal entropy298,
                             int elem_type)
{
    m_Elements->addUniqueElement(symbol, weight, atomicNumber,
                                 entropy298, elem_type);
}

void Phase::addUniqueElement(const XML_Node& e)
{
    m_Elements->addUniqueElement(e);
}

void Phase::addElementsFromXML(const XML_Node& phase)
{
    m_Elements->addElementsFromXML(phase);
}

void Phase::freezeElements()
{
    m_Elements->freezeElements();
}

bool Phase::elementsFrozen()
{
    return m_Elements->elementsFrozen();
}

size_t Phase::addUniqueElementAfterFreeze(const std::string& symbol,
        doublereal weight, int atomicNumber,
        doublereal entropy298, int elem_type)
{
    size_t ii = elementIndex(symbol);
    if (ii != npos) {
        return ii;
    }
    // Check to see that the element isn't really in the list
    m_Elements->m_elementsFrozen = false;
    addUniqueElement(symbol, weight, atomicNumber, entropy298, elem_type);
    m_Elements->m_elementsFrozen = true;
    size_t m_mm = m_Elements->nElements();
    ii = elementIndex(symbol);
    if (ii != m_mm-1) {
        throw CanteraError("Phase::addElementAfterFreeze()", "confused");
    }
    if (m_kk > 0) {
        vector_fp old(m_speciesComp);
        m_speciesComp.resize(m_kk*m_mm, 0.0);
        for (size_t k = 0; k < m_kk; k++) {
            size_t m_old = m_mm - 1;
            for (size_t m = 0; m < m_old; m++) {
                m_speciesComp[k * m_mm + m] =  old[k * (m_old) + m];
            }
            m_speciesComp[k * (m_mm) + (m_mm-1)] = 0.0;
        }
    }
    return ii;
}

void Phase::addSpecies(const std::string& name, const doublereal* comp,
                       doublereal charge, doublereal size)
{
    m_Elements->freezeElements();
    m_speciesNames.push_back(name);
    m_speciesCharge.push_back(charge);
    m_speciesSize.push_back(size);
    size_t ne = m_Elements->nElements();
    // Create a changeable copy of the element composition. We now change
    // the charge potentially
    vector_fp compNew(ne);
    for (size_t m = 0; m < ne; m++) {
        compNew[m] = comp[m];
    }
    double wt = 0.0;
    const vector_fp& aw = m_Elements->atomicWeights();
    if (charge != 0.0) {
        size_t eindex = m_Elements->elementIndex("E");
        if (eindex != npos) {
            doublereal ecomp = compNew[eindex];
            if (fabs(charge + ecomp) > 0.001) {
                if (ecomp != 0.0) {
                    throw CanteraError("Phase::addSpecies",
                        "Input charge and element E compositions differ "
                        "for species " + name);
                } else {
                    // Just fix up the element E composition based on the input
                    // species charge
                    compNew[eindex] = -charge;
                }
            }
        } else {
            addUniqueElementAfterFreeze("E", 0.000545, 0, 0.0,
                                        CT_ELEM_TYPE_ELECTRONCHARGE);
            ne = m_Elements->nElements();
            eindex = m_Elements->elementIndex("E");
            compNew.resize(ne);
            compNew[ne - 1] = - charge;
        }
    }
    for (size_t m = 0; m < ne; m++) {
        m_speciesComp.push_back(compNew[m]);
        wt += compNew[m] * aw[m];
    }
    m_molwts.push_back(wt);
    m_kk++;
}

void Phase::addUniqueSpecies(const std::string& name, const doublereal* comp,
                             doublereal charge, doublereal size)
{
    vector<string>::const_iterator it = m_speciesNames.begin();
    for (size_t k = 0; k < m_kk; k++) {
        if (*it == name) {
            // We have found a match. At this point we could do some
            // compatibility checks. However, let's just return for the moment
            // without specifying any error.
            size_t m_mm = m_Elements->nElements();
            for (size_t i = 0; i < m_mm; i++) {
                if (comp[i] != m_speciesComp[m_kk * m_mm + i]) {
                    throw CanteraError("addUniqueSpecies",
                                       "Duplicate species have different "
                                       "compositions: " + *it);
                }
            }
            if (charge != m_speciesCharge[m_kk]) {
                throw CanteraError("addUniqueSpecies",
                                   "Duplicate species have different "
                                   "charges: " + *it);
            }
            if (size != m_speciesSize[m_kk]) {
                throw CanteraError("addUniqueSpecies",
                                   "Duplicate species have different "
                                   "sizes: " + *it);
            }
            return;
        }
        ++it;
    }
    addSpecies(name, comp, charge, size);
}

void Phase::freezeSpecies()
{
    m_speciesFrozen = true;
    init(molecularWeights());
    size_t kk = nSpecies();
    size_t nv = kk + 2;
    m_kk = nSpecies();
}

void Phase::init(const vector_fp& mw)
{
    m_kk = mw.size();
    m_rmolwts.resize(m_kk);
    m_y.resize(m_kk, 0.0);
    m_ym.resize(m_kk, 0.0);
    copy(mw.begin(), mw.end(), m_molwts.begin());
    for (size_t k = 0; k < m_kk; k++) {
        if (m_molwts[k] < 0.0) {
            throw CanteraError("Phase::init",
                               "negative molecular weight for species number "
                               + int2str(k));
        }

        // Some surface phases may define species representing empty sites
        // hat have zero molecular weight. Give them a very small molecular
        // weight to avoid dividing by zero.
        if (m_molwts[k] < Tiny) {
            m_molwts[k] = Tiny;
        }
        m_rmolwts[k] = 1.0/m_molwts[k];
    }

    // Now that we have resized the State object, let's fill it with a valid
    // mass fraction vector that sums to one. The Phase object should never
    // have a mass fraction vector that doesn't sum to one. We will assume that
    // species 0 has a mass fraction of 1.0 and mass fraction of all other
    // species is 0.0.
    m_y[0] = 1.0;
    m_ym[0] = m_y[0] * m_rmolwts[0];
    m_mmw = 1.0 / m_ym[0];
}

bool Phase::ready() const
{
    return (m_kk > 0 && m_Elements->elementsFrozen() && m_speciesFrozen);
}

} // namespace Cantera
