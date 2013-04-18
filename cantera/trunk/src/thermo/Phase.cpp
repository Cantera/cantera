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
    m_temp(0.001),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1),
    m_speciesFrozen(false),
    m_elementsFrozen(false),
    m_mm(0),
    m_elem_type(0)
{
}

Phase::Phase(const Phase& right) :
    m_kk(0),
    m_ndim(3),
    m_xml(0),
    m_id("<phase>"),
    m_name(""),
    m_temp(0.001),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1),
    m_speciesFrozen(false) ,
    m_elementsFrozen(false),
    m_mm(0),
    m_elem_type(0)
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
    m_speciesNames = right.m_speciesNames;
    m_speciesComp = right.m_speciesComp;
    m_speciesCharge = right.m_speciesCharge;
    m_speciesSize = right.m_speciesSize;

    m_mm             = right.m_mm;
    m_elementsFrozen = right.m_elementsFrozen;
    m_atomicWeights  = right.m_atomicWeights;
    m_atomicNumbers  = right.m_atomicNumbers;
    m_elementNames   = right.m_elementNames;
    m_entropy298     = right.m_entropy298;
    m_elem_type      = right.m_elem_type;
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

Phase::~Phase()
{
    if (m_xml) {
        delete m_xml;
        m_xml = 0;
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

void Phase::setID(const std::string& id_)
{
    m_id = id_;
}

std::string Phase::name() const
{
    return m_name;
}

void Phase::setName(const std::string& nm)
{
    m_name = nm;
}

size_t Phase::nElements() const
{
    return m_mm;
}

void Phase::checkElementIndex(size_t m) const
{
    if (m >= m_mm) {
        throw IndexError("checkElementIndex", "elements", m, m_mm-1);
    }
}

void Phase::checkElementArraySize(size_t mm) const
{
    if (m_mm > mm) {
        throw ArraySizeError("checkElementArraySize", mm, m_mm);
    }
}

string Phase::elementName(size_t m) const
{
    checkElementIndex(m);
    return m_elementNames[m];
}

size_t Phase::elementIndex(const std::string& elementName) const
{
    for (size_t i = 0; i < m_mm; i++) {
        if (m_elementNames[i] == elementName) {
            return i;
        }
    }
    return npos;
}

const vector<string>& Phase::elementNames() const
{
    return m_elementNames;
}

doublereal Phase::atomicWeight(size_t m) const
{
    return m_atomicWeights[m];
}

doublereal Phase::entropyElement298(size_t m) const
{
    AssertThrowMsg(m_entropy298[m] != ENTROPY298_UNKNOWN,
                   "Elements::entropy298",
                   "Entropy at 298 K of element is unknown");
    AssertTrace(m < m_mm);
    return m_entropy298[m];
}

const vector_fp& Phase::atomicWeights() const
{
    return m_atomicWeights;
}

int Phase::atomicNumber(size_t m) const
{
    return m_atomicNumbers[m];
}

int Phase::elementType(size_t m) const
{
    return m_elem_type[m];
}

int Phase::changeElementType(int m, int elem_type)
{
    int old = m_elem_type[m];
    m_elem_type[m] = elem_type;
    return old;
}

doublereal Phase::nAtoms(size_t k, size_t m) const
{
    checkElementIndex(m);
    checkSpeciesIndex(k);
    return m_speciesComp[m_mm * k + m];
}

void Phase::getAtoms(size_t k, double* atomArray) const
{
    for (size_t m = 0; m < m_mm; m++) {
        atomArray[m] = (double) m_speciesComp[m_mm * k + m];
    }
}

size_t Phase::speciesIndex(const std::string& nameStr) const
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
    checkSpeciesIndex(k);
    return m_speciesNames[k];
}

const vector<string>& Phase::speciesNames() const
{
    return m_speciesNames;
}

void Phase::checkSpeciesIndex(size_t k) const
{
    if (k >= m_kk) {
        throw IndexError("checkSpeciesIndex", "species", k, m_kk-1);
    }
}

void Phase::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_kk);
    }
}

std::string Phase::speciesSPName(int k) const
{
    std::string sn = speciesName(k);
    return m_name + ":" + sn;
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
    // Use m_y as a temporary work vector for the non-negative mole fractions
    doublereal norm = 0.0;
    /*
     * sum is calculated below as the unnormalized molecular weight
     */
    doublereal sum = 0;
    for (size_t k = 0; k < m_kk; k++) {
        double xk = std::max(x[k], 0.0); // Ignore negative mole fractions
        m_y[k] = xk;
        norm += xk;
        sum += m_molwts[k] * xk;
    }
    /*
     * Set m_ym_ to the normalized mole fractions divided by the normalized mean molecular weight:
     *         m_ym_k = X_k / (sum_k X_k M_k)
     */
    // transform(m_y.begin(), m_y.end(), m_ym.begin(), timesConstant<double>(1.0/sum));
    const doublereal invSum = 1.0/sum;
    for (size_t k=0; k < m_kk; k++) {
        m_ym[k] = m_y[k]*invSum;
    }
    /*
     * Now set m_y to the normalized mass fractions
     *          m_y =  X_k M_k / (sum_k X_k M_k)
     */
    // transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(), m_y.begin(), multiplies<double>());
    for (size_t k=0; k < m_kk; k++) {
        m_y[k] = m_ym[k] * m_molwts[k];
    }
    /*
     * Calculate the normalized molecular weight
     */
    m_mmw = sum/norm;
    m_stateNum++;
}

void Phase::setMoleFractions_NoNorm(const doublereal* const x)
{
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    doublereal rmmw = 1.0/m_mmw;
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(rmmw));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(),
              m_y.begin(), multiplies<double>());
    m_stateNum++;
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
    compositionMap c = parseCompString(x, speciesNames());
    setMoleFractionsByName(c);
}

void Phase::setMassFractions(const doublereal* const y)
{
    for (size_t k = 0; k < m_kk; k++) {
        m_y[k] = std::max(y[k], 0.0); // Ignore negative mass fractions
    }
    doublereal norm = accumulate(m_y.begin(), m_y.end(), 0.0);
    scale(m_y.begin(), m_y.end(), m_y.begin(), 1.0/norm);

    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(),
              m_ym.begin(), multiplies<double>());
    m_mmw = 1.0 / accumulate(m_ym.begin(), m_ym.end(), 0.0);
    m_stateNum++;
}

void Phase::setMassFractions_NoNorm(const doublereal* const y)
{
    doublereal sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
              multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    m_mmw = 1.0/sum;
    m_stateNum++;
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
    compositionMap c = parseCompString(y, speciesNames());
    setMassFractionsByName(c);
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
    checkSpeciesIndex(k);
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
    checkSpeciesIndex(k);
    return m_ym[k] * m_mmw;
}

doublereal Phase::moleFraction(const std::string& nameSpec) const
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
    checkSpeciesIndex(k);
    return m_y[k];
}

doublereal Phase::massFraction(const std::string& nameSpec) const
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
    checkSpeciesIndex(k);
    return m_y[k] * m_dens * m_rmolwts[k] ;
}

void Phase::getConcentrations(doublereal* const c) const
{
    scale(m_ym.begin(), m_ym.end(), c, m_dens);
}

void Phase::setConcentrations(const doublereal* const conc)
{
    // Use m_y as temporary storage for non-negative concentrations
    doublereal sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        double ck = std::max(conc[k], 0.0); // Ignore negative concentrations
        m_y[k] = ck;
        sum += ck * m_molwts[k];
        norm += ck;
    }
    m_mmw = sum/norm;
    setDensity(sum);
    doublereal rsum = 1.0/sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = m_y[k] * rsum;
        m_y[k] = m_ym[k] * m_molwts[k]; // m_y is now the mass fraction
    }
    m_stateNum++;
}

doublereal Phase::molarDensity() const
{
    return density()/meanMolecularWeight();
}

void Phase::setMolarDensity(const doublereal molar_density)
{
    m_dens = molar_density*meanMolecularWeight();
}

doublereal Phase::molarVolume() const
{
    return 1.0/molarDensity();
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
    if (weight == -12345.0) {
        weight = LookupWtElements(symbol);
        if (weight < 0.0) {
            throw ElementsFrozen("addElement");
        }
    }
    if (m_elementsFrozen) {
        throw ElementsFrozen("addElement");
        return;
    }
    m_atomicWeights.push_back(weight);
    m_elementNames.push_back(symbol);
    if (symbol == "E") {
        m_elem_type.push_back(CT_ELEM_TYPE_ELECTRONCHARGE);
    } else {
        m_elem_type.push_back(CT_ELEM_TYPE_ABSPOS);
    }

    m_mm++;
}

void Phase::addElement(const XML_Node& e)
{
    doublereal weight = fpValue(e["atomicWt"]);
    string symbol = e["name"];
    addElement(symbol, weight);
}

void Phase::addUniqueElement(const std::string& symbol, doublereal weight,
                             int atomic_number, doublereal entropy298,
                             int elem_type)
{
    if (weight == -12345.0) {
        weight =  LookupWtElements(symbol);
        if (weight < 0.0) {
            throw ElementsFrozen("addElement");
        }
    }
    /*
     * First decide if this element has been previously added
     * by conducting a string search. If it unique, add it to
     * the list.
     */
    int ifound = 0;
    int i = 0;
    for (vector<string>::const_iterator it = m_elementNames.begin();
            it < m_elementNames.end(); ++it, ++i) {
        if (*it == symbol) {
            ifound = 1;
            break;
        }
    }
    if (!ifound) {
        if (m_elementsFrozen) {
            throw ElementsFrozen("addElement");
            return;
        }
        m_atomicWeights.push_back(weight);
        m_elementNames.push_back(symbol);
        m_atomicNumbers.push_back(atomic_number);
        m_entropy298.push_back(entropy298);
        if (symbol == "E") {
            m_elem_type.push_back(CT_ELEM_TYPE_ELECTRONCHARGE);
        } else {
            m_elem_type.push_back(elem_type);
        }
        m_mm++;
    } else {
        if (m_atomicWeights[i] != weight) {
            throw CanteraError("AddUniqueElement",
                               "Duplicate Elements (" + symbol + ") have different weights");
        }
    }
}

void Phase::addUniqueElement(const XML_Node& e)
{
    doublereal weight = 0.0;
    if (e.hasAttrib("atomicWt")) {
        weight = fpValue(stripws(e["atomicWt"]));
    }
    int anum = 0;
    if (e.hasAttrib("atomicNumber")) {
        anum = atoi(stripws(e["atomicNumber"]).c_str());
    }
    string symbol = e["name"];
    doublereal entropy298 = ENTROPY298_UNKNOWN;
    if (e.hasChild("entropy298")) {
        XML_Node& e298Node = e.child("entropy298");
        if (e298Node.hasAttrib("value")) {
            entropy298 = fpValueCheck(stripws(e298Node["value"]));
        }
    }
    if (weight != 0.0) {
        addUniqueElement(symbol, weight, anum, entropy298);
    } else {
        addUniqueElement(symbol);
    }
}

void Phase::addElementsFromXML(const XML_Node& phase)
{
    // get the declared element names
    if (! phase.hasChild("elementArray")) {
        throw CanteraError("Elements::addElementsFromXML",
                           "phase xml node doesn't have \"elementArray\" XML Node");
    }
    XML_Node& elements = phase.child("elementArray");
    vector<string> enames;
    ctml::getStringArray(elements, enames);

    // // element database defaults to elements.xml
    string element_database = "elements.xml";
    if (elements.hasAttrib("datasrc")) {
        element_database = elements["datasrc"];
    }

    XML_Node* doc = get_XML_File(element_database);
    XML_Node* dbe = &doc->child("ctml/elementData");

    XML_Node& root = phase.root();
    XML_Node* local_db = 0;
    if (root.hasChild("ctml")) {
        if (root.child("ctml").hasChild("elementData")) {
            local_db = &root.child("ctml/elementData");
        }
    }

    int nel = static_cast<int>(enames.size());
    int i;
    string enm;
    XML_Node* e = 0;
    for (i = 0; i < nel; i++) {
        e = 0;
        if (local_db) {
            //writelog("looking in local database.");
            e = local_db->findByAttr("name",enames[i]);
            //if (!e) writelog(enames[i]+" not found.");
        }
        if (!e) {
            e = dbe->findByAttr("name",enames[i]);
        }
        if (e) {
            addUniqueElement(*e);
        } else {
            throw CanteraError("addElementsFromXML","no data for element "
                               +enames[i]);
        }
    }
}

void Phase::freezeElements()
{
    m_elementsFrozen = true;
}

bool Phase::elementsFrozen()
{
    return m_elementsFrozen;
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
    m_elementsFrozen = false;
    addUniqueElement(symbol, weight, atomicNumber, entropy298, elem_type);
    m_elementsFrozen = true;
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

void Phase::addSpecies(const std::string& name_, const doublereal* comp,
                       doublereal charge_, doublereal size_)
{
    freezeElements();
    m_speciesNames.push_back(name_);
    m_speciesCharge.push_back(charge_);
    m_speciesSize.push_back(size_);
    size_t ne = nElements();
    // Create a changeable copy of the element composition. We now change
    // the charge potentially
    vector_fp compNew(ne);
    for (size_t m = 0; m < ne; m++) {
        compNew[m] = comp[m];
    }
    double wt = 0.0;
    const vector_fp& aw = atomicWeights();
    if (charge_ != 0.0) {
        size_t eindex = elementIndex("E");
        if (eindex != npos) {
            doublereal ecomp = compNew[eindex];
            if (fabs(charge_ + ecomp) > 0.001) {
                if (ecomp != 0.0) {
                    throw CanteraError("Phase::addSpecies",
                                       "Input charge and element E compositions differ "
                                       "for species " + name_);
                } else {
                    // Just fix up the element E composition based on the input
                    // species charge
                    compNew[eindex] = -charge_;
                }
            }
        } else {
            addUniqueElementAfterFreeze("E", 0.000545, 0, 0.0,
                                        CT_ELEM_TYPE_ELECTRONCHARGE);
            ne = nElements();
            eindex = elementIndex("E");
            compNew.resize(ne);
            compNew[ne - 1] = - charge_;
        }
    }
    for (size_t m = 0; m < ne; m++) {
        m_speciesComp.push_back(compNew[m]);
        wt += compNew[m] * aw[m];
    }
    m_molwts.push_back(wt);
    m_kk++;
}

void Phase::addUniqueSpecies(const std::string& name_, const doublereal* comp,
                             doublereal charge_, doublereal size_)
{
    for (size_t k = 0; k < m_kk; k++) {
        if (m_speciesNames[k] == name_) {
            // We have found a match. Do some compatibility checks.
            for (size_t i = 0; i < m_mm; i++) {
                if (comp[i] != m_speciesComp[k * m_mm + i]) {
                    throw CanteraError("addUniqueSpecies",
                                       "Duplicate species have different "
                                       "compositions: " + name_);
                }
            }
            if (charge_ != m_speciesCharge[k]) {
                throw CanteraError("addUniqueSpecies",
                                   "Duplicate species have different "
                                   "charges: " + name_);
            }
            if (size_ != m_speciesSize[k]) {
                throw CanteraError("addUniqueSpecies",
                                   "Duplicate species have different "
                                   "sizes: " + name_);
            }
            return;
        }
    }
    addSpecies(name_, comp, charge_, size_);
}

void Phase::freezeSpecies()
{
    m_speciesFrozen = true;
    init(molecularWeights());
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
        // that have zero molecular weight. Give them a very small molecular
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
    return (m_kk > 0 && m_elementsFrozen && m_speciesFrozen);
}

} // namespace Cantera
