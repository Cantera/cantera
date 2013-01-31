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

namespace Cantera {
template<typename ValAndDerivType>
Phase<ValAndDerivType>::Phase() :
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

template<typename ValAndDerivType>
Phase<ValAndDerivType>::Phase(const Phase& right) :
        m_kk(0),
        m_ndim(3),
        m_xml(0),
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
    // Use the assignment operator to do the actual copying
    *this = operator=(right);
}

template<typename ValAndDerivType>
Phase<ValAndDerivType>& Phase<ValAndDerivType>::operator=(const Phase<ValAndDerivType>& right)
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

    m_mm = right.m_mm;
    m_elementsFrozen = right.m_elementsFrozen;
    m_atomicWeights = right.m_atomicWeights;
    m_atomicNumbers = right.m_atomicNumbers;
    m_elementNames = right.m_elementNames;
    m_entropy298 = right.m_entropy298;
    m_elem_type = right.m_elem_type;
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
    m_id = right.m_id;
    m_name = right.m_name;

    return *this;
}

// Destructor.
template<typename ValAndDerivType>
Phase<ValAndDerivType>::~Phase()
{
    if (m_xml) {
        delete m_xml;
        m_xml = 0;
    }
}

template<typename ValAndDerivType>
XML_Node& Phase<ValAndDerivType>::xml()
{
    return *m_xml;
}

template<typename ValAndDerivType>
std::string Phase<ValAndDerivType>::id() const
{
    return m_id;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setID(const std::string& id)
{
    m_id = id;
}

template<typename ValAndDerivType>
std::string Phase<ValAndDerivType>::name() const
{
    return m_name;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setName(const std::string& nm)
{
    m_name = nm;
}

template<typename ValAndDerivType>
size_t Phase<ValAndDerivType>::nElements() const
{
    return m_mm;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::checkElementIndex(size_t m) const
{
    if (m >= m_mm) {
        throw IndexError("checkElementIndex", "elements", m, m_mm - 1);
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::checkElementArraySize(size_t mm) const
{
    if (m_mm > mm) {
        throw ArraySizeError("checkElementArraySize", mm, m_mm);
    }
}

template<typename ValAndDerivType>
string Phase<ValAndDerivType>::elementName(size_t m) const
{
    checkElementIndex(m);
    return m_elementNames[m];
}

template<typename ValAndDerivType>
size_t Phase<ValAndDerivType>::elementIndex(const std::string& name) const
{
    for (size_t i = 0; i < m_mm; i++) {
        if (m_elementNames[i] == name) {
            return i;
        }
    }
    return npos;
}

template<typename ValAndDerivType>
const vector<string>& Phase<ValAndDerivType>::elementNames() const
{
    return m_elementNames;
}

template<typename ValAndDerivType>
doublereal Phase<ValAndDerivType>::atomicWeight(size_t m) const
{
    return m_atomicWeights[m];
}

template<typename ValAndDerivType>
doublereal Phase<ValAndDerivType>::entropyElement298(size_t m) const
{
    AssertThrowMsg(m_entropy298[m] != ENTROPY298_UNKNOWN, "Elements::entropy298", "Entropy at 298 K of element is unknown");
    AssertTrace(m < m_mm);
    return (m_entropy298[m]);
}

template<typename ValAndDerivType>
const vector_fp& Phase<ValAndDerivType>::atomicWeights() const
{
    return m_atomicWeights;
}

template<typename ValAndDerivType>
int Phase<ValAndDerivType>::atomicNumber(size_t m) const
{
    return m_atomicNumbers[m];
}

template<typename ValAndDerivType>
int Phase<ValAndDerivType>::elementType(size_t m) const
{
    return m_elem_type[m];
}

template<typename ValAndDerivType>
int Phase<ValAndDerivType>::changeElementType(int m, int elem_type)
{
    int old = m_elem_type[m];
    m_elem_type[m] = elem_type;
    return old;
}

template<typename ValAndDerivType>
doublereal Phase<ValAndDerivType>::nAtoms(size_t k, size_t m) const
{
    checkElementIndex(m);
    checkSpeciesIndex(k);
    return m_speciesComp[m_mm * k + m];
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getAtoms(size_t k, double* atomArray) const
{
    for (size_t m = 0; m < m_mm; m++) {
        atomArray[m] = (double) m_speciesComp[m_mm * k + m];
    }
}

template<typename ValAndDerivType>
size_t Phase<ValAndDerivType>::speciesIndex(const std::string& nameStr) const
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

template<typename ValAndDerivType>
string Phase<ValAndDerivType>::speciesName(size_t k) const
{
    checkSpeciesIndex(k);
    return m_speciesNames[k];
}

template<typename ValAndDerivType>
const vector<string>& Phase<ValAndDerivType>::speciesNames() const
{
    return m_speciesNames;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::checkSpeciesIndex(size_t k) const
{
    if (k >= m_kk) {
        throw IndexError("checkSpeciesIndex", "species", k, m_kk - 1);
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_kk);
    }
}

template<typename ValAndDerivType>
std::string Phase<ValAndDerivType>::speciesSPName(int k) const
{
    std::string sn = speciesName(k);
    return (m_name + ":" + sn);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::saveState(vector_fp& state) const
{
    state.resize(nSpecies() + 2);
    saveState(state.size(), &(state[0]));
}
//=====================================================================================================================
template<typename ValAndDerivType>
void Phase<ValAndDerivType>::saveState(size_t lenstate, doublereal* state) const
{
    state[0] = temperature();
    state[1] = density();
    getMassFractions(state + 2);
}

template<>
void Phase<doubleFAD>::saveState(size_t lenstate, doublereal* state) const
{
    state[0] = temperature();
    state[1] = density().val();

    getMassFractions(DATA_PTR(m_vdtmps));
    for (size_t k = 0; k < m_kk; k++) {
        state[k + 2] = m_vdtmps[k].val();
    }

}

//=====================================================================================================================

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::restoreState(const vector_fp& state)
{
    restoreState(state.size(), &state[0]);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::restoreState(size_t lenstate, const doublereal* state)
{
    if (lenstate >= nSpecies() + 2) {
        setMassFractions_NoNorm(state + 2);
        setTemperature(state[0]);
        setDensity(state[1]);
    } else {
        throw ArraySizeError("Phase<ValAndDerivType>::restoreState", lenstate, nSpecies() + 2);
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMoleFractions(const doublereal* const x)
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
//    transform(m_y.begin(), m_y.end(), m_ym.begin(), timesConstant<double>(1.0/sum));
    const doublereal invSum = 1.0 / sum;
    for (size_t k = 0; k < m_kk; k++) {
        m_ym[k] = m_y[k] * invSum;
    }
    /*
     * Now set m_y to the normalized mass fractions
     *          m_y =  X_k M_k / (sum_k X_k M_k)
     */
//    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(), m_y.begin(), multiplies<double>());
    for (size_t k = 0; k < m_kk; k++) {
        m_y[k] = m_ym[k] * m_molwts[k];
    }
    /*
     * Calculate the normalized molecular weight
     */
    m_mmw = sum / norm;
    m_stateNum++;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMoleFractions_NoNorm(const doublereal* const x)
{
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    ValAndDerivType rmmw = 1.0 / m_mmw;
    transform(x, x + m_kk, m_ym.begin(), timesConstantVal<ValAndDerivType>(rmmw));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(), m_y.begin(), multiplies<ValAndDerivType>());
    m_stateNum++;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMoleFractionsByName(compositionMap& xMap)
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

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMoleFractionsByName(const std::string& x)
{
    compositionMap c = parseCompString(x, speciesNames());
    setMoleFractionsByName(c);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMassFractions(const doublereal* const y)
{
    for (size_t k = 0; k < m_kk; k++) {
        m_y[k] = std::max(y[k], 0.0); // Ignore negative mass fractions
    }
    //   doublereal norm = accumulate(m_y.begin(), m_y.end(), 0.0);
    ValAndDerivType norm = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        norm += m_y[k];
    }
    //  scale(m_y.begin(), m_y.end(), m_y.begin(), 1.0 / norm);
    for (size_t k = 0; k < m_kk; k++) {
        m_y[k] /= norm;
    }

    // norm2 = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(), multiplies<ValAndDerivType>());
    ValAndDerivType norm2 = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        norm2 += m_ym[k];
    }
    m_mmw = 1.0 / norm2;
    m_stateNum++;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMassFractions_NoNorm(const doublereal* const y)
{
    ValAndDerivType sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(), multiplies<ValAndDerivType>());

    //  sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    for (size_t k = 0; k < m_kk; k++) {
        sum += m_ym[k];
    }

    m_mmw = 1.0 / sum;
    m_stateNum++;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMassFractionsByName(compositionMap& yMap)
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

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMassFractionsByName(const std::string& y)
{
    compositionMap c = parseCompString(y, speciesNames());
    setMassFractionsByName(c);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TRX(doublereal t, doublereal dens, const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setDensity(dens);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TNX(doublereal t, doublereal n, const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setMolarDensity(n);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TRX(doublereal t, doublereal dens, compositionMap& x)
{
    setMoleFractionsByName(x);
    setTemperature(t);
    setDensity(dens);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TRY(doublereal t, doublereal dens, const doublereal* y)
{
    setMassFractions(y);
    setTemperature(t);
    setDensity(dens);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TRY(doublereal t, doublereal dens, compositionMap& y)
{
    setMassFractionsByName(y);
    setTemperature(t);
    setDensity(dens);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TR(doublereal t, doublereal rho)
{
    setTemperature(t);
    setDensity(rho);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TX(doublereal t, doublereal* x)
{
    setTemperature(t);
    setMoleFractions(x);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_TY(doublereal t, doublereal* y)
{
    setTemperature(t);
    setMassFractions(y);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_RX(doublereal rho, doublereal* x)
{
    setMoleFractions(x);
    setDensity(rho);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setState_RY(doublereal rho, doublereal* y)
{
    setMassFractions(y);
    setDensity(rho);
}

template<typename ValAndDerivType>
doublereal Phase<ValAndDerivType>::molecularWeight(size_t k) const
{
    checkSpeciesIndex(k);
    return m_molwts[k];
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMolecularWeights(vector_fp& weights) const
{
    const vector_fp& mw = molecularWeights();
    if (weights.size() < mw.size()) {
        weights.resize(mw.size());
    }
    copy(mw.begin(), mw.end(), weights.begin());
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMolecularWeights(doublereal* weights) const
{
    const vector_fp& mw = molecularWeights();
    copy(mw.begin(), mw.end(), weights);
}

template<typename ValAndDerivType>
const vector_fp& Phase<ValAndDerivType>::molecularWeights() const
{
    return m_molwts;
}

//=====================================================================================================================

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMoleFractionsByName(compositionMap& x) const
{
    x.clear();
    size_t kk = nSpecies();
    for (size_t k = 0; k < kk; k++) {
        x[speciesName(k)] = Phase<ValAndDerivType>::moleFraction(k);
    }
}

template<>
void Phase<doubleFAD>::getMoleFractionsByName(compositionMap& x) const
{
    x.clear();
    size_t kk = nSpecies();
    for (size_t k = 0; k < kk; k++) {
        doublereal xval = (Phase<doubleFAD>::moleFraction(k)).val();
        x[speciesName(k)] = xval;
    }
}

//=====================================================================================================================

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMoleFractions(ValAndDerivType* const x) const
{
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = m_mmw * m_ym[k];
    }
    //scale(m_ym.begin(), m_ym.end(), x, m_mmw);
}

// HKM This works
template<typename ValAndDerivType>
template<typename OutType>
void Phase<ValAndDerivType>::getMoleFractions(OutType * const x) const
{
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = m_mmw * m_ym[k];
    }
}

// Template specialization to get a vector of doubles out of the doubleFAD templated Phase object
/*
 * This is a little dangerous, because we have stripped out the derivative information just by putting the wrong argument on the command line.
 *
 */
template<>
template<>
void Phase<doubleFAD>::getMoleFractions(doublereal * const x) const
{
    doublereal mmw = m_mmw.val();
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = mmw * m_ym[k].val();
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMoleFractionsNoDeriv(doublereal* const x) const
{
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = m_mmw * m_ym[k];
    }
    //scale(m_ym.begin(), m_ym.end(), x, m_mmw);
}

template<>
void Phase<doubleFAD>::getMoleFractionsNoDeriv(doublereal* const x) const
{
    doublereal mmw = m_mmw.val();
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = mmw * m_ym[k].val();
    }
    //scale(m_ym.begin(), m_ym.end(), x, m_mmw);
}

//=====================================================================================================================
template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::moleFraction(size_t k) const
{
    checkSpeciesIndex(k);
    return m_ym[k] * m_mmw;
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::moleFraction(const std::string& nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return moleFraction(iloc);
    } else {
        return 0.0;
    }
}

template<typename ValAndDerivType>
const ValAndDerivType* Phase<ValAndDerivType>::moleFractdivMMW() const
{
    return &m_ym[0];
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::massFraction(size_t k) const
{
    checkSpeciesIndex(k);
    return m_y[k];
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::massFraction(const std::string& nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return massFractions()[iloc];
    } else {
        return 0.0;
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getMassFractions(ValAndDerivType* const y) const
{
    copy(m_y.begin(), m_y.end(), y);
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::concentration(const size_t k) const
{
    checkSpeciesIndex(k);
    return m_y[k] * m_dens * m_rmolwts[k];
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::getConcentrations(ValAndDerivType* const c) const
{
    // scale(m_ym.begin(), m_ym.end(), c, m_dens);
    for (size_t k = 0; k != m_kk; ++k) {
        c[k] = m_dens * m_ym[k];
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setConcentrations(const doublereal* const conc)
{
    // Use m_y as temporary storage for non-negative concentrations
    doublereal sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        double ck = std::max(conc[k], 0.0); // Ignore negative concentrations
        m_y[k] = ck;
        sum += ck * m_molwts[k];
        norm += ck;
    }
    m_mmw = sum / norm;
    setDensity(sum);
    doublereal rsum = 1.0 / sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = m_y[k] * rsum;
        m_y[k] = m_ym[k] * m_molwts[k]; // m_y is now the mass fraction
    }
    m_stateNum++;
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::molarDensity() const
{
    ValAndDerivType dd = density() / meanMolecularWeight();
    return dd;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::setMolarDensity(const doublereal molarDensity)
{
    m_dens = molarDensity * meanMolecularWeight();
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::molarVolume() const
{
    return 1.0 / molarDensity();
}

template<typename ValAndDerivType>
doublereal Phase<ValAndDerivType>::charge(size_t k) const
{
    return m_speciesCharge[k];
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::chargeDensity() const
{
    size_t kk = nSpecies();
    ValAndDerivType cdens = 0.0;
    for (size_t k = 0; k < kk; k++) {
        cdens += charge(k) * moleFraction(k);
    }
    cdens *= Faraday;
    return cdens;
}

// HKM The concept of using Xi/W as a basic storage needs to be reworked when going to a derivative implementation
template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::mean_X(const ValAndDerivType* const Q) const
{
    ValAndDerivType zz = 0.0;
    return m_mmw * std::inner_product(m_ym.begin(), m_ym.end(), Q, zz);
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::mean_Y(const doublereal* const Q) const
{
    // return dot(m_y.begin(), m_y.end(), Q);
    ValAndDerivType sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += m_y[k] * Q[k];
    }
    return sum;
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::sum_xlogx() const
{
    return m_mmw * Cantera::sum_xlogx_valderiv<ValAndDerivType>(m_ym.begin(), m_ym.end()) + log(m_mmw);
}

template<typename ValAndDerivType>
ValAndDerivType Phase<ValAndDerivType>::sum_xlogQ(doublereal* Q) const
{
    ValAndDerivType sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += (m_ym[k]) * std::log(Q[k] + Tiny);
    }
    return m_mmw * sum;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addElement(const std::string& symbol, doublereal weight)
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

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addElement(const XML_Node& e)
{
    doublereal weight = atof(e["atomicWt"].c_str());
    string symbol = e["name"];
    addElement(symbol, weight);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addUniqueElement(const std::string& symbol, doublereal weight, int atomicNumber, doublereal entropy298,
                                              int elem_type)
{
    if (weight == -12345.0) {
        weight = LookupWtElements(symbol);
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
    for (vector<string>::const_iterator it = m_elementNames.begin(); it < m_elementNames.end(); ++it, ++i) {
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
        m_atomicNumbers.push_back(atomicNumber);
        m_entropy298.push_back(entropy298);
        if (symbol == "E") {
            m_elem_type.push_back(CT_ELEM_TYPE_ELECTRONCHARGE);
        } else {
            m_elem_type.push_back(elem_type);
        }
        m_mm++;
    } else {
        if (m_atomicWeights[i] != weight) {
            throw CanteraError("AddUniqueElement", "Duplicate Elements (" + symbol + ") have different weights");
        }
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addUniqueElement(const XML_Node& e)
{
    doublereal weight = 0.0;
    if (e.hasAttrib("atomicWt")) {
        weight = atof(stripws(e["atomicWt"]).c_str());
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
            entropy298 = atofCheck(stripws(e298Node["value"]).c_str());
        }
    }
    if (weight != 0.0) {
        addUniqueElement(symbol, weight, anum, entropy298);
    } else {
        addUniqueElement(symbol);
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addElementsFromXML(const XML_Node& phase)
{
    // get the declared element names
    if (!phase.hasChild("elementArray")) {
        throw CanteraError("Elements::addElementsFromXML", "phase xml node doesn't have \"elementArray\" XML Node");
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
            e = local_db->findByAttr("name", enames[i]);
            //if (!e) writelog(enames[i]+" not found.");
        }
        if (!e) {
            e = dbe->findByAttr("name", enames[i]);
        }
        if (e) {
            addUniqueElement(*e);
        } else {
            throw CanteraError("addElementsFromXML", "no data for element " + enames[i]);
        }
    }
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::freezeElements()
{
    m_elementsFrozen = true;
}

template<typename ValAndDerivType>
bool Phase<ValAndDerivType>::elementsFrozen()
{
    return m_elementsFrozen;
}

template<typename ValAndDerivType>
size_t Phase<ValAndDerivType>::addUniqueElementAfterFreeze(const std::string& symbol, doublereal weight, int atomicNumber,
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
    if (ii != m_mm - 1) {
        throw CanteraError("Phase<ValAndDerivType>::addElementAfterFreeze()", "confused");
    }
    if (m_kk > 0) {
        vector_fp old(m_speciesComp);
        m_speciesComp.resize(m_kk * m_mm, 0.0);
        for (size_t k = 0; k < m_kk; k++) {
            size_t m_old = m_mm - 1;
            for (size_t m = 0; m < m_old; m++) {
                m_speciesComp[k * m_mm + m] = old[k * (m_old) + m];
            }
            m_speciesComp[k * (m_mm) + (m_mm - 1)] = 0.0;
        }
    }
    return ii;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addSpecies(const std::string& name, const doublereal* comp, doublereal charge, doublereal size)
{
    freezeElements();
    m_speciesNames.push_back(name);
    m_speciesCharge.push_back(charge);
    m_speciesSize.push_back(size);
    size_t ne = nElements();
    // Create a changeable copy of the element composition. We now change
    // the charge potentially
    vector_fp compNew(ne);
    for (size_t m = 0; m < ne; m++) {
        compNew[m] = comp[m];
    }
    double wt = 0.0;
    const vector_fp& aw = atomicWeights();
    if (charge != 0.0) {
        size_t eindex = elementIndex("E");
        if (eindex != npos) {
            doublereal ecomp = compNew[eindex];
            if (fabs(charge + ecomp) > 0.001) {
                if (ecomp != 0.0) {
                    throw CanteraError("Phase<ValAndDerivType>::addSpecies", "Input charge and element E compositions differ "
                            "for species " + name);
                } else {
                    // Just fix up the element E composition based on the input
                    // species charge
                    compNew[eindex] = -charge;
                }
            }
        } else {
            addUniqueElementAfterFreeze("E", 0.000545, 0, 0.0, CT_ELEM_TYPE_ELECTRONCHARGE);
            ne = nElements();
            eindex = elementIndex("E");
            compNew.resize(ne);
            compNew[ne - 1] = -charge;
        }
    }
    for (size_t m = 0; m < ne; m++) {
        m_speciesComp.push_back(compNew[m]);
        wt += compNew[m] * aw[m];
    }
    m_molwts.push_back(wt);
    m_kk++;
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::addUniqueSpecies(const std::string& name, const doublereal* comp, doublereal charge, doublereal size)
{
    for (size_t k = 0; k < m_kk; k++) {
        if (m_speciesNames[k] == name) {
            // We have found a match. Do some compatibility checks.
            for (size_t i = 0; i < m_mm; i++) {
                if (comp[i] != m_speciesComp[k * m_mm + i]) {
                    throw CanteraError("addUniqueSpecies", "Duplicate species have different "
                            "compositions: " + name);
                }
            }
            if (charge != m_speciesCharge[k]) {
                throw CanteraError("addUniqueSpecies", "Duplicate species have different "
                        "charges: " + name);
            }
            if (size != m_speciesSize[k]) {
                throw CanteraError("addUniqueSpecies", "Duplicate species have different "
                        "sizes: " + name);
            }
            return;
        }
    }
    addSpecies(name, comp, charge, size);
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::freezeSpecies()
{
    m_speciesFrozen = true;
    init(molecularWeights());
}

template<typename ValAndDerivType>
void Phase<ValAndDerivType>::init(const vector_fp& mw)
{
    m_kk = mw.size();
    m_rmolwts.resize(m_kk);
    m_y.resize(m_kk, 0.0);
    m_ym.resize(m_kk, 0.0);
    m_vdtmps.resize(m_kk, 0.0);
    copy(mw.begin(), mw.end(), m_molwts.begin());
    for (size_t k = 0; k < m_kk; k++) {
        if (m_molwts[k] < 0.0) {
            throw CanteraError("Phase<ValAndDerivType>::init", "negative molecular weight for species number " + int2str(k));
        }

        // Some surface phases may define species representing empty sites
        // that have zero molecular weight. Give them a very small molecular
        // weight to avoid dividing by zero.
        if (m_molwts[k] < Tiny) {
            m_molwts[k] = Tiny;
        }
        m_rmolwts[k] = 1.0 / m_molwts[k];
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

template<typename ValAndDerivType>
bool Phase<ValAndDerivType>::ready() const
{
    return (m_kk > 0 && m_elementsFrozen && m_speciesFrozen);
}

//! Explicit Instantiation
template class Phase<doublereal> ;

#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class Phase<doubleFAD> ;
#endif
#endif

} // namespace Cantera
