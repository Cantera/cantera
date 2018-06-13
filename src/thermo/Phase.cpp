/**
 *  @file Phase.cpp
 *   Definition file for class Phase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Phase.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace std;

namespace Cantera
{

Phase::Phase() :
    m_kk(0),
    m_ndim(3),
    m_undefinedElementBehavior(UndefElement::add),
    m_xml(new XML_Node("phase")),
    m_id("<phase>"),
    m_temp(0.001),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1),
    m_mm(0),
    m_elem_type(0)
{
}

Phase::~Phase()
{
    if (m_xml) {
        XML_Node* xroot = &m_xml->root();
        delete xroot;
    }
    m_xml = 0;
}

XML_Node& Phase::xml() const
{
    return *m_xml;
}

void Phase::setXMLdata(XML_Node& xmlPhase)
{
    XML_Node* xroot = &xmlPhase.root();
    XML_Node *root_xml = new XML_Node();
    xroot->copy(root_xml);
    if (m_xml) {
       XML_Node *rOld = &m_xml->root();
       delete rOld;
       m_xml = 0;
    }
    m_xml = findXMLPhase(root_xml, xmlPhase.id());
    if (!m_xml) {
        throw CanteraError("Phase::setXMLdata()", "XML 'phase' node not found");
    }
    if (&m_xml->root() != root_xml) {
        throw CanteraError("Phase::setXMLdata()", "Root XML node not found");
    }
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
    checkElementIndex(m);
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
    size_t loc = getValue(m_speciesIndices, toLowerCopy(nameStr), npos);
    if (loc == npos && nameStr.find(':') != npos) {
        std::string pn;
        std::string sn = toLowerCopy(parseSpeciesName(nameStr, pn));
        if (pn == "" || pn == m_name || pn == m_id) {
            return getValue(m_speciesIndices, sn, npos);
        } else {
            return npos;
        }
    } else {
        return loc;
    }
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
    return m_name + ":" + speciesName(k);
}

void Phase::saveState(vector_fp& state) const
{
    state.resize(nSpecies() + 2);
    saveState(state.size(), &state[0]);
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
    compositionChanged();
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
    // sum is calculated below as the unnormalized molecular weight
    doublereal sum = 0;
    for (size_t k = 0; k < m_kk; k++) {
        double xk = std::max(x[k], 0.0); // Ignore negative mole fractions
        m_y[k] = xk;
        norm += xk;
        sum += m_molwts[k] * xk;
    }

    // Set m_ym_ to the normalized mole fractions divided by the normalized mean
    // molecular weight:
    //     m_ym_k = X_k / (sum_k X_k M_k)
    const doublereal invSum = 1.0/sum;
    for (size_t k=0; k < m_kk; k++) {
        m_ym[k] = m_y[k]*invSum;
    }

    // Now set m_y to the normalized mass fractions:
    //     m_y =  X_k M_k / (sum_k X_k M_k)
    for (size_t k=0; k < m_kk; k++) {
        m_y[k] = m_ym[k] * m_molwts[k];
    }

    // Calculate the normalized molecular weight
    m_mmw = sum/norm;
    compositionChanged();
}

void Phase::setMoleFractions_NoNorm(const doublereal* const x)
{
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(1.0/m_mmw));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(),
              m_y.begin(), multiplies<double>());
    compositionChanged();
}

void Phase::setMoleFractionsByName(const compositionMap& xMap)
{
    vector_fp mf(m_kk, 0.0);
    for (const auto& sp : xMap) {
        try {
            mf[m_speciesIndices.at(toLowerCopy(sp.first))] = sp.second;
        } catch (std::out_of_range&) {
            throw CanteraError("Phase::setMoleFractionsByName",
                               "Unknown species '{}'", sp.first);
        }
    }
    setMoleFractions(&mf[0]);
}

void Phase::setMoleFractionsByName(const std::string& x)
{
    setMoleFractionsByName(parseCompString(x));
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
    compositionChanged();
}

void Phase::setMassFractions_NoNorm(const doublereal* const y)
{
    doublereal sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
              multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    m_mmw = 1.0/sum;
    compositionChanged();
}

void Phase::setMassFractionsByName(const compositionMap& yMap)
{
    vector_fp mf(m_kk, 0.0);
    for (const auto& sp : yMap) {
        try {
            mf[m_speciesIndices.at(toLowerCopy(sp.first))] = sp.second;
        } catch (std::out_of_range&) {
            throw CanteraError("Phase::setMassFractionsByName",
                               "Unknown species '{}'", sp.first);
        }
    }
    setMassFractions(&mf[0]);
}

void Phase::setMassFractionsByName(const std::string& y)
{
    setMassFractionsByName(parseCompString(y));
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

void Phase::setState_TRX(doublereal t, doublereal dens, const compositionMap& x)
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

void Phase::setState_TRY(doublereal t, doublereal dens, const compositionMap& y)
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
    weights = molecularWeights();
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

compositionMap Phase::getMoleFractionsByName(double threshold) const
{
    compositionMap comp;
    for (size_t k = 0; k < m_kk; k++) {
        double x = moleFraction(k);
        if (x > threshold) {
            comp[speciesName(k)] = x;
        }
    }
    return comp;
}

compositionMap Phase::getMassFractionsByName(double threshold) const
{
    compositionMap comp;
    for (size_t k = 0; k < m_kk; k++) {
        double x = massFraction(k);
        if (x > threshold) {
            comp[speciesName(k)] = x;
        }
    }
    return comp;
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
    return m_y[k] * m_dens * m_rmolwts[k];
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
    compositionChanged();
}

void Phase::setConcentrationsNoNorm(const double* const conc)
{
    doublereal sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        sum += conc[k] * m_molwts[k];
        norm += conc[k];
    }
    m_mmw = sum/norm;
    setDensity(sum);
    doublereal rsum = 1.0/sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = conc[k] * rsum;
        m_y[k] = m_ym[k] * m_molwts[k];
    }
    compositionChanged();
}

doublereal Phase::elementalMassFraction(const size_t m) const
{
    checkElementIndex(m);
    doublereal Z_m = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        Z_m += nAtoms(k, m) * atomicWeight(m) / molecularWeight(k)
            * massFraction(k);
    }
    return Z_m;
}

doublereal Phase::elementalMoleFraction(const size_t m) const
{
    checkElementIndex(m);
    double denom = 0;
    for (size_t k = 0; k < m_kk; k++) {
        double atoms = 0;
        for (size_t j = 0; j < nElements(); j++) {
            atoms += nAtoms(k, j);
        }
        denom += atoms * moleFraction(k);
    }
    doublereal numerator = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        numerator += nAtoms(k, m) * moleFraction(k);
    }
    return numerator / denom;
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
    doublereal cdens = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        cdens += charge(k)*moleFraction(k);
    }
    return cdens * Faraday;
}

doublereal Phase::mean_X(const doublereal* const Q) const
{
    return m_mmw*std::inner_product(m_ym.begin(), m_ym.end(), Q, 0.0);
}

doublereal Phase::mean_X(const vector_fp& Q) const
{
    return m_mmw*std::inner_product(m_ym.begin(), m_ym.end(), Q.begin(), 0.0);
}

doublereal Phase::sum_xlogx() const
{
    return m_mmw* Cantera::sum_xlogx(m_ym.begin(), m_ym.end()) + log(m_mmw);
}

size_t Phase::addElement(const std::string& symbol, doublereal weight,
                         int atomic_number, doublereal entropy298,
                         int elem_type)
{
    // Look up the atomic weight if not given
    if (weight == 0.0) {
        try {
            weight = getElementWeight(symbol);
        } catch (CanteraError&) {
            // assume this is just a custom element with zero atomic weight
        }
    } else if (weight == -12345.0) {
        weight = getElementWeight(symbol);
    }

    // Try to look up the standard entropy if not given. Fail silently.
    if (entropy298 == ENTROPY298_UNKNOWN) {
        try {
            XML_Node* db = get_XML_File("elements.xml");
            XML_Node* elnode = db->findByAttr("name", symbol);
            if (elnode && elnode->hasChild("entropy298")) {
                entropy298 = fpValueCheck(elnode->child("entropy298")["value"]);
            }
        } catch (CanteraError&) {
        }
    }

    // Check for duplicates
    auto iter = find(m_elementNames.begin(), m_elementNames.end(), symbol);
    if (iter != m_elementNames.end()) {
        size_t m = iter - m_elementNames.begin();
        if (m_atomicWeights[m] != weight) {
            throw CanteraError("Phase::addElement",
                "Duplicate elements ({}) have different weights", symbol);
        } else {
            // Ignore attempt to add duplicate element with the same weight
            return m;
        }
    }

    // Add the new element
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

    // Update species compositions
    if (m_kk) {
        vector_fp old(m_speciesComp);
        m_speciesComp.resize(m_kk*m_mm, 0.0);
        for (size_t k = 0; k < m_kk; k++) {
            size_t m_old = m_mm - 1;
            for (size_t m = 0; m < m_old; m++) {
                m_speciesComp[k * m_mm + m] = old[k * (m_old) + m];
            }
            m_speciesComp[k * (m_mm) + (m_mm-1)] = 0.0;
        }
    }

    return m_mm-1;
}

bool Phase::addSpecies(shared_ptr<Species> spec) {
    if (m_species.find(toLowerCopy(spec->name)) != m_species.end()) {
        throw CanteraError("Phase::addSpecies",
            "Phase '{}' already contains a species named '{}'.",
            m_name, spec->name);
    }
    vector_fp comp(nElements());
    for (const auto& elem : spec->composition) {
        size_t m = elementIndex(elem.first);
        if (m == npos) { // Element doesn't exist in this phase
            switch (m_undefinedElementBehavior) {
            case UndefElement::ignore:
                return false;

            case UndefElement::add:
                addElement(elem.first);
                comp.resize(nElements());
                m = elementIndex(elem.first);
                break;

            case UndefElement::error:
            default:
                throw CanteraError("Phase::addSpecies",
                    "Species '{}' contains an undefined element '{}'.",
                    spec->name, elem.first);
            }
        }
        comp[m] = elem.second;
    }

    m_speciesNames.push_back(spec->name);
    m_species[toLowerCopy(spec->name)] = spec;
    m_speciesIndices[toLowerCopy(spec->name)] = m_kk;
    m_speciesCharge.push_back(spec->charge);
    size_t ne = nElements();

    double wt = 0.0;
    const vector_fp& aw = atomicWeights();
    if (spec->charge != 0.0) {
        size_t eindex = elementIndex("E");
        if (eindex != npos) {
            doublereal ecomp = comp[eindex];
            if (fabs(spec->charge + ecomp) > 0.001) {
                if (ecomp != 0.0) {
                    throw CanteraError("Phase::addSpecies",
                                       "Input charge and element E compositions differ "
                                       "for species " + spec->name);
                } else {
                    // Just fix up the element E composition based on the input
                    // species charge
                    comp[eindex] = -spec->charge;
                }
            }
        } else {
            addElement("E", 0.000545, 0, 0.0, CT_ELEM_TYPE_ELECTRONCHARGE);
            ne = nElements();
            eindex = elementIndex("E");
            comp.resize(ne);
            comp[ne - 1] = - spec->charge;
        }
    }
    for (size_t m = 0; m < ne; m++) {
        m_speciesComp.push_back(comp[m]);
        wt += comp[m] * aw[m];
    }

    // Some surface phases may define species representing empty sites
    // that have zero molecular weight. Give them a very small molecular
    // weight to avoid dividing by zero.
    wt = std::max(wt, Tiny);
    m_molwts.push_back(wt);
    m_rmolwts.push_back(1.0/wt);
    m_kk++;

    // Ensure that the Phase has a valid mass fraction vector that sums to
    // one. We will assume that species 0 has a mass fraction of 1.0 and mass
    // fraction of all other species is 0.0.
    if (m_kk == 1) {
        m_y.push_back(1.0);
        m_ym.push_back(m_rmolwts[0]);
        m_mmw = 1.0 / m_ym[0];
    } else {
        m_y.push_back(0.0);
        m_ym.push_back(0.0);
    }
    invalidateCache();
    return true;
}

void Phase::modifySpecies(size_t k, shared_ptr<Species> spec)
{
    if (speciesName(k) != spec->name) {
        throw CanteraError("Phase::modifySpecies",
            "New species name '{}' does not match existing name '{}'",
                           spec->name, speciesName(k));
    }
    const shared_ptr<Species>& old = m_species[toLowerCopy(spec->name)];
    if (spec->composition != old->composition) {
        throw CanteraError("Phase::modifySpecies",
            "New composition for '{}' does not match existing composition",
            spec->name);
    }
    m_species[toLowerCopy(spec->name)] = spec;
    invalidateCache();
}

shared_ptr<Species> Phase::species(const std::string& name) const
{
    return m_species.at(toLowerCopy(name));
}

shared_ptr<Species> Phase::species(size_t k) const
{
    return species(m_speciesNames[k]);
}

void Phase::ignoreUndefinedElements() {
    m_undefinedElementBehavior = UndefElement::ignore;
}

void Phase::addUndefinedElements() {
    m_undefinedElementBehavior = UndefElement::add;
}

void Phase::throwUndefinedElements() {
    m_undefinedElementBehavior = UndefElement::error;
}

bool Phase::ready() const
{
    return (m_kk > 0);
}

void Phase::invalidateCache() {
    m_cache.clear();
}

void Phase::compositionChanged() {
    m_stateNum++;
}

} // namespace Cantera
