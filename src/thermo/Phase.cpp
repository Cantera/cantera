/**
 *  @file Phase.cpp
 *   Definition file for class Phase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Phase.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace std;

namespace Cantera
{

Phase::Phase() :
    m_kk(0),
    m_ndim(3),
    m_undefinedElementBehavior(UndefElement::add),
    m_caseSensitiveSpecies(false),
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
}

std::string Phase::name() const
{
    return m_name;
}

void Phase::setName(const std::string& name)
{
    m_name = name;
}

size_t Phase::nElements() const
{
    return m_mm;
}

void Phase::checkElementIndex(size_t m) const
{
    if (m >= m_mm) {
        throw IndexError("Phase::checkElementIndex", "elements", m, m_mm-1);
    }
}

void Phase::checkElementArraySize(size_t mm) const
{
    if (m_mm > mm) {
        throw ArraySizeError("Phase::checkElementArraySize", mm, m_mm);
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

size_t Phase::findSpeciesLower(const std::string& name) const
{
    std::string nLower = toLowerCopy(name);
    size_t loc = npos;
    auto it = m_speciesLower.find(nLower);
    if (it == m_speciesLower.end()) {
        return npos;
    } else {
        loc = it->second;
        if (loc == npos) {
            throw CanteraError("Phase::findSpeciesLower",
                               "Lowercase species name '{}' is not unique. "
                               "Set Phase::caseSensitiveSpecies to true to "
                               "enforce case sensitive species names", nLower);
        }
    }
    return loc;
}

size_t Phase::speciesIndex(const std::string& nameStr) const
{
    size_t loc = npos;

    auto it = m_speciesIndices.find(nameStr);
    if (it != m_speciesIndices.end()) {
        return it->second;
    } else if (!m_caseSensitiveSpecies) {
        loc = findSpeciesLower(nameStr);
    }
    return loc;
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
        throw IndexError("Phase::checkSpeciesIndex", "species", k, m_kk-1);
    }
}

void Phase::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("Phase::checkSpeciesArraySize", kk, m_kk);
    }
}

std::string Phase::speciesSPName(int k) const
{
    return m_name + ":" + speciesName(k);
}

std::map<std::string, size_t> Phase::nativeState() const
{
    if (isPure()) {
        if (isCompressible()) {
            return { {"T", 0}, {"D", 1} };
        } else {
            return { {"T", 0}, {"P", 1} };
        }
    } else {
        if (isCompressible()) {
            return { {"T", 0}, {"D", 1}, {"Y", 2} };
        } else {
            return { {"T", 0}, {"P", 1}, {"Y", 2} };
        }
    }
}

vector<std::string> Phase::fullStates() const
{
    if (isPure()) {
        if (isCompressible()) {
            return {"TD", "TP", "UV", "DP", "HP", "SP", "SV"};
        } else {
            return {"TP", "HP", "SP"};
        }
    } else {
        if (isCompressible()) {
            return {"TDX", "TDY", "TPX", "TPY", "UVX", "UVY", "DPX", "DPY",
                    "HPX", "HPY", "SPX", "SPY", "SVX", "SVY"};
        } else {
            return {"TPX", "TPY", "HPX", "HPY", "SPX", "SPY"};
        }
    }
}

vector<std::string> Phase::partialStates() const
{
    if (isPure()) {
        return {};
    } else {
        if (isCompressible()) {
            return {"TD", "TP", "UV", "DP", "HP", "SP", "SV"};
        } else {
            return {"TP", "HP", "SP"};
        }
    }
}

size_t Phase::stateSize() const {
    if (isPure()) {
        return 2;
    } else {
        return nSpecies() + 2;
    }
}

void Phase::saveState(vector_fp& state) const
{
    state.resize(stateSize());
    saveState(state.size(), &state[0]);
}

void Phase::saveState(size_t lenstate, doublereal* state) const
{
    auto native = nativeState();

    // function assumes default definition of nativeState
    state[native.at("T")] = temperature();
    if (isCompressible()) {
        state[native.at("D")] = density();
    } else {
        state[native.at("P")] = pressure();
    }
    if (native.count("X")) {
        getMoleFractions(state + native["X"]);
    } else if (native.count("Y")) {
        getMassFractions(state + native["Y"]);
    }
}

void Phase::restoreState(const vector_fp& state)
{
    restoreState(state.size(),&state[0]);
}

void Phase::restoreState(size_t lenstate, const double* state)
{
    size_t ls = stateSize();
    if (lenstate < ls) {
        throw ArraySizeError("Phase::restoreState",
                             lenstate, ls);
    }

    auto native = nativeState();
    setTemperature(state[native.at("T")]);
    if (isCompressible()) {
        setDensity(state[native.at("D")]);
    } else {
        setPressure(state[native.at("P")]);
    }

    if (native.count("X")) {
        setMoleFractions_NoNorm(state + native["X"]);
    } else if (native.count("Y")) {
        setMassFractions_NoNorm(state + native["Y"]);
    }
    compositionChanged();
}

void Phase::setMoleFractions(const double* const x)
{
    // Use m_y as a temporary work vector for the non-negative mole fractions
    double norm = 0.0;
    // sum is calculated below as the unnormalized molecular weight
    double sum = 0;
    for (size_t k = 0; k < m_kk; k++) {
        double xk = std::max(x[k], 0.0); // Ignore negative mole fractions
        m_y[k] = xk;
        norm += xk;
        sum += m_molwts[k] * xk;
    }

    // Set m_ym to the normalized mole fractions divided by the normalized mean
    // molecular weight:
    //     m_ym_k = X_k / (sum_k X_k M_k)
    const double invSum = 1.0/sum;
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

void Phase::setMoleFractions_NoNorm(const double* const x)
{
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    scale(x, x + m_kk, m_ym.begin(), 1.0/m_mmw);
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(),
              m_y.begin(), multiplies<double>());
    compositionChanged();
}

void Phase::setMoleFractionsByName(const compositionMap& xMap)
{
    vector_fp mf = getCompositionFromMap(xMap);
    setMoleFractions(mf.data());
}

void Phase::setMoleFractionsByName(const std::string& x)
{
    setMoleFractionsByName(parseCompString(x));
}

void Phase::setMassFractions(const double* const y)
{
    for (size_t k = 0; k < m_kk; k++) {
        m_y[k] = std::max(y[k], 0.0); // Ignore negative mass fractions
    }
    double norm = accumulate(m_y.begin(), m_y.end(), 0.0);
    scale(m_y.begin(), m_y.end(), m_y.begin(), 1.0/norm);

    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(),
              m_ym.begin(), multiplies<double>());
    m_mmw = 1.0 / accumulate(m_ym.begin(), m_ym.end(), 0.0);
    compositionChanged();
}

void Phase::setMassFractions_NoNorm(const double* const y)
{
    double sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
              multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    m_mmw = 1.0/sum;
    compositionChanged();
}

void Phase::setMassFractionsByName(const compositionMap& yMap)
{
    vector_fp mf = getCompositionFromMap(yMap);
    setMassFractions(mf.data());
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

const vector_fp& Phase::inverseMolecularWeights() const
{
    return m_rmolwts;
}

void Phase::getCharges(double* charges) const
{
    copy(m_speciesCharge.begin(), m_speciesCharge.end(), charges);
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

void Phase::getMoleFractions(double* const x) const
{
    scale(m_ym.begin(), m_ym.end(), x, m_mmw);
}

double Phase::moleFraction(size_t k) const
{
    checkSpeciesIndex(k);
    return m_ym[k] * m_mmw;
}

double Phase::moleFraction(const std::string& nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return moleFraction(iloc);
    } else {
        return 0.0;
    }
}

const double* Phase::moleFractdivMMW() const
{
    return &m_ym[0];
}

double Phase::massFraction(size_t k) const
{
    checkSpeciesIndex(k);
    return m_y[k];
}

double Phase::massFraction(const std::string& nameSpec) const
{
    size_t iloc = speciesIndex(nameSpec);
    if (iloc != npos) {
        return massFractions()[iloc];
    } else {
        return 0.0;
    }
}

void Phase::getMassFractions(double* const y) const
{
    copy(m_y.begin(), m_y.end(), y);
}

double Phase::concentration(const size_t k) const
{
    checkSpeciesIndex(k);
    return m_y[k] * m_dens * m_rmolwts[k];
}

void Phase::getConcentrations(double* const c) const
{
    scale(m_ym.begin(), m_ym.end(), c, m_dens);
}

void Phase::setConcentrations(const double* const conc)
{
    assertCompressible("setConcentrations");

    // Use m_y as temporary storage for non-negative concentrations
    double sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        double ck = std::max(conc[k], 0.0); // Ignore negative concentrations
        m_y[k] = ck;
        sum += ck * m_molwts[k];
        norm += ck;
    }
    m_mmw = sum/norm;
    setDensity(sum);
    double rsum = 1.0/sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = m_y[k] * rsum;
        m_y[k] = m_ym[k] * m_molwts[k]; // m_y is now the mass fraction
    }
    compositionChanged();
}

void Phase::setConcentrationsNoNorm(const double* const conc)
{
    assertCompressible("setConcentrationsNoNorm");

    double sum = 0.0, norm = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        sum += conc[k] * m_molwts[k];
        norm += conc[k];
    }
    m_mmw = sum/norm;
    setDensity(sum);
    double rsum = 1.0/sum;
    for (size_t k = 0; k != m_kk; ++k) {
        m_ym[k] = conc[k] * rsum;
        m_y[k] = m_ym[k] * m_molwts[k];
    }
    compositionChanged();
}

void Phase::setMolesNoTruncate(const double* const N)
{
    // get total moles
    copy(N, N + m_kk, m_ym.begin());
    double totalMoles = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    // get total mass
    copy(N, N + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_molwts.begin(), m_y.begin(), multiplies<double>());
    double totalMass = accumulate(m_y.begin(), m_y.end(), 0.0);
    // mean molecular weight
    m_mmw = totalMass/totalMoles;
    // mass fractions
    scale(m_y.begin(), m_y.end(), m_y.begin(), 1/totalMass);
    // moles fractions/m_mmw
    scale(m_ym.begin(), m_ym.end(), m_ym.begin(), 1/(m_mmw * totalMoles));
    // composition has changed
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

double Phase::molarDensity() const
{
    return density()/meanMolecularWeight();
}

void Phase::setMolarDensity(const double molar_density)
{
    assertCompressible("setMolarDensity");
    m_dens = molar_density*meanMolecularWeight();
}

double Phase::molarVolume() const
{
    return 1.0/molarDensity();
}

void Phase::setDensity(const double density_)
{
    assertCompressible("setDensity");
    if (density_ > 0.0) {
        m_dens = density_;
    } else {
        throw CanteraError("Phase::setDensity",
            "density must be positive. density = {}", density_);
    }
}

void Phase::assignDensity(const double density_)
{
    if (density_ > 0.0) {
        m_dens = density_;
    } else {
        throw CanteraError("Phase::assignDensity",
            "density must be positive. density = {}", density_);
    }
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
    double sumxlogx = 0;
    for (size_t k = 0; k < m_kk; k++) {
        sumxlogx += m_ym[k] * std::log(std::max(m_ym[k], SmallNumber));
    }
    return m_mmw * sumxlogx + std::log(m_mmw);
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
            const static AnyMap db = AnyMap::fromYamlFile(
                "element-standard-entropies.yaml");
            const AnyMap& elem = db["elements"].getMapWhere("symbol", symbol);
            entropy298 = elem.convert("entropy298", "J/kmol/K", ENTROPY298_UNKNOWN);
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
    if (m_species.find(spec->name) != m_species.end()) {
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

    size_t ne = nElements();
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

    double wt = 0.0;
    for (size_t m = 0; m < ne; m++) {
        wt += comp[m] * aw[m];
    }

    // Some surface phases may define species representing empty sites
    // that have zero molecular weight. Give them a very small molecular
    // weight to avoid dividing by zero.
    wt = std::max(wt, Tiny);

    spec->setMolecularWeight(wt);

    m_molwts.push_back(wt);
    m_rmolwts.push_back(1.0/wt);

    for (size_t m = 0; m < ne; m++) {
        m_speciesComp.push_back(comp[m]);
    }

    m_speciesNames.push_back(spec->name);
    m_species[spec->name] = spec;
    m_speciesIndices[spec->name] = m_kk;
    m_speciesCharge.push_back(spec->charge);

    std::string nLower = toLowerCopy(spec->name);
    if (m_speciesLower.find(nLower) == m_speciesLower.end()) {
        m_speciesLower[nLower] = m_kk;
    } else {
        m_speciesLower[nLower] = npos;
    }
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
    const shared_ptr<Species>& old = m_species[spec->name];
    if (spec->composition != old->composition) {
        throw CanteraError("Phase::modifySpecies",
            "New composition for '{}' does not match existing composition",
            spec->name);
    }
    m_species[spec->name] = spec;
    invalidateCache();
}

void Phase::addSpeciesAlias(const std::string& name, const std::string& alias)
{
    if (speciesIndex(alias) != npos) {
        throw CanteraError("Phase::addSpeciesAlias",
            "Invalid alias '{}': species already exists", alias);
    }
    size_t k = speciesIndex(name);
    if (k != npos) {
        m_speciesIndices[alias] = k;
    } else {
        throw CanteraError("Phase::addSpeciesAlias",
            "Unable to add alias '{}' "
            "(original species '{}' not found).", alias, name);
    }
}

vector<std::string> Phase::findIsomers(const compositionMap& compMap) const
{
    vector<std::string> isomerNames;

    for (const auto& k : m_species) {
        if (k.second->composition == compMap) {
            isomerNames.emplace_back(k.first);
        }
    }

    return isomerNames;
}

vector<std::string> Phase::findIsomers(const std::string& comp) const
{
    return findIsomers(parseCompString(comp));
}

shared_ptr<Species> Phase::species(const std::string& name) const
{
    size_t k = speciesIndex(name);
    if (k != npos) {
        return m_species.at(speciesName(k));
    } else {
        throw CanteraError("Phase::species",
                           "Unknown species '{}'", name);
    }
}

shared_ptr<Species> Phase::species(size_t k) const
{
    checkSpeciesIndex(k);
    return m_species.at(m_speciesNames[k]);
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
    m_stateNum++;
    m_cache.clear();
}

void Phase::setMolecularWeight(const int k, const double mw)
{
    m_molwts[k] = mw;
    m_rmolwts[k] = 1.0/mw;

    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
              multiplies<double>());
    m_mmw = 1.0 / accumulate(m_ym.begin(), m_ym.end(), 0.0);
}

void Phase::compositionChanged() {
    m_stateNum++;
}

vector_fp Phase::getCompositionFromMap(const compositionMap& comp) const
{
    vector_fp X(m_kk);
    for (const auto& sp : comp) {
        size_t loc = speciesIndex(sp.first);
        if (loc == npos) {
            throw CanteraError("Phase::getCompositionFromMap",
                               "Unknown species '{}'", sp.first);
        }
        X[loc] = sp.second;
    }
    return X;
}

void Phase::massFractionsToMoleFractions(const double* Y, double* X) const
{
    double rmmw = 0.0;
    for (size_t k = 0; k != m_kk; ++k) {
        rmmw += Y[k]/m_molwts[k];
    }
    if (rmmw == 0.0) {
        throw CanteraError("Phase::massFractionsToMoleFractions",
                           "no input composition given");
    }
    for (size_t k = 0; k != m_kk; ++k) {
        X[k] = Y[k]/(rmmw*m_molwts[k]);
    }
}

void Phase::moleFractionsToMassFractions(const double* X, double* Y) const
{
    double mmw = dot(X, X+m_kk, m_molwts.data());
    if (mmw == 0.0) {
        throw CanteraError("Phase::moleFractionsToMassFractions",
                           "no input composition given");
    }
    double rmmw = 1.0/mmw;
    for (size_t k = 0; k != m_kk; ++k) {
        Y[k] = X[k]*m_molwts[k]*rmmw;
    }
}

} // namespace Cantera
