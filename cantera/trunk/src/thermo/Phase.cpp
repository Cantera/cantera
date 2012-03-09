/**
 *  @file Phase.cpp
 *   Definition file for class, Phase, which contains functions for setting the
 *   state of a phase, and for referencing species by name
 *   (see \ref phases and class \link Cantera::Phase Phase\endlink).
 */

// Copyright 2001  California Institute of Technology

#include "cantera/base/ct_defs.h"
#include "cantera/thermo/Phase.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/ctexceptions.h"

using namespace std;

namespace Cantera
{

Phase::Phase() :
    Constituents(),
    m_kk(0),
    m_ndim(3),
    m_xml(new XML_Node("phase")),
    m_id("<phase>"),
    m_name(""),
    m_temp(0.0),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1)
{
}

/*
 * Copy Constructor
 *
 * This function just does the default initialization, and
 * then calls the assignment operator.
 */
Phase::Phase(const Phase& right) :
    Constituents(),
    m_kk(0),
    m_ndim(3),
    m_xml(new XML_Node("phase")),
    m_id("<phase>"),
    m_name(""),
    m_temp(0.0),
    m_dens(0.001),
    m_mmw(0.0),
    m_stateNum(-1)
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}

/*
 * Assignment operator
 *
 * This operation is sort of complicated. We have to
 * call the assignment operator for the Constituents and
 * State operators that Phase inherits from. Then,
 * we have to copy our own data, making sure to do a
 * deep copy on the XML_Node data owned by this object.
 */
Phase& Phase::operator=(const Phase& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }
    /*
     * Now call the inherited-classes assignment operators.
     */
    (void) Constituents::operator=(right);
    /*
     * Handle its own data
     */
    m_kk    = right.m_kk;
    m_ndim  = right.m_ndim;
    m_data  = right.m_data;
    m_temp           = right.m_temp;
    m_dens           = right.m_dens;
    m_mmw            = right.m_mmw;
    m_ym             = right.m_ym;
    m_y              = right.m_y;
    m_molwts         = right.m_molwts;
    m_rmolwts        = right.m_rmolwts;
    m_stateNum = -1;

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
        m_xml   = new XML_Node();
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

// Returns the index of a species named 'name' within the Phase object
/*
 * The first species in the phase will have an index 0, and the last one in the
 * phase will have an index of nSpecies() - 1.
 *
 *
 *  A species name may be referred to via three methods:
 *
 *    -   "speciesName"
 *    -   "PhaseId:speciesName"
 *    -   "phaseName:speciesName"
 *    .
 *
 *  The first two methods of naming may not yield a unique species within
 *  complicated assemblies of Cantera Phases.
 *
 * @param nameStr String name of the species. It may also be the phase name
 *                species name combination, separated by a colon.
 * @return     Returns the index of the species. If the name is not found,
 *             the value of -1 is returned.
 */
size_t Phase::speciesIndex(std::string nameStr) const
{
    std::string pn;
    std::string sn = parseSpeciesName(nameStr, pn);
    if (pn == "" || pn == m_name || pn == m_id) {
        return Constituents::speciesIndex(sn);
    }
    return npos;
}

std::string Phase::speciesSPName(int k) const
{
    std::string sn = Constituents::speciesName(k);
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

/** Set the temperature (K), density (kg/m^3), and mole fractions. */
void Phase::setState_TRX(doublereal t, doublereal dens,
                         const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setDensity(dens);
}

void Phase::setState_TNX(doublereal t, doublereal n,
                         const doublereal* x)
{
    setMoleFractions(x);
    setTemperature(t);
    setMolarDensity(n);
}

/** Set the temperature (K), density (kg/m^3), and mole fractions. */
void Phase::setState_TRX(doublereal t, doublereal dens,
                         compositionMap& x)
{
    setMoleFractionsByName(x);
    setTemperature(t);
    setDensity(dens);
}

/** Set the temperature (K), density (kg/m^3), and mass fractions. */
void Phase::setState_TRY(doublereal t, doublereal dens,
                         const doublereal* y)
{
    setMassFractions(y);
    setTemperature(t);
    setDensity(dens);
}

/** Set the temperature (K), density (kg/m^3), and mass fractions. */
void Phase::setState_TRY(doublereal t, doublereal dens,
                         compositionMap& y)
{
    setMassFractionsByName(y);
    setTemperature(t);
    setDensity(dens);
}

/** Set the temperature (K) and density (kg/m^3) */
void Phase::setState_TR(doublereal t, doublereal rho)
{
    setTemperature(t);
    setDensity(rho);
}

/** Set the temperature (K) and mole fractions.  */
void Phase::setState_TX(doublereal t, doublereal* x)
{
    setTemperature(t);
    setMoleFractions(x);
}

/** Set the temperature (K) and mass fractions.  */
void Phase::setState_TY(doublereal t, doublereal* y)
{
    setTemperature(t);
    setMassFractions(y);
}

/** Set the density (kg/m^3) and mole fractions.  */
void Phase::setState_RX(doublereal rho, doublereal* x)
{
    setMoleFractions(x);
    setDensity(rho);
}

/** Set the density (kg/m^3) and mass fractions.  */
void Phase::setState_RY(doublereal rho, doublereal* y)
{
    setMassFractions(y);
    setDensity(rho);
}

/*
 * Copy the vector of molecular weights into vector weights.
 */
void Phase::getMolecularWeights(vector_fp& weights) const
{
    const vector_fp& mw = Constituents::molecularWeights();
    if (weights.size() < mw.size()) {
        weights.resize(mw.size());
    }
    copy(mw.begin(), mw.end(), weights.begin());
}

/*
 * Copy the vector of molecular weights into array weights.
 * @deprecated
 */
void Phase::getMolecularWeights(int iwt, doublereal* weights) const
{
    const vector_fp& mw = Constituents::molecularWeights();
    copy(mw.begin(), mw.end(), weights);
}

/*
 * Copy the vector of molecular weights into array weights.
 */
void Phase::getMolecularWeights(doublereal* weights) const
{
    const vector_fp& mw = Constituents::molecularWeights();
    copy(mw.begin(), mw.end(), weights);
}

/**
 * Return a const reference to the internal vector of
 * molecular weights.
 */
const vector_fp& Phase::molecularWeights() const
{
    return Constituents::molecularWeights();
}


/**
 * Get the mole fractions by name.
 */
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

/**
 *  Finished adding species, prepare to use them for calculation
 *  of mixture properties.
 */
void Phase::freezeSpecies()
{
    Constituents::freezeSpecies();
    init(Constituents::molecularWeights());
    size_t kk = nSpecies();
    size_t nv = kk + 2;
    m_data.resize(nv,0.0);
    m_data[0] = 300.0;
    m_data[1] = 0.001;
    m_data[2] = 1.0;

    m_kk = nSpecies();
}

void Phase::init(const vector_fp& mw)
{
    m_kk = mw.size();
    m_molwts.resize(m_kk);
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
        /*
         * Some surface phases may define species representing
         * empty sites that have zero molecular weight. Give them
         * a very small molecular weight to avoid dividing by
         * zero.
         */
        if (m_molwts[k] < Tiny) {
            m_molwts[k] = Tiny;
        }
        m_rmolwts[k] = 1.0/m_molwts[k];
    }

    /*
     * Now that we have resized the State object, let's fill it with
     * a valid mass fraction vector that sums to one. The State object
     * should never have a mass fraction vector that doesn't sum to one.
     * We will assume that species 0 has a mass fraction of 1.0 and
     * mass fraction of all other species is 0.0.
     */
    m_y[0] = 1.0;
    m_ym[0] = m_y[0] * m_rmolwts[0];
    m_mmw = 1.0 / m_ym[0];
}

bool Phase::ready() const
{
    return (m_kk > 0 && Constituents::ready());
}

} // namespace Cantera
