/**
 * @file MultiPhase.cpp
 * Definitions for the \link Cantera::MultiPhase MultiPhase\endlink
 * object that is used to set up multiphase equilibrium problems (see \ref equilfunctions).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/ChemEquil.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

MultiPhase::MultiPhase() :
    m_temp(298.15),
    m_press(OneBar),
    m_nel(0),
    m_nsp(0),
    m_init(false),
    m_eloc(npos),
    m_Tmin(1.0),
    m_Tmax(100000.0)
{
}

void MultiPhase::addPhases(MultiPhase& mix)
{
    for (size_t n = 0; n < mix.nPhases(); n++) {
        addPhase(mix.m_phase[n], mix.m_moles[n]);
    }
}

void MultiPhase::addPhases(std::vector<ThermoPhase*>& phases,
                           const vector_fp& phaseMoles)
{
    for (size_t n = 0; n < phases.size(); n++) {
        addPhase(phases[n], phaseMoles[n]);
    }
    init();
}

void MultiPhase::addPhase(ThermoPhase* p, doublereal moles)
{
    if (m_init) {
        throw CanteraError("MultiPhase::addPhase",
                           "phases cannot be added after init() has been called.");
    }

    if (!p->compatibleWithMultiPhase()) {
        throw CanteraError("MultiPhase::addPhase", "Phase '{}'' is not "
            "compatible with MultiPhase equilibrium solver", p->name());
    }

    // save the pointer to the phase object
    m_phase.push_back(p);

    // store its number of moles
    m_moles.push_back(moles);
    m_temp_OK.push_back(true);

    // update the total number of species
    m_nsp += p->nSpecies();

    // determine if this phase has new elements for each new element, add an
    // entry in the map from names to index number + 1:

    // iterate over the elements in this phase
    for (size_t m = 0; m < p->nElements(); m++) {
        string ename = p->elementName(m);

        // if no entry is found for this element name, then it is a new element.
        // In this case, add the name to the list of names, increment the
        // element count, and add an entry to the name->(index+1) map.
        if (m_enamemap.find(ename) == m_enamemap.end()) {
            m_enamemap[ename] = m_nel + 1;
            m_enames.push_back(ename);
            m_atomicNumber.push_back(p->atomicNumber(m));

            // Element 'E' (or 'e') is special. Note its location.
            if (ename == "E" || ename == "e") {
                m_eloc = m_nel;
            }

            m_nel++;
        }
    }

    // If the mixture temperature hasn't been set, then set the temperature and
    // pressure to the values for the phase being added. There is no good way to
    // do this. However, this will be overridden later.
    if (m_temp == 298.15 && p->temperature() > 2.0E-3) {
        m_temp = p->temperature();
        m_press = p->pressure();
    }

    // If this is a solution phase, update the minimum and maximum mixture
    // temperatures. Stoichiometric phases are excluded, since a mixture may
    // define multiple stoichiometric phases, each of which has thermo data
    // valid only over a limited range. For example, a mixture might be defined
    // to contain a phase representing water ice and one representing liquid
    // water, only one of which should be present if the mixture represents an
    // equilibrium state.
    if (p->nSpecies() > 1) {
        m_Tmin = std::max(p->minTemp(), m_Tmin);
        m_Tmax = std::min(p->maxTemp(), m_Tmax);
    }
}

void MultiPhase::init()
{
    if (m_init) {
        return;
    }

    // allocate space for the atomic composition matrix
    m_atoms.resize(m_nel, m_nsp, 0.0);
    m_moleFractions.resize(m_nsp, 0.0);
    m_elemAbundances.resize(m_nel, 0.0);

    // iterate over the elements
    //   -> fill in m_atoms(m,k), m_snames(k), m_spphase(k), m_spstart(ip)
    for (size_t m = 0; m < m_nel; m++) {
        size_t k = 0;
        // iterate over the phases
        for (size_t ip = 0; ip < nPhases(); ip++) {
            ThermoPhase* p = m_phase[ip];
            size_t nsp = p->nSpecies();
            size_t mlocal = p->elementIndex(m_enames[m]);
            for (size_t kp = 0; kp < nsp; kp++) {
                if (mlocal != npos) {
                    m_atoms(m, k) = p->nAtoms(kp, mlocal);
                }
                if (m == 0) {
                    m_snames.push_back(p->speciesName(kp));
                    if (kp == 0) {
                        m_spstart.push_back(m_spphase.size());
                    }
                    m_spphase.push_back(ip);
                }
                k++;
            }
        }
    }

    // set the initial composition within each phase to the
    // mole fractions stored in the phase objects
    m_init = true;
    uploadMoleFractionsFromPhases();
    updatePhases();
}

ThermoPhase& MultiPhase::phase(size_t n)
{
    if (!m_init) {
        init();
    }
    m_phase[n]->setTemperature(m_temp);
    m_phase[n]->setMoleFractions_NoNorm(&m_moleFractions[m_spstart[n]]);
    m_phase[n]->setPressure(m_press);
    return *m_phase[n];
}

void MultiPhase::checkPhaseIndex(size_t m) const
{
    if (m >= nPhases()) {
        throw IndexError("MultiPhase::checkPhaseIndex", "phase", m, nPhases()-1);
    }
}

void MultiPhase::checkPhaseArraySize(size_t mm) const
{
    if (nPhases() > mm) {
        throw ArraySizeError("MultiPhase::checkPhaseIndex", mm, nPhases());
    }
}

doublereal MultiPhase::speciesMoles(size_t k) const
{
    size_t ip = m_spphase[k];
    return m_moles[ip]*m_moleFractions[k];
}

doublereal MultiPhase::elementMoles(size_t m) const
{
    doublereal sum = 0.0;
    for (size_t i = 0; i < nPhases(); i++) {
        double phasesum = 0.0;
        size_t nsp = m_phase[i]->nSpecies();
        for (size_t ik = 0; ik < nsp; ik++) {
            size_t k = speciesIndex(ik, i);
            phasesum += m_atoms(m,k)*m_moleFractions[k];
        }
        sum += phasesum * m_moles[i];
    }
    return sum;
}

doublereal MultiPhase::charge() const
{
    doublereal sum = 0.0;
    for (size_t i = 0; i < nPhases(); i++) {
        sum += phaseCharge(i);
    }
    return sum;
}

size_t MultiPhase::speciesIndex(const std::string& speciesName, const std::string& phaseName)
{
    if (!m_init) {
        init();
    }
    size_t p = phaseIndex(phaseName);
    if (p == npos) {
        throw CanteraError("MultiPhase::speciesIndex", "phase not found: " + phaseName);
    }
    size_t k = m_phase[p]->speciesIndex(speciesName);
    if (k == npos) {
        throw CanteraError("MultiPhase::speciesIndex", "species not found: " + speciesName);
    }
    return m_spstart[p] + k;
}

doublereal MultiPhase::phaseCharge(size_t p) const
{
    doublereal phasesum = 0.0;
    size_t nsp = m_phase[p]->nSpecies();
    for (size_t ik = 0; ik < nsp; ik++) {
        size_t k = speciesIndex(ik, p);
        phasesum += m_phase[p]->charge(ik)*m_moleFractions[k];
    }
    return Faraday*phasesum*m_moles[p];
}

void MultiPhase::getChemPotentials(doublereal* mu) const
{
    updatePhases();
    size_t loc = 0;
    for (size_t i = 0; i < nPhases(); i++) {
        m_phase[i]->getChemPotentials(mu + loc);
        loc += m_phase[i]->nSpecies();
    }
}

void MultiPhase::getValidChemPotentials(doublereal not_mu,
                                        doublereal* mu, bool standard) const
{
    updatePhases();
    // iterate over the phases
    size_t loc = 0;
    for (size_t i = 0; i < nPhases(); i++) {
        if (tempOK(i) || m_phase[i]->nSpecies() > 1) {
            if (!standard) {
                m_phase[i]->getChemPotentials(mu + loc);
            } else {
                m_phase[i]->getStandardChemPotentials(mu + loc);
            }
        } else {
            fill(mu + loc, mu + loc + m_phase[i]->nSpecies(), not_mu);
        }
        loc += m_phase[i]->nSpecies();
    }
}

bool MultiPhase::solutionSpecies(size_t k) const
{
    if (m_phase[m_spphase[k]]->nSpecies() > 1) {
        return true;
    } else {
        return false;
    }
}

doublereal MultiPhase::gibbs() const
{
    doublereal sum = 0.0;
    updatePhases();
    for (size_t i = 0; i < nPhases(); i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->gibbs_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::enthalpy() const
{
    doublereal sum = 0.0;
    updatePhases();
    for (size_t i = 0; i < nPhases(); i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->enthalpy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::IntEnergy() const
{
    doublereal sum = 0.0;
    updatePhases();
    for (size_t i = 0; i < nPhases(); i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->intEnergy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::entropy() const
{
    doublereal sum = 0.0;
    updatePhases();
    for (size_t i = 0; i < nPhases(); i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->entropy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::cp() const
{
    doublereal sum = 0.0;
    updatePhases();
    for (size_t i = 0; i < nPhases(); i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->cp_mole() * m_moles[i];
        }
    }
    return sum;
}

void MultiPhase::setPhaseMoleFractions(const size_t n, const doublereal* const x)
{
    if (!m_init) {
        init();
    }
    ThermoPhase* p = m_phase[n];
    p->setState_TPX(m_temp, m_press, x);
    size_t istart = m_spstart[n];
    for (size_t k = 0; k < p->nSpecies(); k++) {
        m_moleFractions[istart+k] = x[k];
    }
}

void MultiPhase::setMolesByName(const compositionMap& xMap)
{
    size_t kk = nSpecies();
    vector_fp moles(kk, 0.0);
    for (size_t k = 0; k < kk; k++) {
        moles[k] = std::max(getValue(xMap, speciesName(k), 0.0), 0.0);
    }
    setMoles(moles.data());
}

void MultiPhase::setMolesByName(const std::string& x)
{
    // build the composition map from the string, and then set the moles.
    compositionMap xx = parseCompString(x, m_snames);
    setMolesByName(xx);
}

void MultiPhase::getMoles(doublereal* molNum) const
{
    // First copy in the mole fractions
    copy(m_moleFractions.begin(), m_moleFractions.end(), molNum);
    doublereal* dtmp = molNum;
    for (size_t ip = 0; ip < nPhases(); ip++) {
        doublereal phasemoles = m_moles[ip];
        ThermoPhase* p = m_phase[ip];
        size_t nsp = p->nSpecies();
        for (size_t ik = 0; ik < nsp; ik++) {
            *(dtmp++) *= phasemoles;
        }
    }
}

void MultiPhase::setMoles(const doublereal* n)
{
    if (!m_init) {
        init();
    }
    size_t loc = 0;
    size_t k = 0;
    for (size_t ip = 0; ip < nPhases(); ip++) {
        ThermoPhase* p = m_phase[ip];
        size_t nsp = p->nSpecies();
        double phasemoles = 0.0;
        for (size_t ik = 0; ik < nsp; ik++) {
            phasemoles += n[k];
            k++;
        }
        m_moles[ip] = phasemoles;
        if (nsp > 1) {
            if (phasemoles > 0.0) {
                p->setState_TPX(m_temp, m_press, n + loc);
                p->getMoleFractions(&m_moleFractions[loc]);
            } else {
                p->getMoleFractions(&m_moleFractions[loc]);
            }
        } else {
            m_moleFractions[loc] = 1.0;
        }
        loc += nsp;
    }
}

void MultiPhase::addSpeciesMoles(const int indexS, const doublereal addedMoles)
{
    vector_fp tmpMoles(m_nsp, 0.0);
    getMoles(tmpMoles.data());
    tmpMoles[indexS] += addedMoles;
    tmpMoles[indexS] = std::max(tmpMoles[indexS], 0.0);
    setMoles(tmpMoles.data());
}

void MultiPhase::setState_TP(const doublereal T, const doublereal Pres)
{
    if (!m_init) {
        init();
    }
    m_temp = T;
    m_press = Pres;
    updatePhases();
}

void MultiPhase::setState_TPMoles(const doublereal T, const doublereal Pres,
                                  const doublereal* n)
{
    m_temp = T;
    m_press = Pres;
    setMoles(n);
}

void MultiPhase::getElemAbundances(doublereal* elemAbundances) const
{
    calcElemAbundances();
    for (size_t eGlobal = 0; eGlobal < m_nel; eGlobal++) {
        elemAbundances[eGlobal] = m_elemAbundances[eGlobal];
    }
}

void MultiPhase::calcElemAbundances() const
{
    size_t loc = 0;
    doublereal spMoles;
    for (size_t eGlobal = 0; eGlobal < m_nel; eGlobal++) {
        m_elemAbundances[eGlobal] = 0.0;
    }
    for (size_t ip = 0; ip < nPhases(); ip++) {
        ThermoPhase* p = m_phase[ip];
        size_t nspPhase = p->nSpecies();
        doublereal phasemoles = m_moles[ip];
        for (size_t ik = 0; ik < nspPhase; ik++) {
            size_t kGlobal = loc + ik;
            spMoles = m_moleFractions[kGlobal] * phasemoles;
            for (size_t eGlobal = 0; eGlobal < m_nel; eGlobal++) {
                m_elemAbundances[eGlobal] += m_atoms(eGlobal, kGlobal) * spMoles;
            }
        }
        loc += nspPhase;
    }
}

doublereal MultiPhase::volume() const
{
    doublereal sum = 0;
    for (size_t i = 0; i < nPhases(); i++) {
        double vol = 1.0/m_phase[i]->molarDensity();
        sum += m_moles[i] * vol;
    }
    return sum;
}

double MultiPhase::equilibrate_MultiPhaseEquil(int XY, doublereal err,
                                               int maxsteps, int maxiter,
                                               int loglevel)
{
    bool strt = false;
    doublereal dta = 0.0;
    if (!m_init) {
        init();
    }

    if (XY == TP) {
        // create an equilibrium manager
        MultiPhaseEquil e(this);
        return e.equilibrate(XY, err, maxsteps, loglevel);
    } else if (XY == HP) {
        double h0 = enthalpy();
        double Tlow = 0.5*m_Tmin; // lower bound on T
        double Thigh = 2.0*m_Tmax; // upper bound on T
        doublereal Hlow = Undef, Hhigh = Undef;
        for (int n = 0; n < maxiter; n++) {
            // if 'strt' is false, the current composition will be used as
            // the starting estimate; otherwise it will be estimated
            MultiPhaseEquil e(this, strt);
            // start with a loose error tolerance, but tighten it as we get
            // close to the final temperature

            try {
                e.equilibrate(TP, err, maxsteps, loglevel);
                double hnow = enthalpy();
                // the equilibrium enthalpy monotonically increases with T;
                // if the current value is below the target, the we know the
                // current temperature is too low. Set
                if (hnow < h0) {
                    if (m_temp > Tlow) {
                        Tlow = m_temp;
                        Hlow = hnow;
                    }
                } else {
                    // the current enthalpy is greater than the target;
                    // therefore the current temperature is too high.
                    if (m_temp < Thigh) {
                        Thigh = m_temp;
                        Hhigh = hnow;
                    }
                }
                double dt;
                if (Hlow != Undef && Hhigh != Undef) {
                    double cpb = (Hhigh - Hlow)/(Thigh - Tlow);
                    dt = (h0 - hnow)/cpb;
                    dta = fabs(dt);
                    double dtmax = 0.5*fabs(Thigh - Tlow);
                    if (dta > dtmax) {
                        dt *= dtmax/dta;
                    }
                } else {
                    double tnew = sqrt(Tlow*Thigh);
                    dt = tnew - m_temp;
                }

                double herr = fabs((h0 - hnow)/h0);

                if (herr < err) {
                    return err;
                }
                double tnew = m_temp + dt;
                if (tnew < 0.0) {
                    tnew = 0.5*m_temp;
                }
                setTemperature(tnew);

                // if the size of Delta T is not too large, use
                // the current composition as the starting estimate
                if (dta < 100.0) {
                    strt = false;
                }

            } catch (CanteraError&) {
                if (!strt) {
                    strt = true;
                } else {
                    double tnew = 0.5*(m_temp + Thigh);
                    if (fabs(tnew - m_temp) < 1.0) {
                        tnew = m_temp + 1.0;
                    }
                    setTemperature(tnew);
                }
            }
        }
        throw CanteraError("MultiPhase::equilibrate_MultiPhaseEquil",
                           "No convergence for T");
    } else if (XY == SP) {
        double s0 = entropy();
        double Tlow = 1.0; // lower bound on T
        double Thigh = 1.0e6; // upper bound on T
        for (int n = 0; n < maxiter; n++) {
            MultiPhaseEquil e(this, strt);

            try {
                e.equilibrate(TP, err, maxsteps, loglevel);
                double snow = entropy();
                if (snow < s0) {
                    Tlow = std::max(Tlow, m_temp);
                } else {
                    Thigh = std::min(Thigh, m_temp);
                }
                double dt = (s0 - snow)*m_temp/cp();
                double dtmax = 0.5*fabs(Thigh - Tlow);
                dtmax = (dtmax > 500.0 ? 500.0 : dtmax);
                dta = fabs(dt);
                if (dta > dtmax) {
                    dt *= dtmax/dta;
                }
                if (dta < 1.0e-4) {
                    return err;
                }
                double tnew = m_temp + dt;
                setTemperature(tnew);

                // if the size of Delta T is not too large, use
                // the current composition as the starting estimate
                if (dta < 100.0) {
                    strt = false;
                }
            } catch (CanteraError&) {
                if (!strt) {
                    strt = true;
                } else {
                    double tnew = 0.5*(m_temp + Thigh);
                    setTemperature(tnew);
                }
            }
        }
        throw CanteraError("MultiPhase::equilibrate_MultiPhaseEquil",
                           "No convergence for T");
    } else if (XY == TV) {
        doublereal v0 = volume();
        bool start = true;
        for (int n = 0; n < maxiter; n++) {
            double pnow = pressure();
            MultiPhaseEquil e(this, start);
            start = false;

            e.equilibrate(TP, err, maxsteps, loglevel);
            double vnow = volume();
            double verr = fabs((v0 - vnow)/v0);

            if (verr < err) {
                return err;
            }
            // find dV/dP
            setPressure(pnow*1.01);
            double dVdP = (volume() - vnow)/(0.01*pnow);
            setPressure(pnow + 0.5*(v0 - vnow)/dVdP);
        }
    } else {
        throw CanteraError("MultiPhase::equilibrate_MultiPhaseEquil",
                           "unknown option");
    }
    return -1.0;
}

void MultiPhase::equilibrate(const std::string& XY, const std::string& solver,
                             double rtol, int max_steps, int max_iter,
                             int estimate_equil, int log_level)
{
    // Save the initial state so that it can be restored in case one of the
    // solvers fails
    vector_fp initial_moleFractions = m_moleFractions;
    vector_fp initial_moles = m_moles;
    double initial_T = m_temp;
    double initial_P = m_press;
    int ixy = _equilflag(XY.c_str());
    if (solver == "auto" || solver == "vcs") {
        try {
            debuglog("Trying VCS equilibrium solver\n", log_level);
            vcs_MultiPhaseEquil eqsolve(this, log_level-1);
            int ret = eqsolve.equilibrate(ixy, estimate_equil, log_level-1,
                                          rtol, max_steps);
            if (ret) {
                throw CanteraError("MultiPhase::equilibrate",
                    "VCS solver failed. Return code: {}", ret);
            }
            debuglog("VCS solver succeeded\n", log_level);
            return;
        } catch (std::exception& err) {
            debuglog("VCS solver failed.\n", log_level);
            debuglog(err.what(), log_level);
            m_moleFractions = initial_moleFractions;
            m_moles = initial_moles;
            m_temp = initial_T;
            m_press = initial_P;
            updatePhases();
            if (solver == "auto") {
            } else {
                throw;
            }
        }
    }

    if (solver == "auto" || solver == "gibbs") {
        try {
            debuglog("Trying MultiPhaseEquil (Gibbs) equilibrium solver\n",
                     log_level);
            equilibrate_MultiPhaseEquil(ixy, rtol, max_steps, max_iter,
                                        log_level-1);
            debuglog("MultiPhaseEquil solver succeeded\n", log_level);
            return;
        } catch (std::exception& err) {
            debuglog("MultiPhaseEquil solver failed.\n", log_level);
            debuglog(err.what(), log_level);
            m_moleFractions = initial_moleFractions;
            m_moles = initial_moles;
            m_temp = initial_T;
            m_press = initial_P;
            updatePhases();
            throw;
        }
    }

    if (solver != "auto") {
        throw CanteraError("MultiPhase::equilibrate",
            "Invalid solver specified: '" + solver + "'");
    }
}

void MultiPhase::setTemperature(const doublereal T)
{
    if (!m_init) {
        init();
    }
    m_temp = T;
    updatePhases();
}

void MultiPhase::checkElementIndex(size_t m) const
{
    if (m >= m_nel) {
        throw IndexError("MultiPhase::checkElementIndex", "elements", m, m_nel-1);
    }
}

void MultiPhase::checkElementArraySize(size_t mm) const
{
    if (m_nel > mm) {
        throw ArraySizeError("MultiPhase::checkElementArraySize", mm, m_nel);
    }
}

std::string MultiPhase::elementName(size_t m) const
{
    return m_enames[m];
}

size_t MultiPhase::elementIndex(const std::string& name) const
{
    for (size_t e = 0; e < m_nel; e++) {
        if (m_enames[e] == name) {
            return e;
        }
    }
    return npos;
}

void MultiPhase::checkSpeciesIndex(size_t k) const
{
    if (k >= m_nsp) {
        throw IndexError("MultiPhase::checkSpeciesIndex", "species", k, m_nsp-1);
    }
}

void MultiPhase::checkSpeciesArraySize(size_t kk) const
{
    if (m_nsp > kk) {
        throw ArraySizeError("MultiPhase::checkSpeciesArraySize", kk, m_nsp);
    }
}

std::string MultiPhase::speciesName(const size_t k) const
{
    return m_snames[k];
}

doublereal MultiPhase::nAtoms(const size_t kGlob, const size_t mGlob) const
{
    return m_atoms(mGlob, kGlob);
}

void MultiPhase::getMoleFractions(doublereal* const x) const
{
    std::copy(m_moleFractions.begin(), m_moleFractions.end(), x);
}

std::string MultiPhase::phaseName(const size_t iph) const
{
    const ThermoPhase* tptr = m_phase[iph];
    return tptr->id();
}

int MultiPhase::phaseIndex(const std::string& pName) const
{
    for (int iph = 0; iph < (int) nPhases(); iph++) {
        if (m_phase[iph]->id() == pName) {
            return iph;
        }
    }
    return -1;
}

doublereal MultiPhase::phaseMoles(const size_t n) const
{
    return m_moles[n];
}

void MultiPhase::setPhaseMoles(const size_t n, const doublereal moles)
{
    m_moles[n] = moles;
}

size_t MultiPhase::speciesPhaseIndex(const size_t kGlob) const
{
    return m_spphase[kGlob];
}

doublereal MultiPhase::moleFraction(const size_t kGlob) const
{
    return m_moleFractions[kGlob];
}

bool MultiPhase::tempOK(const size_t p) const
{
    return m_temp_OK[p];
}

void MultiPhase::uploadMoleFractionsFromPhases()
{
    size_t loc = 0;
    for (size_t ip = 0; ip < nPhases(); ip++) {
        ThermoPhase* p = m_phase[ip];
        p->getMoleFractions(&m_moleFractions[loc]);
        loc += p->nSpecies();
    }
    calcElemAbundances();
}

void MultiPhase::updatePhases() const
{
    size_t loc = 0;
    for (size_t p = 0; p < nPhases(); p++) {
        m_phase[p]->setState_TPX(m_temp, m_press, &m_moleFractions[loc]);
        loc += m_phase[p]->nSpecies();
        m_temp_OK[p] = true;
        if (m_temp < m_phase[p]->minTemp() || m_temp > m_phase[p]->maxTemp()) {
            m_temp_OK[p] = false;
        }
    }
}
}
