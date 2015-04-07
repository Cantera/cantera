/**
 * @file MultiPhase.cpp
 * Definitions for the \link Cantera::MultiPhase MultiPhase\endlink
 * object that is used to set up multiphase equilibrium problems (see \ref equilfunctions).
 */
#include "cantera/equil/MultiPhase.h"
#include "cantera/equil/MultiPhaseEquil.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

MultiPhase::MultiPhase() :
    m_np(0),
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

MultiPhase::MultiPhase(const MultiPhase& right) :
    m_np(0),
    m_temp(298.15),
    m_press(OneBar),
    m_nel(0),
    m_nsp(0),
    m_init(false),
    m_eloc(npos),
    m_Tmin(1.0),
    m_Tmax(100000.0)
{
    operator=(right);
}

MultiPhase::~MultiPhase()
{
}

MultiPhase& MultiPhase::operator=(const MultiPhase& right)
{
    if (&right != this) {
        m_moles = right.m_moles;
        // shallow copy of phase pointers
        m_phase = right.m_phase;
        m_atoms = right.m_atoms;
        m_moleFractions = right.m_moleFractions;
        m_spphase = right.m_spphase;
        m_spstart = right.m_spstart;
        m_enames = right.m_enames;
        m_enamemap = right.m_enamemap;
        m_np = right.m_np;
        m_temp = right.m_temp;
        m_press = right.m_press;
        m_nel = right.m_nel;
        m_nsp = right.m_nsp;
        m_init = right.m_init;
        m_eloc = right.m_eloc;
        m_temp_OK = right.m_temp_OK;
        m_Tmin = right.m_Tmin;
        m_Tmax = right.m_Tmax;
        m_elemAbundances = right.m_elemAbundances;
    }
    return *this;
}

void MultiPhase::
addPhases(MultiPhase& mix)
{
    size_t n;
    for (n = 0; n < mix.m_np; n++) {
        addPhase(mix.m_phase[n], mix.m_moles[n]);
    }
}

void MultiPhase::
addPhases(std::vector<ThermoPhase*>& phases, const vector_fp& phaseMoles)
{
    size_t np = phases.size();
    size_t n;
    for (n = 0; n < np; n++) {
        addPhase(phases[n], phaseMoles[n]);
    }
    init();
}

void MultiPhase::
addPhase(ThermoPhase* p, doublereal moles)
{
    if (m_init) {
        throw CanteraError("addPhase",
                           "phases cannot be added after init() has been called.");
    }

    // save the pointer to the phase object
    m_phase.push_back(p);

    // store its number of moles
    m_moles.push_back(moles);
    m_temp_OK.push_back(true);

    // update the number of phases and the total number of
    // species
    m_np = m_phase.size();
    m_nsp += p->nSpecies();

    // determine if this phase has new elements
    // for each new element, add an entry in the map
    // from names to index number + 1:

    string ename;
    // iterate over the elements in this phase
    size_t m, nel = p->nElements();
    for (m = 0; m < nel; m++) {
        ename = p->elementName(m);

        // if no entry is found for this element name, then
        // it is a new element. In this case, add the name
        // to the list of names, increment the element count,
        // and add an entry to the name->(index+1) map.
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

    // If the mixture temperature hasn't been set, then set the
    // temperature and pressure to the values for the phase being
    // added. There is no good way to do this. However, this will be overridden later.
    if (m_temp == 298.15 && p->temperature() > 2.0E-3) {
        m_temp = p->temperature();
        m_press = p->pressure();
    }

    // If this is a solution phase, update the minimum and maximum
    // mixture temperatures. Stoichiometric phases are excluded,
    // since a mixture may define multiple stoichiometric phases,
    // each of which has thermo data valid only over a limited
    // range. For example, a mixture might be defined to contain a
    // phase representing water ice and one representing liquid
    // water, only one of which should be present if the mixture
    // represents an equilibrium state.
    if (p->nSpecies() > 1) {
        double t = p->minTemp();
        if (t > m_Tmin) {
            m_Tmin = t;
        }
        t = p->maxTemp();
        if (t < m_Tmax) {
            m_Tmax = t;
        }
    }
}

void MultiPhase::init()
{
    if (m_init) {
        return;
    }
    size_t ip, kp, k = 0, nsp, m;
    size_t mlocal;
    string sym;

    // allocate space for the atomic composition matrix
    m_atoms.resize(m_nel, m_nsp, 0.0);
    m_moleFractions.resize(m_nsp, 0.0);
    m_elemAbundances.resize(m_nel, 0.0);

    // iterate over the elements
    //   -> fill in m_atoms(m,k), m_snames(k), m_spphase(k),
    //              m_sptart(ip)
    for (m = 0; m < m_nel; m++) {
        sym = m_enames[m];
        k = 0;
        // iterate over the phases
        for (ip = 0; ip < m_np; ip++) {
            ThermoPhase* p = m_phase[ip];
            nsp = p->nSpecies();
            mlocal = p->elementIndex(sym);
            for (kp = 0; kp < nsp; kp++) {
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

    if (m_eloc != npos) {
        doublereal esum;
        for (k = 0; k < m_nsp; k++) {
            esum = 0.0;
            for (m = 0; m < m_nel; m++) {
                if (m != m_eloc) {
                    esum += m_atoms(m,k) * m_atomicNumber[m];
                }
            }
            //m_atoms(m_eloc, k) += esum;
        }
    }

    /// set the initial composition within each phase to the
    /// mole fractions stored in the phase objects
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
    m_phase[n]->setMoleFractions_NoNorm(DATA_PTR(m_moleFractions) + m_spstart[n]);
    m_phase[n]->setPressure(m_press);
    return *m_phase[n];
}

void MultiPhase::checkPhaseIndex(size_t m) const
{
    if (m >= nPhases()) {
        throw IndexError("checkPhaseIndex", "phase", m, nPhases()-1);
    }
}

void MultiPhase::checkPhaseArraySize(size_t mm) const
{
    if (nPhases() > mm) {
        throw ArraySizeError("checkPhaseIndex", mm, nPhases());
    }
}

doublereal MultiPhase::speciesMoles(size_t k) const
{
    size_t ip = m_spphase[k];
    return m_moles[ip]*m_moleFractions[k];
}

doublereal MultiPhase::elementMoles(size_t m) const
{
    doublereal sum = 0.0, phasesum;
    size_t i, k = 0, ik, nsp;
    for (i = 0; i < m_np; i++) {
        phasesum = 0.0;
        nsp = m_phase[i]->nSpecies();
        for (ik = 0; ik < nsp; ik++) {
            k = speciesIndex(ik, i);
            phasesum += m_atoms(m,k)*m_moleFractions[k];
        }
        sum += phasesum * m_moles[i];
    }
    return sum;
}

doublereal MultiPhase::charge() const
{
    doublereal sum = 0.0;
    size_t i;
    for (i = 0; i < m_np; i++) {
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
    size_t ik, k, nsp = m_phase[p]->nSpecies();
    for (ik = 0; ik < nsp; ik++) {
        k = speciesIndex(ik, p);
        phasesum += m_phase[p]->charge(ik)*m_moleFractions[k];
    }
    return Faraday*phasesum*m_moles[p];
}

void MultiPhase::getChemPotentials(doublereal* mu) const
{
    size_t i, loc = 0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
        m_phase[i]->getChemPotentials(mu + loc);
        loc += m_phase[i]->nSpecies();
    }
}

void MultiPhase::getValidChemPotentials(doublereal not_mu,
                                        doublereal* mu, bool standard) const
{
    size_t i, loc = 0;

    updatePhases();
    // iterate over the phases
    for (i = 0; i < m_np; i++) {
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
    size_t i;
    doublereal sum = 0.0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->gibbs_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::enthalpy() const
{
    size_t i;
    doublereal sum = 0.0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->enthalpy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::IntEnergy() const
{
    size_t i;
    doublereal sum = 0.0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->intEnergy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::entropy() const
{
    size_t i;
    doublereal sum = 0.0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
        if (m_moles[i] > 0.0) {
            sum += m_phase[i]->entropy_mole() * m_moles[i];
        }
    }
    return sum;
}

doublereal MultiPhase::cp() const
{
    size_t i;
    doublereal sum = 0.0;
    updatePhases();
    for (i = 0; i < m_np; i++) {
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

void MultiPhase::setMolesByName(compositionMap& xMap)
{
    size_t kk = nSpecies();
    doublereal x;
    vector_fp moles(kk, 0.0);
    for (size_t k = 0; k < kk; k++) {
        x = xMap[speciesName(k)];
        if (x > 0.0) {
            moles[k] = x;
        }
    }
    setMoles(DATA_PTR(moles));
}

void MultiPhase::setMolesByName(const std::string& x)
{
    // build the composition map from the string, and then set the moles.
    compositionMap xx = parseCompString(x, m_snames);
    setMolesByName(xx);
}

void MultiPhase::getMoles(doublereal* molNum) const
{
    /*
     * First copy in the mole fractions
     */
    copy(m_moleFractions.begin(), m_moleFractions.end(), molNum);
    size_t ik;
    doublereal* dtmp = molNum;
    for (size_t ip = 0; ip < m_np; ip++) {
        doublereal phasemoles = m_moles[ip];
        ThermoPhase* p = m_phase[ip];
        size_t nsp = p->nSpecies();
        for (ik = 0; ik < nsp; ik++) {
            *(dtmp++) *= phasemoles;
        }
    }
}

void MultiPhase::setMoles(const doublereal* n)
{
    if (!m_init) {
        init();
    }
    size_t ip, loc = 0;
    size_t ik, k = 0, nsp;
    doublereal phasemoles;
    for (ip = 0; ip < m_np; ip++) {
        ThermoPhase* p = m_phase[ip];
        nsp = p->nSpecies();
        phasemoles = 0.0;
        for (ik = 0; ik < nsp; ik++) {
            phasemoles += n[k];
            k++;
        }
        m_moles[ip] = phasemoles;
        if (nsp > 1) {
            if (phasemoles > 0.0) {
                p->setState_TPX(m_temp, m_press, n + loc);
                p->getMoleFractions(DATA_PTR(m_moleFractions) + loc);
            } else {
                p->getMoleFractions(DATA_PTR(m_moleFractions) + loc);
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
    getMoles(DATA_PTR(tmpMoles));
    tmpMoles[indexS] += addedMoles;
    if (tmpMoles[indexS] < 0.0) {
        tmpMoles[indexS] = 0.0;
    }
    setMoles(DATA_PTR(tmpMoles));
}

void MultiPhase::setState_TP(const doublereal T, const doublereal Pres)
{
    if (!m_init) {
        init();
    }
    m_temp  = T;
    m_press = Pres;
    updatePhases();
}

void MultiPhase::setState_TPMoles(const doublereal T, const doublereal Pres,
                                  const doublereal* n)
{
    m_temp  = T;
    m_press = Pres;
    setMoles(n);
}

void MultiPhase::getElemAbundances(doublereal* elemAbundances) const
{
    size_t eGlobal;
    calcElemAbundances();
    for (eGlobal = 0; eGlobal < m_nel; eGlobal++) {
        elemAbundances[eGlobal] = m_elemAbundances[eGlobal];
    }
}

void MultiPhase::calcElemAbundances() const
{
    size_t loc = 0;
    size_t eGlobal;
    size_t ik, kGlobal;
    doublereal spMoles;
    for (eGlobal = 0; eGlobal < m_nel; eGlobal++) {
        m_elemAbundances[eGlobal] = 0.0;
    }
    for (size_t ip = 0; ip < m_np; ip++) {
        ThermoPhase* p = m_phase[ip];
        size_t nspPhase = p->nSpecies();
        doublereal phasemoles = m_moles[ip];
        for (ik = 0; ik < nspPhase; ik++) {
            kGlobal = loc + ik;
            spMoles = m_moleFractions[kGlobal] * phasemoles;
            for (eGlobal = 0; eGlobal < m_nel; eGlobal++) {
                m_elemAbundances[eGlobal] += m_atoms(eGlobal, kGlobal) * spMoles;
            }
        }
        loc += nspPhase;
    }
}

doublereal MultiPhase::volume() const
{
    int i;
    doublereal sum = 0;
    for (i = 0; i < int(m_np); i++) {
        double vol = 1.0/m_phase[i]->molarDensity();
        sum += m_moles[i] * vol;
    }
    return sum;
}

doublereal MultiPhase::equilibrate(int XY, doublereal err,
                                   int maxsteps, int maxiter, int loglevel)
{
    bool strt = false;
    doublereal dt;
    doublereal h0;
    int n;
    doublereal hnow, herr = 1.0;
    doublereal snow, serr = 1.0, s0;
    doublereal Tlow = -1.0, Thigh = -1.0;
    doublereal Hlow = Undef, Hhigh = Undef, tnew;
    doublereal dta=0.0, dtmax, cpb;
    MultiPhaseEquil* e = 0;

    if (!m_init) {
        init();
    }
    if (loglevel > 0) {
        beginLogGroup("MultiPhase::equilibrate", loglevel);
    }

    if (XY == TP) {
        if (loglevel > 0) {
            addLogEntry("problem type","fixed T,P");
            addLogEntry("Temperature",temperature());
            addLogEntry("Pressure", pressure());
        }

        // create an equilibrium manager
        e = new MultiPhaseEquil(this);
        try {
            e->equilibrate(XY, err, maxsteps, loglevel);
        } catch (CanteraError& err) {
            err.save();
            if (loglevel > 0) {
                endLogGroup();
            }
            delete e;
            e = 0;
            throw err;
        }
        goto done;
    }

    else if (XY == HP) {
        h0 = enthalpy();
        Tlow = 0.5*m_Tmin;      // lower bound on T
        Thigh = 2.0*m_Tmax;     // upper bound on T
        if (loglevel > 0) {
            addLogEntry("problem type","fixed H,P");
            addLogEntry("H target",fp2str(h0));
        }
        for (n = 0; n < maxiter; n++) {

            // if 'strt' is false, the current composition will be used as
            // the starting estimate; otherwise it will be estimated
            //                if (e) {
            //    cout << "e should be zero, but it is not!" << endl;
            //    delete e;
            // }
            e = new MultiPhaseEquil(this, strt);
            // start with a loose error tolerance, but tighten it as we get
            // close to the final temperature
            if (loglevel > 0) {
                beginLogGroup("iteration "+int2str(n));
            }

            try {
                e->equilibrate(TP, err, maxsteps, loglevel);
                hnow = enthalpy();
                // the equilibrium enthalpy monotonically increases with T;
                // if the current value is below the target, the we know the
                // current temperature is too low. Set
                if (hnow < h0) {
                    if (m_temp > Tlow) {
                        Tlow = m_temp;
                        Hlow = hnow;
                    }
                }
                // the current enthalpy is greater than the target; therefore the
                // current temperature is too high.
                else {
                    if (m_temp < Thigh) {
                        Thigh = m_temp;
                        Hhigh = hnow;
                    }
                }
                if (Hlow != Undef && Hhigh != Undef) {
                    cpb = (Hhigh - Hlow)/(Thigh - Tlow);
                    dt = (h0 - hnow)/cpb;
                    dta = fabs(dt);
                    dtmax = 0.5*fabs(Thigh - Tlow);
                    if (dta > dtmax) {
                        dt *= dtmax/dta;
                    }
                } else {
                    tnew = sqrt(Tlow*Thigh);
                    dt = tnew - m_temp;
                    //cpb = cp();
                }

                herr = fabs((h0 - hnow)/h0);
                if (loglevel > 0) {
                    addLogEntry("T",fp2str(temperature()));
                    addLogEntry("H",fp2str(hnow));
                    addLogEntry("H rel error",fp2str(herr));
                    addLogEntry("lower T bound",fp2str(Tlow));
                    addLogEntry("upper T bound",fp2str(Thigh));
                    endLogGroup(); // iteration
                }


                if (herr < err) { // || dta < 1.0e-4) {
                    if (loglevel > 0) {
                        addLogEntry("T iterations",int2str(n));
                        addLogEntry("Final T",fp2str(temperature()));
                        addLogEntry("H rel error",fp2str(herr));
                    }
                    goto done;
                }
                tnew = m_temp + dt;
                if (tnew < 0.0) {
                    tnew = 0.5*m_temp;
                }
                //dta = fabs(tnew - m_temp);
                setTemperature(tnew);

                // if the size of Delta T is not too large, use
                // the current composition as the starting estimate
                if (dta < 100.0) {
                    strt = false;
                }

            }

            catch (CanteraError& err) {
                err.save();
                if (!strt) {
                    if (loglevel > 0)
                        addLogEntry("no convergence",
                                    "try estimating starting composition");
                    strt = true;
                } else {
                    tnew = 0.5*(m_temp + Thigh);
                    if (fabs(tnew - m_temp) < 1.0) {
                        tnew = m_temp + 1.0;
                    }
                    setTemperature(tnew);
                    if (loglevel > 0)
                        addLogEntry("no convergence",
                                    "trying T = "+fp2str(m_temp));
                }
                if (loglevel > 0) {
                    endLogGroup();
                }
            }
            delete e;
            e = 0;
        }
        if (loglevel > 0) {
            addLogEntry("reached max number of T iterations",int2str(maxiter));
            endLogGroup();
        }
        throw CanteraError("MultiPhase::equilibrate",
                           "No convergence for T");
    } else if (XY == SP) {
        s0 = entropy();
        Tlow = 1.0; // m_Tmin;      // lower bound on T
        Thigh = 1.0e6; // m_Tmax;   // upper bound on T
        if (loglevel > 0) {
            addLogEntry("problem type","fixed S,P");
            addLogEntry("S target",fp2str(s0));
            addLogEntry("min T",fp2str(Tlow));
            addLogEntry("max T",fp2str(Thigh));
        }
        for (n = 0; n < maxiter; n++) {
            delete e;
            e = new MultiPhaseEquil(this, strt);
            if (loglevel > 0) {
                beginLogGroup("iteration "+int2str(n));
            }

            try {
                e->equilibrate(TP, err, maxsteps, loglevel);
                snow = entropy();
                if (snow < s0) {
                    if (m_temp > Tlow) {
                        Tlow = m_temp;
                    }
                } else {
                    if (m_temp < Thigh) {
                        Thigh = m_temp;
                    }
                }
                serr = fabs((s0 - snow)/s0);
                if (loglevel > 0) {
                    addLogEntry("T",fp2str(temperature()));
                    addLogEntry("S",fp2str(snow));
                    addLogEntry("S rel error",fp2str(serr));
                    endLogGroup();
                }
                dt = (s0 - snow)*m_temp/cp();
                dtmax = 0.5*fabs(Thigh - Tlow);
                dtmax = (dtmax > 500.0 ? 500.0 : dtmax);
                dta = fabs(dt);
                if (dta > dtmax) {
                    dt *= dtmax/dta;
                }
                if (herr < err || dta < 1.0e-4) {
                    if (loglevel > 0) {
                        addLogEntry("T iterations",int2str(n));
                        addLogEntry("Final T",fp2str(temperature()));
                        addLogEntry("S rel error",fp2str(serr));
                    }
                    goto done;
                }
                tnew = m_temp + dt;
                setTemperature(tnew);

                // if the size of Delta T is not too large, use
                // the current composition as the starting estimate
                if (dta < 100.0) {
                    strt = false;
                }
            }

            catch (CanteraError& err) {
                err.save();
                if (!strt) {
                    if (loglevel > 0) {
                        addLogEntry("no convergence",
                                    "setting strt to True");
                    }
                    strt = true;
                } else {
                    tnew = 0.5*(m_temp + Thigh);
                    setTemperature(tnew);
                    if (loglevel > 0) {
                        addLogEntry("no convergence",
                                    "trying T = "+fp2str(m_temp));
                    }
                }
                if (loglevel > 0) {
                    endLogGroup();
                }
            }
            delete e;
            e = 0;
        }
        if (loglevel > 0) {
            addLogEntry("reached max number of T iterations",int2str(maxiter));
            endLogGroup();
        }
        throw CanteraError("MultiPhase::equilibrate",
                           "No convergence for T");
    } else if (XY == TV) {
        addLogEntry("problem type","fixed T, V");
        //            doublereal dt = 1.0e3;
        doublereal v0 = volume();
        doublereal dVdP;
        int n;
        bool start = true;
        doublereal vnow, pnow, verr;
        for (n = 0; n < maxiter; n++) {
            pnow = pressure();
            MultiPhaseEquil e(this, start);
            start = false;
            beginLogGroup("iteration "+int2str(n));

            e.equilibrate(TP, err, maxsteps, loglevel);
            vnow = volume();
            verr = fabs((v0 - vnow)/v0);
            addLogEntry("P",fp2str(pressure()));
            addLogEntry("V rel error",fp2str(verr));
            endLogGroup();

            if (verr < err) {
                addLogEntry("P iterations",int2str(n));
                addLogEntry("Final P",fp2str(pressure()));
                addLogEntry("V rel error",fp2str(verr));
                goto done;
            }
            // find dV/dP
            setPressure(pnow*1.01);
            dVdP = (volume() - vnow)/(0.01*pnow);
            setPressure(pnow + 0.5*(v0 - vnow)/dVdP);
        }
    }

    else {
        if (loglevel > 0) {
            endLogGroup();
        }
        throw CanteraError("MultiPhase::equilibrate","unknown option");
    }
    return -1.0;
done:
    delete e;
    e = 0;
    if (loglevel > 0) {
        endLogGroup();
    }
    return err;
}

#ifdef MULTIPHASE_DEVEL
void importFromXML(string infile, string id)
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* x = get_XML_Node(string("#")+id, root);
    if (x.name() != "multiphase")
        throw CanteraError("MultiPhase::importFromXML",
                           "Current XML_Node is not a multiphase element.");
    vector<XML_Node*> phases;
    x.getChildren("phase",phases);
    int np = phases.size();
    int n;
    ThermoPhase* p;
    for (n = 0; n < np; n++) {
        XML_Node& ph = *phases[n];
        srcfile = infile;
        if (ph.hasAttrib("src")) {
            srcfile = ph["src"];
        }
        idstr =  ph["id"];
        p = newPhase(srcfile, idstr);
        if (p) {
            addPhase(p, ph.value());
        }
    }
}
#endif

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
        throw IndexError("checkElementIndex", "elements", m, m_nel-1);
    }
}

void MultiPhase::checkElementArraySize(size_t mm) const
{
    if (m_nel > mm) {
        throw ArraySizeError("checkElementArraySize", mm, m_nel);
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
        throw IndexError("checkSpeciesIndex", "species", k, m_nsp-1);
    }
}

void MultiPhase::checkSpeciesArraySize(size_t kk) const
{
    if (m_nsp > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_nsp);
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
    std::string tmp;
    for (int iph = 0; iph < (int) m_np; iph++) {
        const ThermoPhase* tptr = m_phase[iph];
        tmp = tptr->id();
        if (tmp == pName) {
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
    size_t ip, loc = 0;
    for (ip = 0; ip < m_np; ip++) {
        ThermoPhase* p = m_phase[ip];
        p->getMoleFractions(DATA_PTR(m_moleFractions) + loc);
        loc += p->nSpecies();
    }
    calcElemAbundances();
}

void MultiPhase::updatePhases() const
{
    size_t p, nsp, loc = 0;
    for (p = 0; p < m_np; p++) {
        nsp = m_phase[p]->nSpecies();
        const doublereal* x = DATA_PTR(m_moleFractions) + loc;
        loc += nsp;
        m_phase[p]->setState_TPX(m_temp, m_press, x);
        m_temp_OK[p] = true;
        if (m_temp < m_phase[p]->minTemp()
                || m_temp > m_phase[p]->maxTemp()) {
            m_temp_OK[p] = false;
        }
    }
}
}
