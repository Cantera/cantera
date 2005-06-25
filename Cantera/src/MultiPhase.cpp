#include "MultiPhase.h"
#include "MultiPhaseEquil.h"

#include "ThermoPhase.h"
#include "DenseMatrix.h"
#include "stringUtils.h"

namespace Cantera {


     MultiPhase::MultiPhase() : m_temp(0.0), m_press(0.0), 
                    m_nel(0), m_nsp(0), m_init(false), m_eloc(-1), 
                    m_equil(0), m_Tmin(1.0), m_Tmax(100000.0) {
     }



    void MultiPhase::
    addPhase(phase_t* p, doublereal moles) {

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
        index_t m, nel = p->nElements();
        for (m = 0; m < nel; m++) {
            ename = p->elementName(m);

            // if no entry is found for this element name, then
            // it is a new element. In this case, add the name
            // to the list of names, increment the element count, 
            // and add an entry to the name->(index+1) map.
            if (m_enamemap[ename] == 0) {
                m_enamemap[ename] = m_nel + 1;
                m_enames.push_back(ename);
                m_atomicNumber.push_back(p->atomicNumber(m));
                if (ename == "E" || ename == "e") m_eloc = m_nel;
                m_nel++;
            }
        }
            
        if (m_temp == 0.0 && p->temperature() > 0.0) {
            m_temp = p->temperature();
            m_press = p->pressure();
        }
        //cout << "min, max = " << m_Tmin << " " << m_Tmax << endl;
        if (p->nSpecies() > 1) {
            double t = p->minTemp();
            if (t > m_Tmin) m_Tmin = t;
            t = p->maxTemp();
            if (t < m_Tmax) m_Tmax = t;
            //cout << p->name() << " " << t << " " << m_Tmax << endl;
        }

    }


    /// Process phases and build atomic composition array. After 
    /// init() has been called, no more phases may be added.
    void MultiPhase::init() {
        if (m_init) return;
        index_t ip, kp, k = 0, nsp, m;
        int mlocal;
        string sym;

        // allocate space for the atomic composition matrix
        m_atoms.resize(m_nel, m_nsp, 0.0);
        m_moleFractions.resize(m_nsp, 0.0);

        // iterate over the elements
        for (m = 0; m < m_nel; m++) {
            sym = m_enames[m];
            k = 0;
            // iterate over the phases
            for (ip = 0; ip < m_np; ip++) {
                phase_t* p = m_phase[ip];
                nsp = p->nSpecies();
                mlocal = p->elementIndex(sym);    
                for (kp = 0; kp < nsp; kp++) {
                    if (mlocal >= 0) {
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

        if (m_eloc >= 0) {
            doublereal esum;
            for (k = 0; k < m_nsp; k++) {
                esum = 0.0;
                for (m = 0; m < m_nel; m++) {
                    if (int(m) != m_eloc)
                        esum += m_atoms(m,k) * m_atomicNumber[m];
                }
                //m_atoms(m_eloc, k) += esum;
            }
        }

        /// set the initial composition within each phase to the
        /// mole fractions stored in the phase objects
        m_init = true;
        updateMoleFractions();
    }


    /// Return a reference to phase n. The state of phase n is
    /// also updated to match the state stored locally in the 
    /// mixture object.
    MultiPhase::phase_t& MultiPhase::phase(index_t n) {
        if (!m_init) init();
        m_phase[n]->setState_TPX(m_temp, m_press, 
            m_moleFractions.begin() + m_spstart[n]);
        return *m_phase[n];
    }

    /// Moles of species \c k.
    doublereal MultiPhase::speciesMoles(index_t k) {
        if (!m_init) init();
        index_t ip = m_spphase[k];
        return m_moles[ip]*m_moleFractions[k];
    }

    /// Total moles of element m, summed over all
    /// phases
    doublereal MultiPhase::elementMoles(index_t m) {
        doublereal sum = 0.0, phasesum;
        index_t i, k = 0, ik, nsp;
        if (!m_init) init();
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

    /// Total charge, summed over all phases
    doublereal MultiPhase::charge() {
        doublereal sum = 0.0;
        index_t i;
        for (i = 0; i < m_np; i++) {
            sum += phaseCharge(i);
        }
        return sum;
    }

    /// Charge of one phase
    doublereal MultiPhase::phaseCharge(index_t p) {
        doublereal phasesum = 0.0;
        int ik, k, nsp = m_phase[p]->nSpecies();
        for (ik = 0; ik < nsp; ik++) {
            k = speciesIndex(ik, p);
            phasesum += m_phase[p]->charge(ik)*m_moleFractions[k];
        }
        return Faraday*phasesum*m_moles[p];
    }


    /// Chemical potentials. Write into array \c mu the chemical
    /// potentials of all species [J/kmol].
    void MultiPhase::getChemPotentials(doublereal* mu) {
        index_t i, loc = 0;
        updatePhases();            
        for (i = 0; i < m_np; i++) {
            m_phase[i]->getChemPotentials(mu + loc);
            loc += m_phase[i]->nSpecies();
        }
    }

    /// Chemical potentials. Write into array \c mu the chemical
    /// potentials of all species [J/kmol].
    void MultiPhase::getValidChemPotentials(doublereal not_mu,
        doublereal* mu, bool standard) {
        index_t i, loc = 0;
        updatePhases();            
        for (i = 0; i < m_np; i++) {
            if (tempOK(i) || m_phase[i]->nSpecies() > 1) {
                if (!standard)
                    m_phase[i]->getChemPotentials(mu + loc);
                else
                    m_phase[i]->getStandardChemPotentials(mu + loc);
            }
            else
                fill(mu + loc, mu + loc + m_phase[i]->nSpecies(), not_mu);
            loc += m_phase[i]->nSpecies();
        }
    }


    /// Chemical potentials. Write into array \c mu the chemical
    /// potentials of all species [J/kmol].
    void MultiPhase::getStandardChemPotentials(doublereal* mu) {
        index_t i, loc = 0;
        updatePhases();
        for (i = 0; i < m_np; i++) {
            m_phase[i]->getStandardChemPotentials(mu + loc);
            loc += m_phase[i]->nSpecies();
        }
    }

    bool MultiPhase::solutionSpecies(index_t k) {
        if (m_phase[m_spphase[k]]->nSpecies() > 1)
            return true;
        else
            return false;
    }

    doublereal MultiPhase::gibbs() {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->gibbs_mole() * m_moles[i];
        return sum;
    }

    doublereal MultiPhase::enthalpy() {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->enthalpy_mole() * m_moles[i];
        return sum;
    }

    doublereal MultiPhase::entropy() {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->entropy_mole() * m_moles[i];
        return sum;
    }

    doublereal MultiPhase::cp() {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->cp_mole() * m_moles[i];
        return sum;
    }

    void MultiPhase::updateMoleFractions() {
        if (!m_init) init();
        // save the current mole fractions for each phase
        index_t ip, loc = 0;
        for (ip = 0; ip < m_np; ip++) {
            phase_t* p = m_phase[ip];
            p->getMoleFractions(m_moleFractions.begin() + loc);
            loc += p->nSpecies();
        }
    }

    void MultiPhase::setPhaseMoleFractions(index_t n, doublereal* x) {
        phase_t* p = m_phase[n];
        p->setState_TPX(m_temp, m_press, x);
    }

    void MultiPhase::setMolesByName(compositionMap& xMap) {
        if (!m_init) init();
        int kk = nSpecies();
        doublereal x;
        vector_fp mf(kk, 0.0);
        for (int k = 0; k < kk; k++) {
            x = xMap[speciesName(k)];
            if (x > 0.0) mf[k] = x;
        }
        setMoles(mf.begin());
    }

    void MultiPhase::setMolesByName(const string& x) {
        compositionMap xx;
        if (!m_init) init();
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) { 
            xx[speciesName(k)] = -1.0;
        }
        parseCompString(x, xx);
        setMolesByName(xx); 
    }

    void MultiPhase::setMoles(doublereal* n) {
        if (!m_init) init();
        index_t ip, loc = 0;
        index_t ik, k = 0, nsp;
        doublereal phasemoles;
        for (ip = 0; ip < m_np; ip++) {
            phase_t* p = m_phase[ip];
            nsp = p->nSpecies();
            phasemoles = 0.0;
            for (ik = 0; ik < nsp; ik++) {
                phasemoles += n[k];
                k++;
            }
            m_moles[ip] = phasemoles;
            if (nsp > 1) {
                p->setState_TPX(m_temp, m_press, n + loc);
                p->getMoleFractions(m_moleFractions.begin() + loc);
            }
            else {
                m_moleFractions[loc] = 1.0;
            }
            loc += p->nSpecies();
        }
    }

    doublereal MultiPhase::volume() {
        int i;
        doublereal sum = 0;
        for (i = 0; i < int(m_np); i++) {
            sum += m_moles[i]/m_phase[i]->molarDensity();
        }
        return sum;
    }

    void MultiPhase::updatePhases() {
        if (!m_init) init();
        index_t p, nsp, loc = 0;
        for (p = 0; p < m_np; p++) {
            nsp = m_phase[p]->nSpecies();
            doublereal* x = m_moleFractions.begin() + loc;
            loc += nsp;
            m_phase[p]->setState_TPX(m_temp, m_press, x);
            m_temp_OK[p] = true;
            if (m_temp < m_phase[p]->minTemp() 
                || m_temp > m_phase[p]->maxTemp()) m_temp_OK[p] = false;
        }
    }            

    doublereal MultiPhase::equilibrate(int XY, doublereal err, 
        int maxsteps, int maxiter, int loglevel) {
        doublereal error;
        bool strt = false;
        doublereal dt;
        doublereal h0;
        int n;
        bool start, once;
        doublereal ferr, hnow, herr = 1.0;
        doublereal Tlow = -1.0, Thigh = -1.0;
        doublereal hlow = 0.0, hhigh = 0.0, slope, tnew;
        doublereal dta, dtmax;

        if (!m_init) init();
        if (loglevel > 0) {
            beginLogGroup("MultiPhase::equilibrate");
        }
        if (XY == TP) {
            if (loglevel > 0) {
                addLogEntry("problem type","fixed T,P");
            }
            // create an equilibrium manager 
            MultiPhaseEquil e(this);
            error = e.equilibrate(XY, err, maxsteps, loglevel-1);
            if (loglevel > 0)  e.printInfo();
            goto done;
        }
        else if (XY == HP) {
            dt = 1.0e2;
            h0 = enthalpy();
            start = true;
            Tlow = m_Tmin;      // lower bound on T
            Thigh = m_Tmax;   // upper bound on T
            hlow = 0.0;
            hhigh = 0.0;
            once = true;
            if (loglevel > 0) {
                addLogEntry("problem type","fixed H,P");
                addLogEntry("H target",fp2str(h0));
                addLogEntry("min T",fp2str(Tlow));
                addLogEntry("max T",fp2str(Thigh));
            }
            ferr = 0.1;
            for (n = 0; n < maxiter; n++) {
                MultiPhaseEquil e(this, strt);
                start = false;
                if (loglevel > 0) {
                    beginLogGroup("iteration "+int2str(n));
                }
                try {
                    error = e.equilibrate(TP, err, maxsteps, loglevel-1);

                    hnow = enthalpy();
                    if (hnow < h0) {
                        if (m_temp > Tlow) {
                            Tlow = m_temp;
                        }
                    }
                    else {
                        if (m_temp < Thigh) {
                            Thigh = m_temp;
                        }
                    }
                    herr = fabs((h0 - hnow)/h0);
                    if (loglevel > 0) {
                        addLogEntry("T",fp2str(temperature()));
                        addLogEntry("H",fp2str(hnow));
                        addLogEntry("H rel error",fp2str(herr));
                        endLogGroup();
                    }
                    dt = (h0 - hnow)/cp();
                    dtmax = 0.5*(Thigh - Tlow);
                    dta = fabs(dt);
                    if (dta > dtmax) dt *= dtmax/dta;
                    if (herr < err || dta < 1.0e-4) {
                        if (loglevel > 0) {
                            addLogEntry("T iterations",int2str(n));
                            addLogEntry("Final T",fp2str(temperature()));
                            addLogEntry("H rel error",fp2str(herr));
                        }
                        goto done;
                    }
                    tnew = m_temp + dt;
                    setTemperature(tnew);
                    if (dta < 100.0) strt = false;
                }
                catch (CanteraError e) {
                    if (!strt) {
                        if (loglevel > 0) 
                            addLogEntry("no convergence","setting strt to True");
                        strt = true;
                    }
                    else {
                        tnew = 0.5*(m_temp + Thigh);
                        setTemperature(tnew);
                        if (loglevel > 0) 
                            addLogEntry("no convergence",
                                "trying T = "+fp2str(m_temp));
                            
                    }
                }
            }
            throw CanteraError("MultiPhase::equilibrate",
                "No convergence for T");
        }
        else if (XY == SP) {
            if (loglevel > 0) {
                addLogEntry("problem type","fixed S,P");
            }
            doublereal dt = 1.0e3;
            doublereal s0 = entropy();
            int n;
            bool start = true;
            doublereal ferr, snow, serr, tnew;
            for (n = 0; n < maxiter; n++) {
                MultiPhaseEquil e(this, start);
                ferr = 0.1;
                start = false;
                if (fabs(dt) < 1.0) ferr = err;
                if (loglevel > 1) {
                    beginLogGroup("iteration "+int2str(n));
                }
                error = e.equilibrate(TP, ferr, maxsteps, loglevel-1);
                snow = entropy();
                tnew = exp(0.5*(s0 - snow)/cp())*temperature();
                serr = fabs((s0 - snow)/s0);
                if (loglevel > 1) {
                    addLogEntry("T",fp2str(temperature()));
                    addLogEntry("S rel error",fp2str(serr));
                    endLogGroup();
                }
                if (serr < err) {
                    if (loglevel > 0) {
                        addLogEntry("T iterations",int2str(n));
                        addLogEntry("Final T",fp2str(temperature()));
                        addLogEntry("S rel error",fp2str(serr));
                    }
                    goto done;
                }
                setTemperature(tnew);
            }
        }
        else if (XY == TV) {
            if (loglevel > 0) {
                addLogEntry("problem type","fixed T, V");
            }
            doublereal dt = 1.0e3;
            doublereal v0 = volume();
            doublereal dVdP;
            int n;
            bool start = true;
            doublereal error, ferr, vnow, pnow, verr, tnew;
            for (n = 0; n < maxiter; n++) {
                pnow = pressure();
                MultiPhaseEquil e(this, start);
                start = false;
                if (loglevel > 1) {
                    beginLogGroup("iteration "+int2str(n));
                }
                error = e.equilibrate(TP, err, maxsteps, loglevel-1);
                vnow = volume();
                verr = fabs((v0 - vnow)/v0);
                if (loglevel > 1) {
                    addLogEntry("P",fp2str(pressure()));
                    addLogEntry("V rel error",fp2str(verr));
                    endLogGroup();
                }
                if (verr < err) {
                    if (loglevel > 0) {
                        addLogEntry("P iterations",int2str(n));
                        addLogEntry("Final P",fp2str(pressure()));
                        addLogEntry("V rel error",fp2str(verr));
                    }
                    goto done;
                }
                // find dV/dP
                setPressure(pnow*1.01);
                dVdP = (volume() - vnow)/(0.01*pnow);
                setPressure(pnow + 0.5*(v0 - vnow)/dVdP);
            }
        }

        else {
            if (loglevel > 0)  endLogGroup();
            throw CanteraError("MultiPhase::equilibrate","unknown option");
        }
        return -1.0;
done:
        if (loglevel > 0)  {
            endLogGroup();
        }
        return err;
    }
}

