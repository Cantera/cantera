#include "MultiPhase.h"
#include "MultiPhaseEquil.h"

#include "ThermoPhase.h"
#include "DenseMatrix.h"
#include "stringUtils.h"

namespace Cantera {

    /// Constructor.
     MultiPhase::MultiPhase() : m_temp(0.0), m_press(0.0), 
                    m_nel(0), m_nsp(0), m_init(false), m_eloc(-1), 
                    m_Tmin(1.0), m_Tmax(100000.0) {
     }

    void MultiPhase::
    addPhases(MultiPhase& mix) {
        index_t n;
        for (n = 0; n < mix.m_np; n++) {
            addPhase(mix.m_phase[n], mix.m_moles[n]);
        }
    }

    void MultiPhase::
    addPhases(phase_list& phases, const vector_fp& phaseMoles) {
        index_t np = phases.size();
        index_t n;
        for (n = 0; n < np; n++) {
            addPhase(phases[n], phaseMoles[n]);
        }
        init();
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
            if (m_enamemap.find(ename) == m_enamemap.end()) {
                m_enamemap[ename] = m_nel + 1;
                m_enames.push_back(ename);
                m_atomicNumber.push_back(p->atomicNumber(m));

                // Element 'E' (or 'e') is special. Note its location.
                if (ename == "E" || ename == "e") m_eloc = m_nel;

                m_nel++;
            }
        }

        // If the mixture temperature hasn't been set, then set the
        // temperature and pressure to the values for the phase being
        // added.
        if (m_temp == 0.0 && p->temperature() > 0.0) {
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
            if (t > m_Tmin) m_Tmin = t;
            t = p->maxTemp();
            if (t < m_Tmax) m_Tmax = t;
        }
    }


    /// Process phases and build atomic composition array. This method
    /// must be called after all phases are added, before doing
    /// anything else with the mixture. After init() has been called,
    /// no more phases may be added.
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
            DATA_PTR(m_moleFractions) + m_spstart[n]);
        return *m_phase[n];
    }

    /// Moles of species \c k.
    doublereal MultiPhase::speciesMoles(index_t k) const {
        index_t ip = m_spphase[k];
        return m_moles[ip]*m_moleFractions[k];
    }

    /// Total moles of element m, summed over all
    /// phases
    doublereal MultiPhase::elementMoles(index_t m) const {
        doublereal sum = 0.0, phasesum;
        index_t i, k = 0, ik, nsp;
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
    doublereal MultiPhase::charge() const {
        doublereal sum = 0.0;
        index_t i;
        for (i = 0; i < m_np; i++) {
            sum += phaseCharge(i);
        }
        return sum;
    }

    /// Net charge of one phase (Coulombs). The net charge is computed as
    /// \f[ Q_p = N_p \sum_k F z_k X_k \f]
    /// where the sum runs only over species in phase \a p.
    /// @param p index of the phase for which the charge is desired.   
    doublereal MultiPhase::phaseCharge(index_t p) const {
        doublereal phasesum = 0.0;
        int ik, k, nsp = m_phase[p]->nSpecies();
        for (ik = 0; ik < nsp; ik++) {
            k = speciesIndex(ik, p);
            phasesum += m_phase[p]->charge(ik)*m_moleFractions[k];
        }
        return Faraday*phasesum*m_moles[p];
    }


    /// Get the chemical potentials of all species in all phases. 
    void MultiPhase::getChemPotentials(doublereal* mu) const {
        index_t i, loc = 0;
        updatePhases();            
        for (i = 0; i < m_np; i++) {
            m_phase[i]->getChemPotentials(mu + loc);
            loc += m_phase[i]->nSpecies();
        }
    }

    /// Get chemical potentials of species with valid thermo
    /// data. This method is designed for use in computing chemical
    /// equilibrium by Gibbs minimization. For solution phases (more
    /// than one species), this does the same thing as
    /// getChemPotentials. But for stoichiometric phases, this writes
    /// into array \a mu the user-specified value \a not_mu instead of
    /// the chemical potential if the temperature is outside the range
    /// for which the thermo data for the one species in the phase are
    /// valid. The need for this arises since many condensed phases
    /// have thermo data fit only for the temperature range for which
    /// they are stable. For example, in the NASA database, the fits
    /// for H2O(s) are only done up to 0 C, the fits for H2O(L) are
    /// only done from 0 C to 100 C, etc. Using the polynomial fits outside
    /// the range for which the fits were done can result in spurious
    /// chemical potentials, and can lead to condensed phases
    /// appearing when in fact they should be absent.
    ///
    /// By setting \a not_mu to a large positive value, it is possible
    /// to force routines which seek to minimize the Gibbs free energy
    /// of the mixture to zero out any phases outside the temperature
    /// range for which their thermo data are valid.
    ///
    /// If this method is called with \a standard set to true, then
    /// the composition-independent standard chemical potentials are
    /// returned instead of the composition-dependent chemical
    /// potentials.
    ///
    void MultiPhase::getValidChemPotentials(doublereal not_mu,
        doublereal* mu, bool standard) const {
        index_t i, loc = 0;

        updatePhases();           
        // iterate over the phases 
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

    /// True if species \a k belongs to a solution phase.
    bool MultiPhase::solutionSpecies(index_t k) const {
        if (m_phase[m_spphase[k]]->nSpecies() > 1)
            return true;
        else
            return false;
    }

    /// The Gibbs free energy of the mixture (J).
    doublereal MultiPhase::gibbs() const {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->gibbs_mole() * m_moles[i];
        return sum;
    }

    /// The enthalpy of the mixture (J).
    doublereal MultiPhase::enthalpy() const {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->enthalpy_mole() * m_moles[i];
        return sum;
    }

    /// The entropy of the mixture (J/K).
    doublereal MultiPhase::entropy() const {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->entropy_mole() * m_moles[i];
        return sum;
    }

    /// The specific heat at constant pressure and composition (J/K).
    /// Note that this does not account for changes in composition of
    /// the mixture with temperature.
    doublereal MultiPhase::cp() const {
        index_t i;
        doublereal sum = 0.0;
        updatePhases();
        for (i = 0; i < m_np; i++) 
            sum += m_phase[i]->cp_mole() * m_moles[i];
        return sum;
    }



    /// Set the mole fractions of phase \a n to the values in 
    /// array \a x.  
    void MultiPhase::setPhaseMoleFractions(index_t n, doublereal* x) {
        phase_t* p = m_phase[n];
        p->setState_TPX(m_temp, m_press, x);
    }

    /// Set the species moles using a map. The map \a xMap maps
    /// species name strings to mole numbers. Mole numbers that are
    /// less than or equal to zero will be set to zero.
    void MultiPhase::setMolesByName(compositionMap& xMap) {
        int kk = nSpecies();
        doublereal x;
        vector_fp moles(kk, 0.0);
        for (int k = 0; k < kk; k++) {
            x = xMap[speciesName(k)];
            if (x > 0.0) moles[k] = x;
        }
        setMoles(DATA_PTR(moles));
    }


    /// Set the species moles using a string. Unspecified species are
    /// set to zero.
    void MultiPhase::setMolesByName(const string& x) {
        compositionMap xx;

        // add an entry in the map for every species, with value -1.0.
        // Function parseCompString (stringUtils.cpp) uses the names
        // in the map to specify the allowed species.
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) { 
            xx[speciesName(k)] = -1.0;
        }

        // build the composition map from the string, and then set the
        // moles.
        parseCompString(x, xx);
        setMolesByName(xx); 
    }

    /// Set the species moles to the values in array \a n. The state
    /// of each phase object is also updated to have the specified
    /// composition and the mixture temperature and pressure.
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
                p->getMoleFractions(DATA_PTR(m_moleFractions) + loc);
            }
            else {
                m_moleFractions[loc] = 1.0;
            }
            loc += p->nSpecies();
        }
    }

    /// The total mixture volume [m^3].
    doublereal MultiPhase::volume() const {
        int i;
        doublereal sum = 0;
        for (i = 0; i < int(m_np); i++) {
            sum += m_moles[i]/m_phase[i]->molarDensity();
        }
        return sum;
    }

    doublereal MultiPhase::equilibrate(int XY, doublereal err, 
        int maxsteps, int maxiter, int loglevel) {
        doublereal error;
        bool strt = false;
        doublereal dt;
        doublereal h0;
        int n;
        bool start;
        doublereal ferr, hnow, herr = 1.0;
        doublereal snow, serr = 1.0, s0;
        doublereal Tlow = -1.0, Thigh = -1.0;
        doublereal Hlow = Undef, Hhigh = Undef, tnew;
        doublereal dta=0.0, dtmax, cpb;
        MultiPhaseEquil* e = 0;

        if (!m_init) init();
        beginLogGroup("MultiPhase::equilibrate", loglevel);

        if (XY == TP) {
            addLogEntry("problem type","fixed T,P");
            addLogEntry("Temperature",temperature());
            addLogEntry("Pressure", pressure());


            // create an equilibrium manager 
            e = new MultiPhaseEquil(this);
            try {
                error = e->equilibrate(XY, err, maxsteps);
            }
            catch (CanteraError err) {
                endLogGroup();
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
            addLogEntry("problem type","fixed H,P");
            addLogEntry("H target",fp2str(h0));

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
                beginLogGroup("iteration "+int2str(n));

                try {
                    error = e->equilibrate(TP, err, maxsteps);
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
                        if (dta > dtmax) dt *= dtmax/dta;  
                    }
                    else {
                        tnew = sqrt(Tlow*Thigh);
                        dt = tnew - m_temp;
                        //cpb = cp();
                    }

                    herr = fabs((h0 - hnow)/h0);
                    addLogEntry("T",fp2str(temperature()));
                    addLogEntry("H",fp2str(hnow));
                    addLogEntry("H rel error",fp2str(herr));
                    addLogEntry("lower T bound",fp2str(Tlow));
                    addLogEntry("upper T bound",fp2str(Thigh));
                    endLogGroup(); // iteration


                    if (herr < err) { // || dta < 1.0e-4) {
                        addLogEntry("T iterations",int2str(n));
                        addLogEntry("Final T",fp2str(temperature()));
                        addLogEntry("H rel error",fp2str(herr));
                        goto done;
                    }
                    tnew = m_temp + dt;
                    if (tnew < 0.0) tnew = 0.5*m_temp;
                    //dta = fabs(tnew - m_temp);
                    setTemperature(tnew);

                    // if the size of Delta T is not too large, use
                    // the current composition as the starting estimate
                    if (dta < 100.0) strt = false;

                }

                catch (CanteraError err) {
                    if (!strt) {
                        addLogEntry("no convergence",
                            "try estimating starting composition");
                        strt = true;
                    }
                    else {
                        tnew = 0.5*(m_temp + Thigh);
                        if (fabs(tnew - m_temp) < 1.0) tnew = m_temp + 1.0;
                        setTemperature(tnew);
                        addLogEntry("no convergence",
                            "trying T = "+fp2str(m_temp));
                    }
                    endLogGroup();
                }
                delete e;
                e = 0;
            }
            addLogEntry("reached max number of T iterations",int2str(maxiter));
            endLogGroup();
            throw CanteraError("MultiPhase::equilibrate",
                "No convergence for T");
        }
        else if (XY == SP) {
            s0 = entropy();
            start = true;
            Tlow = 1.0; // m_Tmin;      // lower bound on T
            Thigh = 1.0e6; // m_Tmax;   // upper bound on T
            addLogEntry("problem type","fixed S,P");
            addLogEntry("S target",fp2str(s0));
            addLogEntry("min T",fp2str(Tlow));
            addLogEntry("max T",fp2str(Thigh));
            
            for (n = 0; n < maxiter; n++) {
                if (e) delete e;
                e = new MultiPhaseEquil(this, strt);
                ferr = 0.1;
                if (fabs(dt) < 1.0) ferr = err;
                //start = false;
                beginLogGroup("iteration "+int2str(n));
                
                try {
                    error = e->equilibrate(TP, err, maxsteps);
                    snow = entropy();
                    if (snow < s0) {
                        if (m_temp > Tlow) Tlow = m_temp;
                    }
                    else {
                        if (m_temp < Thigh) Thigh = m_temp;
                    }
                    serr = fabs((s0 - snow)/s0);
                    addLogEntry("T",fp2str(temperature()));
                    addLogEntry("S",fp2str(snow));
                    addLogEntry("S rel error",fp2str(serr));
                    endLogGroup();
                    
                    dt = (s0 - snow)*m_temp/cp();
                    dtmax = 0.5*fabs(Thigh - Tlow);
                    dtmax = (dtmax > 500.0 ? 500.0 : dtmax);
                    dta = fabs(dt);
                    if (dta > dtmax) dt *= dtmax/dta;
                    if (herr < err || dta < 1.0e-4) {
                        addLogEntry("T iterations",int2str(n));
                        addLogEntry("Final T",fp2str(temperature()));
                        addLogEntry("S rel error",fp2str(serr));
                        goto done;
                    }
                    tnew = m_temp + dt;
                    setTemperature(tnew);

                    // if the size of Delta T is not too large, use
                    // the current composition as the starting estimate
                    if (dta < 100.0) strt = false;
                }

                catch (CanteraError err) {
                    if (!strt) {
                        addLogEntry("no convergence",
                            "setting strt to True");
                        strt = true;
                    }
                    else {
                        tnew = 0.5*(m_temp + Thigh);
                        setTemperature(tnew);
                        addLogEntry("no convergence",
                            "trying T = "+fp2str(m_temp));
                            
                    }
                    endLogGroup();
                }
                delete e;
                e = 0;
            }
            addLogEntry("reached max number of T iterations",int2str(maxiter));
            endLogGroup();
            throw CanteraError("MultiPhase::equilibrate",
                "No convergence for T");
        }

//         else if (XY == SP) {
//             if (loglevel > 0) {
//                 addLogEntry("problem type","fixed S,P");
//             }
//             doublereal dt = 1.0e3;
//             doublereal s0 = entropy();
//             int n;
//             bool start = true;
//             doublereal ferr, snow, serr, tnew;
//             for (n = 0; n < maxiter; n++) {
//                 e = new MultiPhaseEquil(this, start);
//                 ferr = 0.1;
//                 start = false;
//                 if (fabs(dt) < 1.0) ferr = err;
//                 if (loglevel > 1) {
//                     beginLogGroup("iteration "+int2str(n));
//                 }
//                 try {
//                     error = e->equilibrate(TP, ferr, maxsteps, loglevel-1);
//                     snow = entropy();
//                     tnew = exp(0.5*(s0 - snow)/cp())*temperature();
//                     serr = fabs((s0 - snow)/s0);
//                     if (loglevel > 1) {
//                         addLogEntry("T",fp2str(temperature()));
//                         addLogEntry("S rel error",fp2str(serr));
//                         endLogGroup();
//                     }
//                     if (serr < err) {
//                         if (loglevel > 0) {
//                             addLogEntry("T iterations",int2str(n));
//                             addLogEntry("Final T",fp2str(temperature()));
//                             addLogEntry("S rel error",fp2str(serr));
//                         }
//                         goto done;
//                     }
//                     setTemperature(tnew);
//                 }
//                 catch (CanteraError err) {
//                     delete e;
//                     if (!strt) {
//                         if (loglevel > 0) 
//                             addLogEntry("no convergence",
//                                 "setting strt to True");
//                         strt = true;
//                     }
//                     else {
//                         tnew = 0.5*(m_temp + Thigh);
//                         setTemperature(tnew);
//                         if (loglevel > 0) 
//                             addLogEntry("no convergence",
//                                 "trying T = "+fp2str(m_temp));   
//                     }
//                 }
//                 endLogGroup();
//             }
//             if (loglevel > 0) write_logfile("equil_err.html");
//             throw CanteraError("MultiPhase::equilibrate",
//                 "No convergence for T");
//         }
        else if (XY == TV) {
            addLogEntry("problem type","fixed T, V");
            //            doublereal dt = 1.0e3;
            doublereal v0 = volume();
            doublereal dVdP;
            int n;
            bool start = true;
            doublereal error, vnow, pnow, verr;
            for (n = 0; n < maxiter; n++) {
                pnow = pressure();
                MultiPhaseEquil e(this, start);
                start = false;
                beginLogGroup("iteration "+int2str(n));
                
                error = e.equilibrate(TP, err, maxsteps);
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
            endLogGroup();
            throw CanteraError("MultiPhase::equilibrate","unknown option");
        }
        return -1.0;
done:
        delete e;
        e = 0;
        endLogGroup();
        return err;
    }

#ifdef MULTIPHASE_DEVEL
    void importFromXML(string infile, string id) {
        XML_Node* root = get_XML_File(infile);
        if (id == "-") id = "";
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
            if (ph.hasAttrib("src")) srcfile = ph["src"];
            idstr =  ph["id"];
            p = newPhase(srcfile, idstr);
            if (p) {
                addPhase(p, ph.value());
            }
        }
    }
#endif

    //-------------------------------------------------------------
    //
    // protected methods
    //
    //-------------------------------------------------------------


    /// Update the locally-stored species mole fractions. 
    void MultiPhase::updateMoleFractions() {
        index_t ip, loc = 0;
        for (ip = 0; ip < m_np; ip++) {
            phase_t* p = m_phase[ip];
            p->getMoleFractions(DATA_PTR(m_moleFractions) + loc);
            loc += p->nSpecies();
        }
    }


    /// synchronize the phase objects with the mixture state. This
    /// method sets each phase to the mixture temperature and
    /// pressure, and sets the phase mole fractions based on the
    /// mixture mole numbers. 
    void MultiPhase::updatePhases() const {
        index_t p, nsp, loc = 0;
        for (p = 0; p < m_np; p++) {
            nsp = m_phase[p]->nSpecies();
            const doublereal* x = DATA_PTR(m_moleFractions) + loc;
            loc += nsp;
            m_phase[p]->setState_TPX(m_temp, m_press, x);
            m_temp_OK[p] = true;
            if (m_temp < m_phase[p]->minTemp() 
                || m_temp > m_phase[p]->maxTemp()) m_temp_OK[p] = false;
        }
    }            

}

