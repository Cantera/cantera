#include "MultiPhase.h"
#include "MultiPhaseEquil.h"

#include "ThermoPhase.h"
#include "DenseMatrix.h"
#include "stringUtils.h"

namespace Cantera {

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
                m_nel++;
            }
        }
            
        if (m_temp == 0.0 && p->temperature() > 0.0) {
            m_temp = p->temperature();
            m_press = p->pressure();
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
        doublereal* mu) {
        index_t i, loc = 0;
        updatePhases();            
        for (i = 0; i < m_np; i++) {
            if (tempOK(i) || m_phase[i]->nSpecies() > 1) 
                m_phase[i]->getChemPotentials(mu + loc);
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
        int maxsteps) {
        cout << "in equil" << endl;
        init();
        if (m_equil == 0) {
            m_equil = new MultiPhaseEquil(this);
        }
        m_equil->equilibrate(XY, err, maxsteps);
        delete m_equil;
        m_equil = 0;
    }

}

