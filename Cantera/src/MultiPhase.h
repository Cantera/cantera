#ifndef CT_MULTIPHASE_H
#define CT_MULTIPHASE_H

#include "ThermoPhase.h"
#include "DenseMatrix.h"
#include "stringUtils.h"
#include <iostream>


namespace Cantera {

    /// A class for multiphase mixtures. The mixture can contain any
    /// number of phases of any type.  All phases have the same
    /// temperature and pressure, and a specified number of moles.
    /// The phases do not need to have the same elements. For example,
    /// a mixture might consist of a gaseous phase with elements (H,
    /// C, O, N), a solid carbon phase containing only element C,
    /// etc. A master element set will be constructed for the mixture
    /// that is the union of the elements of each phase.

    class MultiPhase {

    public:

        typedef size_t       index_t;
        typedef ThermoPhase  phase_t;
        typedef DenseMatrix  array_t;

        /// Constructor. The constructor takes no arguments, since
        /// phases are added using method addPhase.
        MultiPhase() : m_temp(0.0), m_press(0.0), 
                       m_nel(0), m_nsp(0), m_init(false) {}

        /// Destructor. Does nothing. Class MultiPhase does not take 
        /// "ownership" (i.e. responsibility for destroying) the
        /// phase objects.  
        virtual ~MultiPhase() {}

        /// Add a phase to the mixture. 
        /// @param p pointer to the phase object
        /// @param moles total number of moles of all species in this phase
        void addPhase(phase_t* p, doublereal moles) {

            if (m_init) {
                throw CanteraError("addPhase",
                    "phases cannot be added after init() has been called.");
            }

            // save the pointer to the phase object
            m_phase.push_back(p);

            // store its number of moles
            m_moles.push_back(moles);

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
            //init();
        }

        int nElements() { return int(m_nel); }
        string elementName(int m) { return m_enames[m]; }
        int elementIndex(string name) { return m_enamemap[name] - 1;}
        
        int nSpecies() { return int(m_nsp); }
        string speciesName(int k) { return m_snames[k]; }
        doublereal nAtoms(int k, int m) {
            if (!m_init) init();
            return m_atoms(m,k); 
        }

        /// Species mole fractions. Write the array of species mole
        /// fractions into array \c x. The mole fractions are
        /// normalized to sum to one in each phase.
        void getMoleFractions(doublereal* x) {
            copy(m_moleFractions.begin(), m_moleFractions.end(), x);
        }


        /// Process phases and build atomic composition array. After 
        /// init() has been called, no more phases may be added.
        void init() {
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

        /// Moles of phase n.
        doublereal phaseMoles(index_t n) {
            return m_moles[n];
        }

        /// Set the number of moles of phase with index n.
        void setPhaseMoles(index_t n, doublereal moles) {
            m_moles[n] = moles;
        }

        /// Return a reference to phase n. The state of phase n is
        /// also updated to match the state stored locally in the 
        /// mixture object.
        phase_t& phase(index_t n) {
            if (!m_init) init();
            m_phase[n]->setState_TPX(m_temp, m_press, 
                m_moleFractions.begin() + m_spstart[n]);
            return *m_phase[n];
        }

        /// Return a const reference to phase n.
        //const phase_t& phase(index_t n) const {
            //            if (!m_init) init();
        //  m_phase[n]->setState_TPX(m_temp, 
        //      m_press, m_moleFractions.begin() + m_spstart[n]);
        //  return *m_phase[n];
        //}

        /// Moles of species \c k.
        doublereal speciesMoles(index_t k) {
            if (!m_init) init();
            index_t ip = m_spphase[k];
            return m_moles[ip]*m_moleFractions[k];
        }

        /// Index of the species belonging to phase number \c p
        /// with index \c k within the phase.
        int speciesIndex(index_t k, index_t p) {
            return m_spstart[p] + k;
        }

        /// Total moles of element m, summed over all
        /// phases
        doublereal elementMoles(index_t m) {
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
        void getChemPotentials(doublereal* mu) {
            index_t i, loc = 0;
            updatePhases();            
            for (i = 0; i < m_np; i++) {
                m_phase[i]->getChemPotentials(mu + loc);
                loc += m_phase[i]->nSpecies();
            }
        }

        /// Chemical potentials. Write into array \c mu the chemical
        /// potentials of all species [J/kmol].
        void getStandardChemPotentials(doublereal* mu) {
            index_t i, loc = 0;
            updatePhases();
            for (i = 0; i < m_np; i++) {
                m_phase[i]->getStandardChemPotentials(mu + loc);
                loc += m_phase[i]->nSpecies();
            }
        }

        /// Temperature [K].
        doublereal temperature() {
            return m_temp;
        }

        /// Set the temperature [K].
        void setTemperature(doublereal T) {
            m_temp = T;
            updatePhases();
        }

        doublereal pressure() {
            return m_press;
        }

        void setPressure(doublereal P) {
            m_press = P;
            updatePhases();
        }

        doublereal gibbs() {
            index_t i;
            doublereal sum = 0.0;
            updatePhases();
            for (i = 0; i < m_np; i++) 
                sum += m_phase[i]->gibbs_mole() * m_moles[i];
            return sum;
        }

        index_t nPhases() {
            return m_np;
        }

        bool solutionSpecies(index_t k) {
            if (m_phase[m_spphase[k]]->nSpecies() > 1)
                return true;
            else
                return false;
        }

        index_t speciesPhaseIndex(index_t k) {
            return m_spphase[k];
        }

        doublereal moleFraction(index_t k) {
            return m_moleFractions[k];
        }

        void updateMoleFractions() {
            if (!m_init) init();
            // save the current mole fractions for each phase
            index_t ip, loc = 0;
            for (ip = 0; ip < m_np; ip++) {
                phase_t* p = m_phase[ip];
                p->getMoleFractions(m_moleFractions.begin() + loc);
                loc += p->nSpecies();
            }
        }

        void setPhaseMoleFractions(index_t n, doublereal* x) {
            phase_t* p = m_phase[n];
            p->setState_TPX(m_temp, m_press, x);
        }

        void setMolesByName(compositionMap& xMap) {
            int kk = nSpecies();
            doublereal x;
            vector_fp mf(kk, 0.0);
            for (int k = 0; k < kk; k++) {
                x = xMap[speciesName(k)];
                if (x > 0.0) mf[k] = x;
            }
            setMoles(mf.begin());
        }

        void setMolesByName(const string& x) {
            compositionMap xx;
            int kk = nSpecies();
            for (int k = 0; k < kk; k++) { 
                xx[speciesName(k)] = -1.0;
            }
            parseCompString(x, xx);
            setMolesByName(xx); 
        }

        void setMoles(doublereal* n) {
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

    protected:

        /// Set the states of the phase objects to the locally-stored
        /// state.  Note that if individual phases have T and P different
        /// than that stored locally, the phase T and P will be modified.
        void updatePhases() {
            if (!m_init) init();
            index_t p, nsp, loc = 0;
            for (p = 0; p < m_np; p++) {
                nsp = m_phase[p]->nSpecies();
                doublereal* x = m_moleFractions.begin() + loc;
                loc += nsp;
                m_phase[p]->setState_TPX(m_temp, m_press, x);
            }
        }
            
        vector_fp m_moles;
        vector<phase_t*> m_phase;
        array_t m_atoms;
        vector_fp m_moleFractions;
        vector_int m_spphase;
        vector_int m_spstart;
        vector<string> m_enames;
        vector<string> m_snames;
        map<string, int> m_enamemap;
        index_t  m_np;
        doublereal m_temp;
        doublereal m_press;
        index_t m_nel; 
        index_t m_nsp;
        bool m_init;
    };

    inline std::ostream& operator<<(std::ostream& s, Cantera::MultiPhase& x) {
        size_t ip;
        for (ip = 0; ip < x.nPhases(); ip++) {
            s << "*************** Phase " << ip << " *****************" << endl;
            s << "Moles: " << x.phaseMoles(ip) << endl;
                
            s << report(x.phase(ip)) << endl;
        }
        return s;
    }
}

#endif
