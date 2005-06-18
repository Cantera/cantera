#ifndef CT_MULTIPHASE_H
#define CT_MULTIPHASE_H

#include "ct_defs.h"
#include "DenseMatrix.h"
#include "ThermoPhase.h"

namespace Cantera {

    class MultiPhaseEquil;

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
                       m_nel(0), m_nsp(0), m_init(false), m_eloc(-1), 
                       m_equil(0) {}

        /// Destructor. Does nothing. Class MultiPhase does not take 
        /// "ownership" (i.e. responsibility for destroying) the
        /// phase objects.  
        virtual ~MultiPhase() {}

        /// Add a phase to the mixture. 
        /// @param p pointer to the phase object
        /// @param moles total number of moles of all species in this phase
        void addPhase(phase_t* p, doublereal moles);

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
        void init();

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
        phase_t& phase(index_t n);

        /// Moles of species \c k.
        doublereal speciesMoles(index_t k);

        /// Index of the species belonging to phase number \c p
        /// with index \c k within the phase.
        int speciesIndex(index_t k, index_t p) {
            return m_spstart[p] + k;
        }


        doublereal charge();
        doublereal phaseCharge(index_t p);

        /// Total moles of element m, summed over all
        /// phases
        doublereal elementMoles(index_t m);

        /// Chemical potentials. Write into array \c mu the chemical
        /// potentials of all species [J/kmol].
        void getChemPotentials(doublereal* mu);

        /// Valid chemical potentials. Write into array \c mu the
        /// chemical potentials of all species with thermo data valid
        /// for the current temperature [J/kmol]. For other species,
        /// set the chemical potential to the value \c not_mu.
        void getValidChemPotentials(doublereal not_mu, doublereal* mu);

        /// Chemical potentials. Write into array \c mu the chemical
        /// potentials of all species [J/kmol].
        void getStandardChemPotentials(doublereal* mu);

        /// Temperature [K].
        doublereal temperature() {
            return m_temp;
        }

        doublereal equilibrate(int XY, doublereal err = 1.0e-9, 
            int maxsteps = 1000, int maxiter = 200, int loglevel = 0);


        /// Set the temperature [K].
        void setTemperature(doublereal T) {
            m_temp = T;
            updatePhases();
        }

        doublereal pressure() {
            return m_press;
        }

        doublereal volume();

        void setPressure(doublereal P) {
            m_press = P;
            updatePhases();
        }

        doublereal enthalpy();
        doublereal entropy();
        doublereal gibbs();
        doublereal cp();

        index_t nPhases() {
            return m_np;
        }

        bool solutionSpecies(index_t k);

        index_t speciesPhaseIndex(index_t k) {
            return m_spphase[k];
        }

        doublereal moleFraction(index_t k) {
            return m_moleFractions[k];
        }

        void updateMoleFractions();

        void setPhaseMoleFractions(index_t n, doublereal* x);

        void setMolesByName(compositionMap& xMap);

        void setMolesByName(const string& x);

        void setMoles(doublereal* n);

        bool tempOK(index_t p) {
            return m_temp_OK[p];
        }

    protected:

        /// Set the states of the phase objects to the locally-stored
        /// state.  Note that if individual phases have T and P different
        /// than that stored locally, the phase T and P will be modified.
        void updatePhases();
            
        vector_fp m_moles;
        vector<phase_t*> m_phase;
        array_t m_atoms;
        vector_fp m_moleFractions;
        vector_int m_spphase;
        vector_int m_spstart;
        vector<string> m_enames;
        vector_int m_atomicNumber;
        vector<string> m_snames;
        map<string, int> m_enamemap;
        index_t  m_np;
        doublereal m_temp;
        doublereal m_press;
        index_t m_nel; 
        index_t m_nsp;
        bool m_init;
        int m_eloc;
        vector<bool> m_temp_OK;
        MultiPhaseEquil* m_equil;
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
