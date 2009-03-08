/**
 * @file MultiPhase.h
 *
 *  $Author: hkmoffa $
 *  $Date: 2006/10/20 21:20:27 $
 *  $Revision: 1.15 $
 */
#ifndef CT_MULTIPHASE_H
#define CT_MULTIPHASE_H

#include "ct_defs.h"
#include "DenseMatrix.h"
#include "ThermoPhase.h"

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

        // some typedefs for convenience
        typedef size_t       index_t;
        typedef ThermoPhase  phase_t;
        typedef DenseMatrix  array_t;
        typedef vector<phase_t*> phase_list;

        /// Constructor. The constructor takes no arguments, since
        /// phases are added using method addPhase.
        MultiPhase();

        /// Destructor. Does nothing. Class MultiPhase does not take 
        /// "ownership" (i.e. responsibility for destroying) the
        /// phase objects.  
        virtual ~MultiPhase() {}

        
        void addPhases(phase_list& phases, const vector_fp& phaseMoles);

        /// Add all phases present in 'mix' to this mixture.
        void addPhases(MultiPhase& mix);

        /// Add a phase to the mixture. 
        /// @param p pointer to the phase object
        /// @param moles total number of moles of all species in this phase
        void addPhase(phase_t* p, doublereal moles);

        /// Number of elements.
        int nElements() const { return int(m_nel); }

        /// Name of element \a m.
        string elementName(int m) const { return m_enames[m]; }

        /// Index of element with name \a name.
        int elementIndex(string name) const { return m_enamemap[name] - 1;}
        
        /// Number of species, summed over all phases.
        int nSpecies() const { return int(m_nsp); }

        /// Name of species with index \a k. 
        string speciesName(int k) const { return m_snames[k]; }

        /// Number of atoms of element \a m in species \a k.
        doublereal nAtoms(int k, int m) {
            if (!m_init) init();
            return m_atoms(m,k); 
        }

        /// Species mole fractions. Write the array of species mole
        /// fractions into array \c x. The mole fractions are
        /// normalized to sum to one in each phase.
        void getMoleFractions(doublereal* x) const {
            copy(m_moleFractions.begin(), m_moleFractions.end(), x);
        }

        /// Process phases and build atomic composition array. After 
        /// init() has been called, no more phases may be added.
        void init();

        /// Moles of phase n.
        doublereal phaseMoles(index_t n) const {
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
        doublereal speciesMoles(index_t k) const;

        /// Index of the species belonging to phase number \c p
        /// with index \c k within the phase.
        int speciesIndex(index_t k, index_t p) const {
            return m_spstart[p] + k;
        }

        /// Minimum temperature for which all solution phases have
        /// valid thermo data. Stoichiometric phases are not
        /// considered, since they may have thermo data only valid for
        /// conditions for which they are stable.
        doublereal minTemp() const { return m_Tmin; }

        /// Maximum temperature for which all solution phases have
        /// valid thermo data. Stoichiometric phases are not
        /// considered, since they may have thermo data only valid for
        /// conditions for which they are stable.
        doublereal maxTemp() const { return m_Tmax; }

        /// Total charge (Coulombs). 
        doublereal charge() const;

        /// Charge (Coulombs) of phase with index \a p.
        doublereal phaseCharge(index_t p) const;

        /// Total moles of element \a m, summed over all phases.
        doublereal elementMoles(index_t m) const;

        /// Chemical potentials. Write into array \a mu the chemical
        /// potentials of all species [J/kmol]. The chemical
        /// potentials are related to the activities by
        /// \f[ \mu_k = \mu_k^0(T, P) + RT \ln a_k. \f].
        void getChemPotentials(doublereal* mu) const;

        /// Valid chemical potentials. Write into array \a mu the
        /// chemical potentials of all species with thermo data valid
        /// for the current temperature [J/kmol]. For other species,
        /// set the chemical potential to the value \a not_mu. If \a
        /// standard is set to true, then the values returned are
        /// standard chemical potentials.
        void getValidChemPotentials(doublereal not_mu, doublereal* mu,
            bool standard = false) const;

        /// Temperature [K].
        doublereal temperature() const { return m_temp; }

        /// Set the mixture to a state of chemical equilibrium.
        /// @param XY Integer flag specifying properties to hold fixed.
        /// @param err Error tolerance for \f$\Delta \mu/RT \f$ for
        /// all reactions. Also used as the relative error tolerance
        /// for the outer loop.
        /// @param maxsteps Maximum number of steps to take in solving
        /// the fixed TP problem.
        /// @param maxiter Maximum number of "outer" iterations for
        /// problems holding fixed something other than (T,P).
        /// @param loglevel Level of diagnostic output, written to a
        /// file in HTML format.
        doublereal equilibrate(int XY, doublereal err = 1.0e-9, 
            int maxsteps = 1000, int maxiter = 200, int loglevel = -99);


        /// Set the temperature [K].
        void setTemperature(doublereal T) {
            m_temp = T;
            updatePhases();
        }

        /// Pressure [Pa].
        doublereal pressure() const {
            return m_press;
        }

        /// Volume [m^3].
        doublereal volume() const;

        /// Set the pressure [Pa].
        void setPressure(doublereal P) {
            m_press = P;
            updatePhases();
        }

        /// Enthalpy [J].
        doublereal enthalpy() const;

        /// Entropy [J/K].
        doublereal entropy() const;

        /// Gibbs function [J].
        doublereal gibbs() const;

        /// Heat capacity at constant pressure [J/K].
        doublereal cp() const;

        /// Number of phases.
        index_t nPhases() const {
            return m_np;
        }

        /// Return true is species \a k is a species in a
        /// multicomponent solution phase.
        bool solutionSpecies(index_t k) const;

        index_t speciesPhaseIndex(index_t k) const{
            return m_spphase[k];
        }

        doublereal moleFraction(index_t k) const{
            return m_moleFractions[k];
        }

        void setPhaseMoleFractions(index_t n, doublereal* x);

        void setMolesByName(compositionMap& xMap);

        void setMolesByName(const string& x);

        void setMoles(doublereal* n);

        /// Return true if the phase \a p has valid thermo data for
        /// the current temperature.
        bool tempOK(index_t p) const {
            return m_temp_OK[p];
        }

    protected:

        // These methods are meant for internal use.

        /// update the locally-stored composition to match the current
        /// compositions of the phase objects.
        void updateMoleFractions();

        /// Set the states of the phase objects to the locally-stored
        /// state.  Note that if individual phases have T and P different
        /// than that stored locally, the phase T and P will be modified.
        void updatePhases() const;
           
        /**
         * Vector of the number of moles in each phase. 
         * Length = m_np, number of phases.
         */
        vector_fp m_moles;

        /**
	 * Vector of the ThermoPhase Pointers.
	 */
        vector<phase_t*> m_phase;
        array_t m_atoms;
      /**
       * Locally storred vector of mole fractions of all species 
       * comprising the MultiPhase object.
       */
        vector_fp m_moleFractions;
        vector_int m_spphase;
        vector_int m_spstart;
        vector<string> m_enames;
        vector_int m_atomicNumber;
        vector<string> m_snames;
        mutable map<string, int> m_enamemap;
        /**
	 *   Number of phases in the MultiPhase object
	 */
        index_t  m_np;
        doublereal m_temp;
        doublereal m_press;
      /**
       * Number of distinct elements in all of the phases
       */
        index_t m_nel; 
      /**
       * Number of distinct species in all of the phases
       */
        index_t m_nsp;
        bool m_init;
        int m_eloc;
        mutable vector<bool> m_temp_OK;
        doublereal m_Tmin, m_Tmax;
    };

    inline std::ostream& operator<<(std::ostream& s, Cantera::MultiPhase& x) {
        size_t ip;
        for (ip = 0; ip < x.nPhases(); ip++) {
            if (x.phase(ip).name() != "") {
                s << "*************** " << x.phase(ip).name() << " *****************" << endl;
            }
            else {
                s << "*************** Phase " << ip << " *****************" << endl;
            }
            s << "Moles: " << x.phaseMoles(ip) << endl;
                
            s << report(x.phase(ip)) << endl;
        }
        return s;
    }
}

#endif
