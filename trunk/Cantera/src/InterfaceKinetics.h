/**
 * @file InterfaceKinetics.h
 *
 * $Author: hkmoffa $
 * $Revision: 1.18 $
 * $Date: 2006/06/13 00:57:35 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_IFACEKINETICS_H
#define CT_IFACEKINETICS_H

#include <fstream>
#include <math.h>
#include <map>
#include <stdlib.h>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"
#include "RateCoeffMgr.h"
#include "ReactionStoichMgr.h"

namespace Cantera {

    // forward references

    class ReactionData;
    class InterfaceKineticsData;
    class ThermoPhase;
    class SurfPhase;
    class ImplicitSurfChem;


    /**
     * Holds mechanism-specific data. 
     */
    class InterfaceKineticsData {
    public:
        InterfaceKineticsData() :
            m_ROP_ok(false), 
            m_temp(0.0), m_logtemp(0.0)
            {}
        virtual ~InterfaceKineticsData(){}

        doublereal m_logp0, m_logc0;
        array_fp m_ropf, m_ropr, m_ropnet;
        //array_fp m_rfn_low, m_rfn_high;
        bool m_ROP_ok;

        doublereal m_temp, m_logtemp;
        vector_fp m_rfn;
        vector_fp m_rkcn;
    };


    ///
    ///  A kinetics manager for heterogeneous reaction mechanisms. The
    ///  reactions are assumed to occur at a 2D interface between two
    ///  3D phases.
    ///
    class InterfaceKinetics : public Kinetics {

    public:

        /**
	 * Constructor 
	 *
	 * @param thermo The optional parameter may be used to initialize
	 *               the object with one ThermoPhase object.
	 *               HKM Note -> Since the interface kinetics
	 *               object will probably require multiple thermophase
	 *               objects, this is probably not a good idea
	 *               to have this parameter.
	 */
        InterfaceKinetics(thermo_t* thermo = 0);


        /// Destructor.
        virtual ~InterfaceKinetics();

        virtual int ID() { return cInterfaceKinetics; }
        virtual int type() { return cInterfaceKinetics; }

	/**
	 * Set the electric potential in the nth phase.
         * @deprecated
	 *
	 * @param n phase Index in this kinetics object. 
	 * @param V Electric potential (volts)
	 */
//         void setElectricPotential(int n, doublereal V) {
//             thermo(n).setElectricPotential(V);
//             m_redo_rates = true;
// 	}


        ///
        ///  @name Reaction Rates Of Progress
        ///
        //@{ 


        virtual void getFwdRatesOfProgress(doublereal* fwdROP) { 
            updateROP(); 
            copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }

        virtual void getRevRatesOfProgress(doublereal* revROP) { 
            updateROP(); 
            copy(m_kdata->m_ropr.begin(), m_kdata->m_ropr.end(), revROP);
        }

        virtual void getNetRatesOfProgress(doublereal* netROP) { 
            updateROP(); 
            copy(m_kdata->m_ropnet.begin(), m_kdata->m_ropnet.end(), netROP);
        }

	virtual void getEquilibriumConstants(doublereal* kc);


	virtual void getDeltaGibbs( doublereal* deltaG);

	/**
	 * Return the vector of values for the reactions change in
	 * enthalpy.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaEnthalpy( doublereal* deltaH);

	/**
	 * Return the vector of values for the reactions change in
	 * entropy.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaEntropy(doublereal* deltaS);

	/**
	 * Return the vector of values for the reaction 
	 * standard state gibbs free energy change.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSGibbs(doublereal* deltaG);

	/**
	 * Return the vector of values for the change in the
	 * standard state enthalpies of reaction.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSEnthalpy(doublereal* deltaH);

	/**
	 * Return the vector of values for the change in the
	 * standard state entropies for each reaction.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaSSEntropy(doublereal* deltaS);


	//@}
        /**
         * @name Species Production Rates
         */
        //@{

        /**
         * Species creation rates [kmol/m^2/s]. Return the species 
         * creation rates in array cdot, which must be
         * dimensioned at least as large as the total number of
         * species in all phases of the kinetics
	 * model
         *  
         */ 
        virtual void getCreationRates(doublereal* cdot) {
            updateROP();
            m_rxnstoich.getCreationRates(m_kk, &m_kdata->m_ropf[0], 
                &m_kdata->m_ropr[0], cdot); 
        }

        /**
         * Species destruction rates [kmol/m^2/s]. Return the species 
         * destruction rates in array ddot, which must be
         * dimensioned at least as large as the total number of
         * species in all phases of the kinetics
	 * model
         *  
         */ 
        virtual void getDestructionRates(doublereal* ddot) {
            updateROP();
            m_rxnstoich.getDestructionRates(m_kk, &m_kdata->m_ropf[0], 
                &m_kdata->m_ropr[0], ddot); 
        }

	/**
         * Species net production rates [kmol/m^2/s]. Return the species
         * net production rates (creation - destruction) in array
         * wdot, which must be dimensioned at least as large as the
         * total number of species in all phases of the kinetics
	 * model
         */ 
        virtual void getNetProductionRates(doublereal* net) {
            updateROP();
             m_rxnstoich.getNetProductionRates(m_kk, 
					      &m_kdata->m_ropnet[0],
					      net); 
        }

        //@}
        /**
         * @name Reaction Mechanism Informational Query Routines
         */
        //@{

	/**
	 * Stoichiometric coefficient of species k as a reactant in
         * reaction i.  
	 */
        virtual doublereal reactantStoichCoeff(int k, int i) const {
            return m_rrxn[k][i];
        }

        /**
         * Stoichiometric coefficient of species k as a product in
         * reaction i.  
         */
        virtual doublereal productStoichCoeff(int k, int i) const {
            return m_prxn[k][i];
        }

        /**
         * Flag specifying the type of reaction. The legal values and
         * their meaning are specific to the particular kinetics
         * manager.
         */
	virtual int reactionType(int i) const {
            return m_index[i].first;
        }

	/**
         * True if reaction i has been declared to be reversible. If
         * isReversible(i) is false, then the reverse rate of progress
         * for reaction i is always zero.
         */
        virtual bool isReversible(int i) {
            if (find(m_revindex.begin(), m_revindex.end(), i) 
                < m_revindex.end()) return true;
            else return false;
        }

	/**
	 * Return a string representing the reaction.
	 */
	virtual string reactionString(int i) const {
            return m_rxneqn[i];
        }


        virtual void getFwdRateConstants(doublereal* kfwd);
        virtual void getRevRateConstants(doublereal* krev, 
					 bool doIrreversible = false);
	virtual void getActivationEnergies(doublereal *E);

	//@}
        /**
         * @name Reaction Mechanism Construction
         */
        //@{

	/**
	 * Prepare the class for the addition of reactions. This function
	 * must be called after instantiation of the class, but before
	 * any reactions are actually added to the mechanism.
	 * This function calculates m_kk the number of species in all
	 * phases participating in the reaction mechanism. We don't know
	 * m_kk previously, before all phases have been added. 
         */
        virtual void init();

        /**
	 *  Add a single reaction to the mechanism.
	 */
        virtual void addReaction(const ReactionData& r);

        /**
	 * Finish adding reactions and prepare for use. This function
	 * must be called after all reactions are entered into the mechanism
	 * and before the mechanism is used to calculate reaction rates. 
	 */
        virtual void finalize();
        virtual bool ready() const;


        void updateROP();


        //const vector<grouplist_t>& reactantGroups(int i)
        //    { return m_rgroups[i]; }
        //const vector<grouplist_t>& productGroups(int i)
        //    { return m_pgroups[i]; }

        void _update_rates_T();
        void _update_rates_phi();
        void _update_rates_C();

        void advanceCoverages(doublereal tstep);
        void checkPartialEquil();
        vector_fp m_grt;

    protected:

	/**
	 * m_kk here is the number of species in all of the phases
	 * that participate in the kinetics mechanism.
	 */
        int                                 m_kk;
        vector_int m_revindex;

        Rate1<SurfaceArrhenius>                    m_rates;        
        //Rate1<Arrhenius>                    m_rates;        
        bool                                m_redo_rates;

	/**
	 * Vector of information about reactions in the
	 * mechanism.
	 * The key is the reaction index (0 < i < m_ii).
	 * The first pair is the reactionType of the reaction.
	 * The second pair is ...
	 */
        mutable map<int, pair<int, int> >   m_index;

        vector<int> m_irrev;

        //        StoichManagerN                      m_reactantStoich;
        //StoichManagerN                      m_revProductStoich;
        //StoichManagerN                      m_irrevProductStoich;

        //StoichManagerN                      m_globalReactantStoich;
        ReactionStoichMgr                   m_rxnstoich;

        int m_nirrev;

	/**
	 * Number of reversible reactions in the mechanism
	 */
        int m_nrev;

        //        map<int, vector<grouplist_t> >      m_rgroups;
        //map<int, vector<grouplist_t> >      m_pgroups;

        vector<int>                         m_rxntype;

	/**
	 *  m_rrxn is a vector of maps. m_rrxn has a length
	 *  equal to the total number of species in the kinetics
	 *  object. For each species, there exists a map, with the 
	 *  reaction number being the key, and the
	 *  reactant stoichiometric coefficient being the value.
	 *  HKM -> mutable because search sometimes creates extra
	 *         entries. To be fixed in future...
	 */
        mutable vector<map<int, doublereal> >     m_rrxn;

	/**
	 *  m_rrxn is a vector of maps. m_rrxn has a length
	 *  equal to the total number of species in the kinetics
	 *  object. For each species, there exists a map, with the 
	 *  reaction number being the key, and the
	 *  product stoichiometric coefficient being the value.
	 */
        mutable vector<map<int, doublereal> >     m_prxn;


        vector<string> m_rxneqn;

	/**
	 * Temporary data storage used in calculating the rates of
	 * of reactions.
	 */
        InterfaceKineticsData* m_kdata;

	/**
	 * An array of generalized concentrations
         * \f$ C_k \f$ that are defined such that \f$ a_k = C_k /
         * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration/
         * These generalized concentrations are used
         * by this kinetics manager class to compute the forward and
         * reverse rates of elementary reactions. The "units" for the
	 * concentrations of each phase depend upon the implementation
	 * of kinetics within that phase.
	 * The order of the species within the vector is based on
	 * the order of listed ThermoPhase objects in the class, and the
	 * order of the species within each ThermoPhase class.
	 */
        vector_fp m_conc;

        vector_fp m_mu0;
        vector_fp m_phi;
        vector_fp m_pot;
        vector_fp m_rwork;
        vector_fp m_E;

        SurfPhase*                             m_surf;
        ImplicitSurfChem*                      m_integrator;

    private:

        int reactionNumber(){ return m_ii;}
        void addElementaryReaction(const ReactionData& r);
        void addGlobalReaction(const ReactionData& r);
        void installReagents(const ReactionData& r);

        //void installGroups(int irxn, const vector<grouplist_t>& r,
        //    const vector<grouplist_t>& p);
        void updateKc();

        void registerReaction(int rxnNumber, int type, int loc) {
            m_index[rxnNumber] = pair<int, int>(type, loc);
        }
        void applyButlerVolmerCorrection(doublereal* kf);
        bool m_finalized;
        bool m_has_coverage_dependence;
    };
}

#endif
