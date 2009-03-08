/**
 *  @file Phase.h
 */

/*
 * $Author: hkmoffa $
 * $Revision: 1.13 $
 * $Date: 2006/06/23 20:35:16 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PHASE_H
#define CT_PHASE_H

#include "State.h"
#include "Constituents.h"
#include "vec_functions.h"

#include "ctml.h"
using namespace ctml;

namespace Cantera {
    

    /** 
     * @defgroup phases Phases of Matter
     *
     * These classes are used to represent phases of matter.
     */


    /**
     * Base class for phases of matter. Class Phase derives from both
     * Constituents and State. In addition to the methods of those two
     * classes, it implements methods that allow referencing a species
     * by name.
     * @ingroup phases
     */
    class Phase : public Constituents, public State {    

    public:

        /// Default constructor.
        Phase() : m_kk(-1), m_ndim(3), m_index(-1), 
                  m_xml(new XML_Node("phase")), 
                  m_id("<phase>"), m_name("") {}

        /// Destructor.
        virtual ~Phase(){ 
	    delete m_xml;
	    m_xml = 0;
	}

	/**
	 * Copy Constructor
	 */
	Phase(const Phase &c);

	/**
	 * Assignment operator
	 */
	const Phase &operator=(const Phase &c);
        
        XML_Node& xml() { return *m_xml; }
        string id() const { return m_id; }
        void setID(string id) {m_id = id;} 

        string name() const { return m_name; }
        void setName(string nm) { m_name = nm; }

        int index() const { return m_index; }
        void setIndex(int m) { m_index = m; }

        /** 
         * Write to vector 'state' the current internal state.  
         * @param state output vector. Will be resized to nSpecies() + 2 on
         * return.
         */
        void saveState(vector_fp& state) const;

        /** 
         * Write to array 'state' the current internal state.
         * @param lenstate length of the state array. Must be >= nSpecies() + 2
	 * @param state output vector. Must be of length  nSpecies() + 2 or
         *              greater.
         */
        void saveState(int lenstate, doublereal* state) const;

        /**
         * Restore a state saved on a previous call to saveState.
         */
        void restoreState(const vector_fp& state);

        void restoreState(int lenstate, const doublereal* state);

        /**
         * Set the species mole fractions by name. 
         * @param xMap map from species names to mole fraction values.
         * Species not listed by name in \c xMap are set to zero.
         */
        void setMoleFractionsByName(compositionMap& xMap);

        void setMoleFractionsByName(const string& x);

        /**
         * Set the species mass fractions by name. 
         * @param yMap map from species names to mass fraction values.
         * Species not listed by name in \c yMap are set to zero.
         */
        void setMassFractionsByName(compositionMap& yMap);

        void setMassFractionsByName(const string& x);

        /** Set the temperature (K), density (kg/m^3), and mole fractions. */
        void setState_TRX(doublereal t, doublereal dens, const doublereal* x);

        /** Set the temperature (K), density (kg/m^3), and mole fractions. */
        void setState_TRX(doublereal t, doublereal dens, compositionMap& x);

        /** Set the temperature (K), density (kg/m^3), and mass fractions. */
        void setState_TRY(doublereal t, doublereal dens, const doublereal* y);

        /** Set the temperature (K), density (kg/m^3), and mass fractions. */
        void setState_TRY(doublereal t, doublereal dens, compositionMap& y);

        /** Set the temperature (K), molar density (kmol/m^3), and mole fractions. */
        void setState_TNX(doublereal t, doublereal n, const doublereal* x);
    
        /** Set the temperature (K) and density (kg/m^3) */
        void setState_TR(doublereal t, doublereal rho);
    
        /** Set the temperature (K) and mole fractions.  */
        void setState_TX(doublereal t, doublereal* x);

        /** Set the temperature (K) and mass fractions.  */
        void setState_TY(doublereal t, doublereal* y);

        /** Set the density (kg/m^3) and mole fractions.  */
        void setState_RX(doublereal rho, doublereal* x);

        /** Set the density (kg/m^3) and mass fractions.  */
        void setState_RY(doublereal rho, doublereal* y);

        /**
         * Copy the vector of molecular weights into vector weights.
         */
        void getMolecularWeights(vector_fp& weights);

        /**
         * Copy the vector of molecular weights into array weights.
         */
        void getMolecularWeights(int iwt, doublereal* weights);

        /**
         * Copy the vector of molecular weights into array weights.
         */
        void getMolecularWeights(doublereal* weights);

        /**
         * Return a const reference to the internal vector of
         * molecular weights.
         */
        const array_fp& molecularWeights();

        /**
         * Get the mole fractions by name. 
         */
        void getMoleFractionsByName(compositionMap& x);

        doublereal moleFraction(int k) const;
        
        doublereal moleFraction(string name) const;

        doublereal massFraction(int k) const;

        doublereal massFraction(string name) const;

        /**
         * Charge density [C/m^3].
         */
        doublereal chargeDensity() const;

        /// Number of spatial dimensions (1, 2, or 3)
        int nDim() {return m_ndim;}
        void setNDim(int ndim) {m_ndim = ndim;}

        /** 
         *  Finished adding species, prepare to use them for calculation
         *  of mixture properties.
         */
        virtual void freezeSpecies();
 
        virtual bool ready() const;


    protected:

	/**
	 * m_kk = Number of species in the phase.  @internal m_kk is a
	 * member of both the State and Constituents classes.
	 * Therefore, to avoid multiple inheritance problems, we need
	 * to restate it in here, so that the declarations in the two
	 * base classes become hidden.
	 */
	int m_kk;
	/**
	 * m_ndim is the dimensionality of the phase.  Volumetric
	 * phases have dimensionality 3 and surface phases have
	 * dimensionality 2.
	 */
        int m_ndim;
	/**
	 * m_index is the index of the phase
	 * 
	 */
        int m_index;

    private:

        vector_fp m_data;
        XML_Node* m_xml;
        string m_id;
        string m_name;
    };

    typedef Phase phase_t;
}

#endif
