/**
 *  @file Phase.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
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
       *
       * @param  c       Reference to the class to be used in the copy
       */
      Phase(const Phase &c);
      
      /**
       * Assignment operator
       *
       * @param  c       Reference to the class to be used in the copy
       */
      const Phase &operator=(const Phase &c);
        
      //! Returns a reference to the XML_Node storred for the phase
      /*!
       *  The XML_Node for the phase contains all of the input data used
       *  to set up the model for the phase, during its initialization.
       */
        XML_Node& xml() { return *m_xml; }

      //! Return the string id for the phase
        std::string id() const { return m_id; }

      //! Set the string id for the phase
      /*!
       * @param id String id of the phase
       */
      void setID(std::string id) {m_id = id;} 

      //! Return the name of the phase
        std::string name() const { return m_name; }

      //! Sets the string name for the phase
      /*!
       * @param nm String name of the phase
       */
      void setName(std::string nm) { m_name = nm; }

      //! Returns the index of the phase
        int index() const { return m_index; }

      //! Sets the index of the phase
      /*!
       * @param m Integer index of the phase
       */
      void setIndex(int m) { m_index = m; }

      //! Save the current internal state of the phase
      /*!
       * Write to vector 'state' the current internal state.
       *
       * @param state output vector. Will be resized to nSpecies() + 2 on return.
       */
      void saveState(vector_fp& state) const;
       
      //! Write to array 'state' the current internal state.
      /*!
       * @param lenstate length of the state array. Must be >= nSpecies() + 2
       * @param state    output vector. Must be of length  nSpecies() + 2 or
       *                 greater.
       */
      void saveState(int lenstate, doublereal* state) const;
      
      //!Restore a state saved on a previous call to saveState.
      /*!
       * @param state State vector containing the previously saved state.
       */
      void restoreState(const vector_fp& state);

      //! Restore the state of the phase from a previously saved state vector.
      /*!
       *  @param lenstate   Length of the state vector
       *  @param state      Vector of state conditions.
       */
        void restoreState(int lenstate, const doublereal* state);

        /**
         * Set the species mole fractions by name. 
         * @param xMap map from species names to mole fraction values.
         * Species not listed by name in \c xMap are set to zero.
         */
        void setMoleFractionsByName(compositionMap& xMap);

      //! Set the mole fractions of a group of species by name
      /*!
       * The string x is in the form of a composition map
       * Species which are not listed by name in the composition
       * map are set to zero.
       *
       * @param x string x in the form of a composition map
       */
        void setMoleFractionsByName(const std::string& x);

        /**
         * Set the species mass fractions by name. 
         * @param yMap map from species names to mass fraction values.
         * Species not listed by name in \c yMap are set to zero.
         */
        void setMassFractionsByName(compositionMap& yMap);

    
      //! Set the species mass fractions by name. 
      /*!
       * Species not listed by name in \c x are set to zero.
       *
       * @param x  String containing a composition map
       */
      void setMassFractionsByName(const std::string& x);

      //! Set the internally storred temperature (K), density, and mole fractions.  
      /*!
       * Note, the mole fractions are always set first, before the density
       *
       * @param t     Temperature in kelvin
       * @param dens  Density (kg/m^3)
       * @param x     vector of species mole fractions.
       *              Length is equal to m_kk
       */
      void setState_TRX(doublereal t, doublereal dens, const doublereal* x);


      //! Set the internally storred temperature (K), density, and mole fractions.  
      /*!
       * Note, the mole fractions are always set first, before the density
       *
       * @param t     Temperature in kelvin
       * @param dens  Density (kg/m^3)
       * @param x     Composition Map containing the mole fractions.
       *              Species not included in the map are assumed to have
       *              a zero mole fraction.
       */
      void setState_TRX(doublereal t, doublereal dens, compositionMap& x);

      //! Set the internally storred temperature (K), density, and mass fractions.  
      /*!
       * Note, the mass fractions are always set first, before the density
       *
       * @param t     Temperature in kelvin
       * @param dens  Density (kg/m^3)
       * @param y     vector of species mass fractions.
       *              Length is equal to m_kk
       */
      void setState_TRY(doublereal t, doublereal dens, const doublereal* y);

      //! Set the internally storred temperature (K), density, and mass fractions.  
      /*!
       * Note, the mass fractions are always set first, before the density
       *
       * @param t     Temperature in kelvin
       * @param dens  Density (kg/m^3)
       * @param y     Composition Map containing the mass fractions.
       *              Species not included in the map are assumed to have
       *              a zero mass fraction.
       */
      void setState_TRY(doublereal t, doublereal dens, compositionMap& y);

      //! Set the internally storred temperature (K), molar density (kmol/m^3), and mole fractions.  
      /*!
       * Note, the mole fractions are always set first, before the molar density
       *
       * @param t     Temperature in kelvin
       * @param n     molar density (kmol/m^3)
       * @param x     vector of species mole fractions.
       *              Length is equal to m_kk
       */
      void setState_TNX(doublereal t, doublereal n, const doublereal* x);
    
      //! Set the internally storred temperature (K) and density (kg/m^3)
      /*!
       * @param t     Temperature in kelvin
       * @param rho   Density (kg/m^3)
       */
      void setState_TR(doublereal t, doublereal rho);
      
      //! Set the internally storred temperature (K) and mole fractions.  
      /*!
       * @param t   Temperature in kelvin
       * @param x   vector of species mole fractions.
       *            Length is equal to m_kk
       */
      void setState_TX(doublereal t, doublereal* x);

      //! Set the internally storred temperature (K) and mass fractions.  
      /*!
       * @param t   Temperature in kelvin
       * @param y   vector of species mass fractions.
       *            Length is equal to m_kk
       */
      void setState_TY(doublereal t, doublereal* y);

      //! Set the density (kg/m^3) and mole fractions. 
      /*!
       * @param rho  Density (kg/m^3)
       * @param x    vector of species mole fractions.
       *             Length is equal to m_kk
       */
      void setState_RX(doublereal rho, doublereal* x);

      //! Set the density (kg/m^3) and mass fractions.  
      /*!
       * @param rho  Density (kg/m^3)
       * @param y    vector of species mass fractions.
       *             Length is equal to m_kk
       */
      void setState_RY(doublereal rho, doublereal* y);

      /**
       * Copy the vector of molecular weights into vector weights.
       *
       * @param weights Output vector of molecular weights (kg/kmol)
       */
      void getMolecularWeights(vector_fp& weights);

      /**
       * Copy the vector of molecular weights into array weights.
       *
       * @param iwt      Unused. 
       * @param weights  Output array of molecular weights (kg/kmol)
       *
       * @deprecated
       */
      void getMolecularWeights(int iwt, doublereal* weights);

      /**
       * Copy the vector of molecular weights into array weights.
       *
       * @param weights  Output array of molecular weights (kg/kmol)
       */
      void getMolecularWeights(doublereal* weights);

      /**
       * Return a const reference to the internal vector of
       * molecular weights.
       */
      const array_fp& molecularWeights();

      /**
       * Get the mole fractions by name. 
       *
       * @param x  Output composition map containing the
       *           species mole fractions.
       */
      void getMoleFractionsByName(compositionMap& x);

      //! Return the mole fraction of a single species
      /*!
       * @param  k  String name of the species
       *
       * @return Mole fraction of the species
       */
      doublereal moleFraction(int k) const;
        
      //! Return the mole fraction of a single species
      /*!
       * @param  name  String name of the species
       *
       * @return Mole fraction of the species
       */
      doublereal moleFraction(std::string name) const;

      //! Return the mass fraction of a single species
      /*!
       * @param  k  String name of the species
       *
       * @return Mass Fraction of the species
       */
      doublereal massFraction(int k) const;

      //! Return the mass fraction of a single species
      /*!
       * @param  name  String name of the species
       *
       * @return Mass Fraction of the species
       */
      doublereal massFraction(std::string name) const;

        /**
         * Charge density [C/m^3].
         */
        doublereal chargeDensity() const;

      /// Returns the number of spatial dimensions (1, 2, or 3)
      int nDim() {return m_ndim;}

      //! Set the number of spatial dimensions (1, 2, or 3)
      /*!
       *  The number of spatial dimensions is used for vector involving
       *  directions.
       *
       * @param ndim   Input number of dimensions.
       */
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
        std::string m_id;
        std::string m_name;
    };

  //! typedef for the base Phase class
    typedef Phase phase_t;
}

#endif
