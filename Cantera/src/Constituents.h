
///  @file Constituents.h
///  Header file for class Constituents


/*  $Author$ 
 *  $Date$
 *  $Revision$
 *
 *  $Log$
 *  Revision 1.6  2004-06-09 00:59:24  dggoodwin
 *  moved functions to read XML input files out of importCTML.cpp and into specific classes
 *
 *  Revision 1.5  2003/11/12 18:58:17  dggoodwin
 *  *** empty log message ***
 *
 *  Revision 1.4  2003/09/03 18:15:50  hkmoffa
 *  Added a vector get for the atoms in a species.
 *
 *  Revision 1.3  2003/07/21 16:02:53  hkmoffa
 *  Took out a double nested @name that gave a warning to doxygen
 *
 *  Revision 1.1.1.1  2003/04/14 17:57:51  dggoodwin
 *  Initial import.
 *
 *  Revision 1.24  2002/12/19 15:19:32  dgg
 *  added log block, replaced include statement for Elements.h with
 *  forward reference for Elements class
 *
 */

//  Copyright 2001  California Institute of Technology


#ifndef CT_CONSTIT_H
#define CT_CONSTIT_H


#include "ct_defs.h"
using namespace std;

#include "SpeciesThermo.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "xml.h"

namespace Cantera {

    class Elements;


    /************** DEFINITIONS OF ERRORS *****************************/

    class SpeciesRangeError : public CanteraError {
    public:
        SpeciesRangeError(string func, int k, int kmax) :
            CanteraError(func, "Species index " + int2str(k) + 
                " outside valid range of 0 to " + int2str(kmax-1)) {}
    };

    /******************************************************************/


    /// Class Constituents manages a set of elements and
    /// species. Class Constituents is designed to provide information
    /// about the elements and species in a phase - names, index
    /// numbers (location in arrays), atomic or molecular weights,
    /// etc. No computations are performed by the methods of this
    /// class. The set of elements must include all those that compose
    /// the species, but may include additional elements. The species
    /// all must belong to the same phase.

    class Constituents {

    public:

      /// Constructor.
      Constituents(Elements* ptr_Elements = 0);

      /// Destructor. 
      ~Constituents();

      /// @name Element Information
      //@{
        
      /// Name of the element with index m.  @param m Element
      /// index. If m < 0 or m >= nElements() an exception is thrown.
      string elementName(int m) const;


      /// Index of element named 'name'. The index is an integer
      /// assigned to each element in the order it was added,
      /// beginning with 0 for the first element.  If 'name' is not
      /// the name of an element in the set, then the value -1 is
      /// returned.
      int elementIndex(string name) const;


      /// Atomic weight of element m. 
      doublereal atomicWeight(int m) const;


      /// Return a read-only reference to the vector of element names.
      const vector<string>& elementNames() const;


      /// Return a read-only reference to the vector of atomic weights.
      const array_fp& atomicWeights() const;


      /// Number of elements.
      int nElements() const;

      //@}



      /// @name Adding Elements and Species
      /// These methods are used to add new elements or species.
      /// These are not usually called by user programs.
      /// 
      /// Since species are checked to insure that they are only 
      /// composed of declared elements, it is necessary to first
      /// add all elements before adding any species. 

      //@{

      /// Add an element. 
      /// @param symbol Atomic symbol string.
      /// @param weight Atomic mass in amu.
      void addElement(const string& symbol, doublereal weight);

      /// Add an element from an XML specification.
      void addElement(const XML_Node& e);

      void addUniqueElement(const string& symbol, doublereal weight);

      void addUniqueElement(const XML_Node& e);

      void addElementsFromXML(const XML_Node& phase);
   
      /// Prohibit addition of more elements, and prepare to add
      /// species.
      void freezeElements();

      /// True if freezeElements has been called.
      bool elementsFrozen();

      //@}
      
	/// Returns the number of species in the phase
        int nSpecies() const { return m_kk; }

        /// Molecular weight of species k.
        doublereal molecularWeight(int k) const;

	/// Molar mass. Preferred name for molecular weight.
	doublereal molarMass(int k) const {
	  return molecularWeight(k);
	}

        /**
	 * Return a const reference to the vector of molecular weights
	 * of the species
	 */
        const array_fp& molecularWeights() const;

        
	///   Electrical charge of one species k molecule, divided by
	///   the magnitude of the electron charge ( \f$ e = 1.602
	///   \times 10^{-19}\f$ Coulombs). Dimensionless.
        doublereal charge(int k) const;

       /**
         * @name Adding Species
         * These methods are used to add new species.
         * They are not usually called by user programs.
         */
        //@{
        void addSpecies(const string& name, const doublereal* comp,
			doublereal charge = 0.0, doublereal size = 1.0);

	void addUniqueSpecies(const string& name, const doublereal* comp,
			      doublereal charge = 0.0, 
			      doublereal size = 1.0);
        /**
         * Index of species named 'name'. The first species added
         * will have index 0, and the last one index nSpecies() - 1.
         */
        int speciesIndex(string name) const;
        /// Name of the species with index k
        string speciesName(int k) const;
        /// Return a const referernce to the vector of species names
        const vector<string>& speciesNames() const;
        /**
	 * size():
	 *   This routine returns the size of species k
	 */
        doublereal size(int k) const { return m_speciesSize[k]; }

#ifdef INCL_DEPRECATED_METHODS
        /**
         * Return a SpeciesData structure containing species data. 
         */
        SpeciesData species(int k) const {
            if (k < 0 || k >= nSpecies()) 
                throw SpeciesRangeError("Constituents::charge",k,nSpecies());
            SpeciesData s;
            s.name = m_speciesNames[k];
            s.phase = m_speciesPhase[k];
            int offset = m_mm * k;
            s.atoms.resize(m_mm);
            for (int m = 0; m < m_mm; m++) {
                s.atoms[m] = m_speciesComp[offset + m];
            }
            s.charge = m_speciesCharge[k];
            s.molecularWeight = m_weight[k];
            return s;
        }
#endif

        /**
         * Prohibit addition of more species, and prepare for
         * calculations with this set of elements and species.
         */
        void freezeSpecies();

        /// True if freezeSpecies has been called.
        bool speciesFrozen() { return m_speciesFrozen; }

        /// Remove all elements and species
        void clear();

        //@}

        /// True if both elements and species have been frozen
        bool ready() const;

        /// Number of atoms of element m in species k.
        doublereal nAtoms(int k, int m) const;

	/**
	 * Get a vector containing the atomic composition 
	 * of species k
	 */
	void getAtoms(int k, double *atomArray) const;

	/// Copy constructor and assignment operator
        Constituents::Constituents(const Constituents& right);
        Constituents& operator=(const Constituents& right);

    protected:
    
        int                            m_kk;
        vector_fp                      m_weight;

        bool                           m_speciesFrozen;

        /*
         * Pointer to the element object corresponding to this
         * phase. Normally, this will be the default Element object
         * common to all phases.
         */
        Elements *                     m_Elements;

        vector<string>                 m_speciesNames;
        vector_fp                      m_speciesComp;
        /**
	 * m_speciesCharge: Vector of species charges
	 *           length = m_kk
	 */
        vector_fp                      m_speciesCharge;
        /**
	 * m_speciesSize(): Vector of species sizes.
	 *           length m_kk
	 *           This is used in some equations of state
	 *           which employ the constant partial molar
	 *           volume approximation. It's so fundamental
	 *           we've put it at the Constituents class level
	 */
        vector_fp                      m_speciesSize;

    private:

    };  
    
    
}  // namespace

#endif
