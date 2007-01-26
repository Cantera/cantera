/*!
 * @file Constituents.h
 *
 *  Header file for class Constituents
 */

/*  $Author$ 
 *  $Date$
 *  $Revision$
 *
 */

//  Copyright 2001  California Institute of Technology


#ifndef CT_CONSTIT_H
#define CT_CONSTIT_H


#include "ct_defs.h"
//using namespace std;

#include "SpeciesThermo.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "xml.h"

namespace Cantera {

    class Elements;

    /************** DEFINITIONS OF ERRORS *****************************/

    //! Specific fatal error indicating that the index of a species is out of range.
    /*!
     *
     *  @ingroup errorhandling
     */
    class SpeciesRangeError : public CanteraError {
    public:
      //! Constructor
      /*!
       *  @param func Function where the error occurred.
       *  @param k    current species index value
       *  @param kmax Maximum permissible species index value. The
       *              minimum permissible species index value is assumed to be 0
       *
       */
      SpeciesRangeError(std::string func, int k, int kmax) :
	CanteraError(func, "Species index " + int2str(k) + 
		     " outside valid range of 0 to " + int2str(kmax-1)) {}
    };

    /******************************************************************/


    /// Class %Constituents manages a set of elements and
    /// species. Class %Constituents is designed to provide information
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

	///  This copy constructor just calls the assignment operator
	///  for this class.
        Constituents(const Constituents& right);

	/// Assignment operator
        Constituents& operator=(const Constituents& right);

	/// @name Element Information
	// @{
        
	/// Name of the element with index m.
	///   This is a passthrough routine to the Element object.
	///   \param m  Element index. 
	///   \exception If m < 0 or m >= nElements(), the
	///          exception, ElementRangeError, is thrown.
	std::string elementName(int m) const;


	/// Index of element named 'name'.
	/// The index is an integer
	/// assigned to each element in the order it was added,
	/// beginning with 0 for the first element. 
	/// @param name  name of the element 
	///
	/// If 'name' is not
	/// the name of an element in the set, then the value -1 is
	/// returned.
	int elementIndex(std::string name) const;


	/// Atomic weight of element m. 
	doublereal atomicWeight(int m) const;

	/// Atomic number of element m. 
	int atomicNumber(int m) const;

	/// Return a read-only reference to the vector of element names.
	const std::vector<std::string>& elementNames() const;


	/// Return a read-only reference to the vector of atomic weights.
	const vector_fp& atomicWeights() const;


	/// Number of elements.
	int nElements() const;
       
	// @}



      /// @name Adding Elements and Species
      /// These methods are used to add new elements or species.
      /// These are not usually called by user programs.
      /// 
      /// Since species are checked to insure that they are only 
      /// composed of declared elements, it is necessary to first
      /// add all elements before adding any species. 

      //@{

      /// Add an element. 
      /// @param symbol Atomic symbol std::string.
      /// @param weight Atomic mass in amu.
      void addElement(const std::string& symbol, doublereal weight);

      /// Add an element from an XML specification.
      void addElement(const XML_Node& e);

      void addUniqueElement(const std::string& symbol, doublereal weight);

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
        const vector_fp& molecularWeights() const;

        
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
        void addSpecies(const std::string& name, const doublereal* comp,
			doublereal charge = 0.0, doublereal size = 1.0);

	void addUniqueSpecies(const std::string& name, const doublereal* comp,
			      doublereal charge = 0.0, 
			      doublereal size = 1.0);
        /**
         * Index of species named 'name'. The first species added
         * will have index 0, and the last one index nSpecies() - 1.
         */
        int speciesIndex(std::string name) const;
        /// Name of the species with index k
        std::string speciesName(int k) const;
        /// Return a const referernce to the vector of species names
        const std::vector<std::string>& speciesNames() const;
        /**
	 * size():
	 *   This routine returns the size of species k
	 */
        doublereal size(int k) const { return m_speciesSize[k]; }

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


    protected:
    
        //! Number of species in the phase.
      int                            m_kk;
      //! Vector of molecular weights of the species
      /*!
       * This vector has length m_kk.
       * The units of the vector are kg kmol-1.
       */
      vector_fp                      m_weight;

      //! Boolean indicating whether the number of species has been frozen.
      /*!
       * During the construction of the phase, this is false. After
       * construction of the the phase, this is true.
       */
      bool                           m_speciesFrozen;

        /*!
         * Pointer to the element object corresponding to this
         * phase. Normally, this will be the default Element object
         * common to all phases.
         */
        Elements *                     m_Elements;

      //! Vector of the species names
      std::vector<std::string>                 m_speciesNames;

      //! Atomic composition of the species.
      /*!
       * the number of atoms of i in species k is equal to
       * m_speciesComp[k * m_mm + i]
       * The length of this vector is equal to m_kk * m_mm
       */
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
