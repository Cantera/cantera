/**
 * @file Constituents.h
 *  Header file  Class \link Cantera::Constituents Constitutents\endlink which 
 *  manages a set of elements and species (see \ref phases).
 */

/*  
 *  $Date: 2009/02/23 20:41:03 $
 *  $Revision: 1.7 $
 */

//  Copyright 2001  California Institute of Technology


#ifndef CT_CONSTIT_H
#define CT_CONSTIT_H


#include "ct_defs.h"

#include "SpeciesThermo.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "xml.h"
#include "Elements.h"

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


  //! Class %Constituents manages a set of elements and  species.
  /*!
   * Class %Constituents is designed to provide information
   * about the elements and species in a phase - names, index
   * numbers (location in arrays), atomic or molecular weights,
   * etc. No computations are performed by the methods of this
   * class. The set of elements must include all those that compose
   * the species, but may include additional elements. The species
   * all must belong to the same phase.
   *
   * @ingroup phases
   */
  class Constituents {

  public:

    //! Constructor.
    /*!
     *  Constructor sets all base variable types to zero. Also, it
     *  sets the pointer to the Elements object for this object.
     *
     * @param ptr_Elements
     *   The default is that a new Elements object is created, so this
     *   Constituents object is independent of any other object. But if
     *   ptr_Elements is supplied, it will be used. This way, a class
     *   implementing a multi-phase mixture is responsible for
     *   maintaining the global elements list for the mixture, and no
     *   static global element list is required.
     */
    Constituents(Elements* ptr_Elements = 0);
       
    /// Destructor. 
    ~Constituents();

    ///  This copy constructor just calls the assignment operator
    ///  for this class.
    /*!
     * @param right     reference to the object to be copied.
     */
    Constituents(const Constituents& right);

    /// Assignment operator
    /*!
     *  @param right    Reference to the object to be copied.
     */
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
    /*!
     * @param m  Element index
     */
    doublereal atomicWeight(int m) const;

    /// Entropy of the element in its standard state at 298 K and 1 bar
    /*!
     * @param m  Element index
     */
    doublereal entropyElement298(int m) const;

    /// Atomic number of element m.
    /*!
     *  @param m Element index
     */
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

    //! Add an element. 
    /*!
     *  @param symbol Atomic symbol std::string.
     *  @param weight Atomic mass in amu.
     */
    void addElement(const std::string& symbol, doublereal weight);
      
    //! Add an element from an XML specification.
    /*!
     * @param e Reference to the XML_Node where the element is described.
     */
    void addElement(const XML_Node& e);

    //! Add an element, checking for uniqueness
    /*!
     * The uniqueness is checked by comparing the string symbol. If
     * not unique, nothing is done.
     *
     * @param symbol  String symbol of the element
     * @param weight  Atomic weight of the element (kg kmol-1).
     * @param atomicNumber Atomic number of the element (unitless)
     * @param entropy298 Entropy of the element at 298 K and 1 bar
     *                   in its most stable form. The default is
     *                   the value ENTROPY298_UNKNOWN, which is 
     *                   interpreted as an unknown, and if used
     *                   will cause Cantera to throw an error.
     */
    void addUniqueElement(const std::string& symbol, doublereal weight,
			  int atomicNumber = 0,
			  doublereal entropy298 = ENTROPY298_UNKNOWN);

    //! Adde an element, checking for uniqueness
    /*!
     * The uniqueness is checked by comparing the string symbol. If
     * not unique, nothing is done.
     *
     * @param e Reference to the XML_Node where the element is described.
     */
    void addUniqueElement(const XML_Node& e);

    //! Add all elements referenced in an XML_Node tree
    /*!
     * @param phase Reference to the top  XML_Node of a phase
     */
    void addElementsFromXML(const XML_Node& phase);
   
    /// Prohibit addition of more elements, and prepare to add species.
    void freezeElements();

    /// True if freezeElements has been called.
    bool elementsFrozen();

    //@}
      
    /// Returns the number of species in the phase
    int nSpecies() const { return m_kk; }

    //! Molecular weight of species \c k.
    /*!
     * @param k   index of species \c k
     * @return
     *      Returns the molecular weight of species \c k.
     */
    doublereal molecularWeight(int k) const;

    //! Return the Molar mass of species \c k
    /*!
     * Preferred name for molecular weight.
     *
     * @param k  index for species
     * @return
     *      Return the molar mass of species k kg/kmol.
     */
    doublereal molarMass(int k) const {
      return molecularWeight(k);
    }

    /**
     * Return a const reference to the vector of molecular weights
     * of the species
     */
    const vector_fp& molecularWeights() const;
       
    /*!
     *   Electrical charge of one species k molecule, divided by
     *   the magnitude of the electron charge ( \f$ e = 1.602
     *   \times 10^{-19}\f$ Coulombs). Dimensionless.
     *
     * @param k species index
     */
    doublereal charge(int k) const;

    /**
     * @name Adding Species
     * These methods are used to add new species.
     * They are not usually called by user programs.
     */
    //@{
    void addSpecies(const std::string& name, const doublereal* comp,
		    doublereal charge = 0.0, doublereal size = 1.0);

    //! Add a species to the phase, checking for uniqueness of the name
    /*!
     * This routine checks for uniqueness of the string name. It only
     * adds the species if it is unique.
     *
     * @param name    String name of the species
     * @param  comp   Double vector containing the elemental composition of the
     *                species.
     * @param charge  Charge of the species. Defaults to zero.
     * @param size    Size of the species (meters). Defaults to 1 meter.
     */
    void addUniqueSpecies(const std::string& name, const doublereal* comp,
			  doublereal charge = 0.0, 
			  doublereal size = 1.0);
      
    //! Index of species named 'name'
    /*!
     * The first species added
     * will have index 0, and the last one index nSpecies() - 1.
     *
     * @param name String name of the species
     * @return 
     *    Returns the index of the species.
     */
    int speciesIndex(std::string name) const;

    //! Name of the species with index k
    /*!
     * @param k index of the species
     */
    std::string speciesName(int k) const;
       
    /// Return a const referernce to the vector of species names
    const std::vector<std::string>& speciesNames() const;
      
    //!  This routine returns the size of species k
    /*!
     * @param k index of the species 
     * @return 
     *      Returns the size of the species. Units are meters.
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

    //! Number of atoms of element \c m in species \c k.
    /*!
     * @param k    species index
     * @param m    element index
     */
    doublereal nAtoms(int k, int m) const;

    //! Get a vector containing the atomic composition  of species k
    /*!
     * @param  k         species index
     * @param atomArray  vector containing the atomic number in the species.
     *                   Length: m_mm
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
