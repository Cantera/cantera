/***********************************************************************
 *  $RCSfile: Elements.h,v $
 *  $Author: dggoodwin $ 
 *  $Date: 2005/06/18 17:01:08 $
 *  $Revision: 1.5 $
 ***********************************************************************/
//  Copyright 2001  California Institute of Technology

/** @file Elements.h 
 *  Header file for class, Elements.
 *
 *  This file contains the declarations for the elements class.
 */

#ifndef CT_ELEMENTS_H
#define CT_ELEMENTS_H

#undef USE_DGG_CODE

#include "ct_defs.h"
//#include "ctexceptions.h"

namespace Cantera {

    class XML_Node;
    class ElementRangeError;

    /** Elements Class: Object contains the elements that make up species.
     * 
     * Class Elements manages the elements that are part of a
     * chemistry specification.  This class may support calculations
     * employing Multiple phases. In this case, a single Elements object may
     * be shared by more than one Constituents class. Reactions between
     * the phases may then be described using stoichiometry base on the
     * same Elements class object.
     * 
     * The member functions return information about the elements described
     * in a particular instantiation of the class.
     */
    class Elements {

    public:
	/// Default constructor for the elements class
        Elements();
        ~Elements();

	static double LookupWtElements(const string &);

        /// Atomic weight of element m. 
        doublereal atomicWeight(int m) const { return m_atomicWeights[m]; }

        /// Atomic number of element m. 
        int atomicNumber(int m) const { return m_atomicNumbers[m]; }

        /// vector of element atomic weights
        const vector_fp& atomicWeights() const { return m_atomicWeights; }

        /** 
	 * Inline function that returns the number of elements in the object.
	 * 
	 *  @return 
	 *    \c int: The number of elements in the object.
	 */
        int nElements() const { return m_mm; }

        /** Function that returns the index of an element.
	 *
         * Index of element named \c name. The index is an integer
         * assigned to each element in the order it was added,
         * beginning with 0 for the first element.  If \c name is not
         * the name of an element in the set, then the value -1 is
         * returned.
	 *
	 * @param name String containing the index.
         */
        int elementIndex(string name) const;

        /*
         * Name of the element with index \c m.  @param m Element 11111
         * index. If m < 0 or m >= nElements() an exception is thrown.
         */
        string elementName(int m) const;

	/* elementNames() */
        /**   Returns a string vector containing the element names  
	 *   
         * Returns a read-only reference to the vector of element names.
	 * @return <tt> const vector<string>& </tt>: The vector contains
	 *         the element names in their indexed order.
         */
        const vector<string>& elementNames() const {
            return m_elementNames;
        }

        /// Add an element.
        void addElement(const string& symbol, 
			doublereal weight = -12345.0);
        void addElement(const XML_Node& e);

	/*
	 * Add an element only if the element hasn't been added before.
	 * This is accomplished via a string match on symbol.
	 */
	void addUniqueElement(const string& symbol, 
            doublereal weight = -12345.0, int atomicNumber = 0);
        void addUniqueElement(const XML_Node& e);

        void addElementsFromXML(const XML_Node& phase);

        /**
         * Prohibit addition of more elements, and prepare to add
         * species.
         */
        void freezeElements();
    
        /// True if freezeElements has been called.
        bool elementsFrozen() { return m_elementsFrozen; }

        /// Remove all elements
        void clear();

        /// True if both elements and species have been frozen
        bool ready() const;

	Elements(const Elements& right);
        Elements& operator=(const Elements& right);


	void subscribe();
	int unsubscribe();
	int reportSubscriptions() const;

    protected:
    
	/******************************************************************/
	/*      Description of DATA in the Object                         */
	/******************************************************************/
	/* n_mm: */
	/**    Number of elements.
	 *
	 */
        int                            m_mm;

	/* m_elementsFrozen: */
	/**   boolean indicating completion of object
	 *
	 *    If this is true, then no elements may be added to the
	 *    object.
	 */
        bool                           m_elementsFrozen;

	/**
	 *  Vector of element atomic weights:
	 *
	 *   units = kg / kmol
	 */
        vector_fp                      m_atomicWeights;

	/**
	 *  Vector of element atomic numbers:
	 *
	 */
        vector_int                      m_atomicNumbers;
	
	/** Vector of strings containing the names of the elements
	 *        
	 *  Note, a string search is the primary way to identify elements.
	 */
        vector<string>                 m_elementNames;

	/**
	 * Number of Constituents Objects that use this object
	 *
	 * Number of Constituents Objects that require this Elements object
	 * to complete its definition.
	 * The destructor checks to see that this is equal to zero.
	 *  when the element object is released.
	 */
	int                            numSubscribers;

	/********* GLOBAL STATIC SECTION *************/
  
    public:
	/** Vector of pointers to Elements Objects
	 *
	 */
	static vector<Elements *> Global_Elements_List;

    };  
    
}  // namespace

#endif
