/**
 *  @file Elements.h
 *  Contains the LookupWtElements function and the definitions of element
 *  constraint types.
 */
//  Copyright 2001  California Institute of Technology

#ifndef CT_ELEMENTS_H
#define CT_ELEMENTS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

/*!
  *  @name Types of Element Constraint Equations
  *
  *   There may be several different types of element constraints handled
  *   by the equilibrium program and by Cantera in other contexts.
  *   These defines are used to assign each constraint to one category.
  *   @{
  */

//! An element constraint that is current turned off
#define CT_ELEM_TYPE_TURNEDOFF       -1

//! Normal element constraint consisting of positive coefficients for the
//! formula matrix.
/*!
 * All species have positive coefficients within the formula matrix.
 * With this constraint, we may employ various strategies to handle
 * small values of the element number successfully.
 */
#define CT_ELEM_TYPE_ABSPOS           0

//! This refers to conservation of electrons
/*!
 * Electrons may have positive or negative values in the Formula matrix.
 */
#define CT_ELEM_TYPE_ELECTRONCHARGE   1

//! This refers to a charge neutrality of a single phase
/*!
 * Charge neutrality may have positive or negative values in the Formula matrix.
 */
#define CT_ELEM_TYPE_CHARGENEUTRALITY 2

//! Constraint associated with maintaining a fixed lattice stoichiometry in a solid
/*!
 * The constraint may have positive or negative values. The lattice 0 species will
 * have negative values while higher lattices will have positive values
 */
#define CT_ELEM_TYPE_LATTICERATIO 3

//! Constraint associated with maintaining frozen kinetic equilibria in
//! some functional groups within molecules
/*!
 *  We seek here to say that some functional groups or ionic states should be
 *  treated as if they are separate elements given the time scale of the problem.
 *  This will be abs positive constraint. We have not implemented any examples yet.
 *  A requirement will be that we must be able to add and subtract these constraints.
 */
#define CT_ELEM_TYPE_KINETICFROZEN 4

//! Constraint associated with the maintenance of a surface phase
/*!
 *  We don't have any examples of this yet either. However, surfaces only exist
 *  because they are interfaces between bulk layers. If we want to treat surfaces
 *  within thermodynamic systems we must come up with a way to constrain their total
 *  number.
 */
#define CT_ELEM_TYPE_SURFACECONSTRAINT 5

//! Other constraint equations
/*!
 * currently there are none
 */
#define CT_ELEM_TYPE_OTHERCONSTRAINT  6
//@}

//! Number indicating we don't know the entropy of the element in its most
//! stable state at 298.15 K and 1 bar.
#define ENTROPY298_UNKNOWN -123456789.

//! Function to look up an atomic weight
//! This function looks up the argument string in the database above and
//! returns the associated molecular weight.
//! The data are from the periodic table.
//!
//! Note: The idea behind this function is to provide a unified source for the
//! element atomic weights. This helps to ensure that mass is conserved.
//!     @param ename  String, Only the first 3 characters are significant
//!     @return The atomic weight of the element
//!     @exception CanteraError If a match is not found, throws a CanteraError
double LookupWtElements(const std::string& ename);

class XML_Node;

//! Object containing the elements that make up species in a phase.
/*!
 * Class Elements manages the elements that are part of a
 * chemistry specification.  This class may support calculations
 * employing Multiple phases. In this case, a single Elements object may
 * be shared by more than one Constituents class. Reactions between
 * the phases may then be described using stoichiometry base on the
 * same Elements class object.
 *
 * The member functions return information about the elements described
 * in a particular instantiation of the class.
 *
 * @ingroup phases
 * @deprecated. Functionality is now part of class Phase. To be removed after
 *     Cantera 2.2.
 */
class Elements
{

public:

    //! Default constructor for the elements class
    Elements();

    //! Default destructor for the elements class
    ~Elements();


    //! copy constructor
    /*!
     *   This copy constructor just calls the assignment operator for this
     *   class. It sets the number of subscribers to zer0.
     *
     * @param right   Reference to the object to be copied.
     */
    Elements(const Elements& right);

    //! Assignment operator
    /*!
     *   This is the assignment operator for the Elements class.
     *   Right now we pretty much do a straight uncomplicated
     *   assignment. However, subscribers are not mucked with, as they
     *   have to do with the address of the object to be subscribed to
     *
     * @param right   Reference to the object to be copied.
     */
    Elements& operator=(const Elements& right);


    //!  Static function to look up an atomic weight
    /*!
     *   This static function looks up the argument string in the
     *   database above and returns the associated molecular weight.
     *   The data are from the periodic table.
     *
     *  Note: The idea behind this function is to provide a unified
     *        source for the element atomic weights. This helps to
     *        ensure that mass is conserved.
     *
     *  @param ename  String, Only the first 3 characters are significant
     *
     *  @return
     *    Return value contains the atomic weight of the element
     *    If a match for the string is not found, a value of -1.0 is
     *    returned.
     *
     *  @exception CanteraError
     *    If a match is not found, a CanteraError is thrown as well
     */
    static double LookupWtElements(const std::string& ename);
    /// Atomic weight of element m.
    /*!
     *  @param m element index
     */
    doublereal atomicWeight(int m) const {
        return m_atomicWeights[m];
    }

    /// Atomic number of element m.
    /*!
     *  @param m element index
     */
    int atomicNumber(int m) const {
        return m_atomicNumbers[m];
    }

    //! Entropy at 298.15 K and 1 bar of stable state
    //! of the element
    /*!
     *   units J kmol-1 K-1
     *
     *  @param m  Element index
     */
    doublereal entropyElement298(int m) const;

    //! Return the element constraint type
    /*!
     * Possible types include:
     *
     * CT_ELEM_TYPE_ABSPOS           0
     * CT_ELEM_TYPE_ELECTRONCHARGE   1
     * CT_ELEM_TYPE_CHARGENEUTRALITY 2
     * CT_ELEM_TYPE_LATTICERATIO 3
     * CT_ELEM_TYPE_KINETICFROZEN 4
     * CT_ELEM_TYPE_SURFACECONSTRAINT 5
     * CT_ELEM_TYPE_OTHERCONSTRAINT  6
     *
     * The default is  CT_ELEM_TYPE_ABSPOS
     *
     *  @param m  Element index
     *
     *  @return Returns the element type
     */
    int elementType(int m) const;

    //! Change the element type of the mth constraint
    /*!
     *  Reassigns an element type
     *
     *  @param m  Element index
     *  @param elem_type New elem type to be assigned
     *
     *  @return Returns the old element type
     */
    int changeElementType(int m, int elem_type);

    /// vector of element atomic weights
    const vector_fp& atomicWeights() const {
        return m_atomicWeights;
    }
    /**
     * Inline function that returns the number of elements in the object.
     *
     *  @return
     *    \c int: The number of elements in the object.
     */
    int nElements() const {
        return m_mm;
    }

    //! Function that returns the index of an element.
    /*!
     * Index of element named \c name. The index is an integer
     * assigned to each element in the order it was added,
     * beginning with 0 for the first element.  If \c name is not
     * the name of an element in the set, then the value -1 is
     * returned.
     *
     * @param name String containing the index.
     */
    int elementIndex(const std::string& name) const;

    //! Name of the element with index \c m.
    /*!
     * @param m Element index. If m < 0 or m >= nElements() an exception is thrown.
     */
    std::string elementName(int m) const;

    //!   Returns a string vector containing the element names
    /*!
     * Returns a read-only reference to the vector of element names.
     * @return <tt> const vector<string>& </tt>: The vector contains
     *         the element names in their indexed order.
     */
    const std::vector<std::string>& elementNames() const {
        return m_elementNames;
    }

    //! Add an element to the current set of elements in the current object.
    /*!
     *   The default weight is a special value, which will cause the
     *   routine to look up the actual weight via a string lookup.
     *
     *   There are two interfaces to this routine. The XML interface
     *   looks up the required parameters for the regular interface
     *   and then calls the base routine.
     *
     * @param symbol string symbol for the element.
     * @param weight Atomic weight of the element. If no argument
     *               is provided, a lookup is attempted.
     */
    void addElement(const std::string& symbol,
                    doublereal weight = -12345.0);
    //! Add an element to the current set of elements in the current object.
    /*!
     * @param   e   Reference to the XML_Node containing the element information
     *              The node name is the element symbol and the atomWt attribute
     *              is used as the atomic weight.
     */
    void addElement(const XML_Node& e);

    //!  Add an element only if the element hasn't been added before.
    /*!
     * This is accomplished via a string match on symbol.
     *
     * @param symbol string symbol for the element.
     * @param weight Atomic weight of the element. If no argument
     *               is provided, a lookup is attempted.
     * @param atomicNumber defaults to 0
     * @param entropy298  Value of the entropy at 298 and 1 bar of the
     *                    element in its most stable form.
     *                    The default is to specify an ENTROPY298_UNKNOWN value,
     *                    which will cause a throw error if its ever
     *                    needed.
     *  @param elem_type  New elem type to be assigned.
     *                    The default is a regular element, CT_ELEM_TYPE_ABSPOS
     */
    void addUniqueElement(const std::string& symbol,
                          doublereal weight = -12345.0, int atomicNumber = 0,
                          doublereal entropy298 = ENTROPY298_UNKNOWN, int elem_type = CT_ELEM_TYPE_ABSPOS);

    //! Add an element to the current set of elements in the current object.
    /*!
     * @param   e   Reference to the XML_Node containing the element information
     *              The node name is the element symbol and the atomWt attribute
     *              is used as the atomic weight.
     */
    void addUniqueElement(const XML_Node& e);

    //! Add multiple elements from a XML_Node phase description
    /*!
     *  @param phase XML_Node reference to a phase
     */
    void addElementsFromXML(const XML_Node& phase);

    //! Prohibit addition of more elements, and prepare to add species.
    void freezeElements();

    //! True if freezeElements has been called.
    bool elementsFrozen() const;

    /// Remove all elements
    void clear();

    /// True if both elements and species have been frozen
    bool ready() const;

    //! subscribe to this object
    /*!
     *  Increment by one the number of subscriptions to this object.
     */
    void subscribe();

    //! unsubscribe to this object
    /*!
     *  decrement by one the number of subscriptions to this object.
     */
    int unsubscribe();

    //! report the number of subscriptions
    int reportSubscriptions() const;

protected:

    /******************************************************************/
    /*      Description of DATA in the Object                         */
    /******************************************************************/

    //!   Number of elements.
    int                            m_mm;

    /* m_elementsFrozen: */
    /**   boolean indicating completion of object
     *
     *    If this is true, then no elements may be added to the
     *    object.
     */
    bool m_elementsFrozen;

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
    std::vector<std::string>                 m_elementNames;

    //! Entropy at 298.15 K and 1 bar of stable state
    /*!
     *   units J kmol-1
     */
    vector_fp m_entropy298;

    //! Vector of element types
    vector_int m_elem_type;
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
    static std::vector<Elements*> Global_Elements_List;

    friend class Constituents;
};


}  // namespace

#endif
