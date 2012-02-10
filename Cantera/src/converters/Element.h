/**
 *  @file Element.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_ELEMENT_H
#define CKR_ELEMENT_H

#include <string>
#include <vector>
#include <ostream>


namespace ckr
{

/**
 *   A class for elements.
 *   Example usage:
 *   @code
 *   Element e;
 *   e.name = "He";
 *   e.atomicWeight = 4.0026;
 *   e.comment = "helium";
 *   @endcode
 */
class Element
{
public:

    /// Construct a new empty Element object
    Element() :
        name("<empty>"),
        atomicWeight(0.0),
        valid(0),
        index(-1),
        weightFromDB(false),
        comment("")
    {}


    /// Construct a new empty Element object
    Element(const std::string& nm, double wt) :
        name(nm),
        atomicWeight(wt),
        valid(0),
        index(-1),
        weightFromDB(false),
        comment("")
    {}


    /// Destructor
    ~Element() {}

    std::string name;            //!<  Element name
    double atomicWeight;         //!<  Atomic weight in amu
    int valid;                   //!<  flag returned by validation routines
    int index;                   //!<  index number
    bool weightFromDB;           //!<  true if atomic weight is not specified
    std::string comment;         //!<  comment in input file


    /**
     * Compare two Element instances for equality based on name.
     * Primarily for internal use.
     */
    bool operator==(const Element& e) const {
        return (name == e.name);
    }
    bool operator!=(const Element& e) const {
        return !(*this == e);
    }
    friend std::ostream& operator<<(std::ostream& s, const Element& e) {
        s << e.name;
        if (!e.weightFromDB) {
            s << "/" << e.atomicWeight << "/";
        }
        if (e.comment != "") {
            s << " !" << e.comment << std::endl;
        } else {
            s << " ";
        }
        return s;
    }
};

/// a list (vector) of Elements
typedef std::vector<Element> elementList;

}


#endif

