/**
 *  @file converters/Group.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_GROUP_H
#define CKR_GROUP_H

#include <string>
#include <vector>

using namespace std;

namespace ckr {

/**
 *   A class for groups.  
 */
class Group {
public:

    /// Construct a new empty Group object
    Group() : name("<empty>"), index(-1) {}

    Group(const string& nm) : name(nm), index(-1) {}

    /// Destructor
    ~Group() {}

    string name;                 //!<  name
    int index;                   //!<  index number
    map<string, double> comp;    //!<  elemental composition

    /**
     * Compare two Group instances for equality based on name. 
     * Primarily for internal use.
     */
    bool operator==(const Group& g) const {
        return (name == g.name);
    }
    bool operator!=(const Group& g) const {
        return !(*this == g);
    }
};

/// a list (vector) of Groups
typedef vector<Group>      groupList;

}


#endif









