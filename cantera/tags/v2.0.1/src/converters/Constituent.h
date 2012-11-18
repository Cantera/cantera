/**
 *  @file Constituent.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_CONSTITUENT_H
#define CKR_CONSTITUENT_H

#include <string>
#include <vector>

namespace ckr
{

/**
 * A class for components of a composite object. The only members are
 * a string identifying the component, and a double specifying the
 * amount.
 */
class Constituent
{
public:
    std::string    name;        //!< The name of the object.
    double         number;      //!< The number of units (molecules, etc.).
};

}

#endif
