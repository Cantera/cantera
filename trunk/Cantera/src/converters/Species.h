/**
 *  @file Species.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_SPECIES_H
#define CKR_SPECIES_H

#include "ckr_defs.h"
#include "Constituent.h"
#include <map>
#include <vector>

//#include "Cantera.h"

using namespace std;

namespace ckr {

/**
 *   Holds species data read in from entries in the THERMO section of
 *   the input file. 
 */
class Species {
public:

    /// Construct an empty Species object
    Species() : 
        name ("<empty>"), 
        id ("<none>"), 
        phase (""), 
        tlow(0.0), 
        tmid(0.0), 
        thigh(0.0), 
        valid(0),
        index(-1)
        {}

    /// Destructor
    ~Species() {}

    /// Assignment operator
    Species& operator=(const Species& s) {
        if (&s == this) return *this;
        name = s.name;
        id = s.id;
        phase = s.phase;
        tlow = s.tlow;
        tmid = s.tmid;
        thigh = s.thigh;
        elements = s.elements;
        comp = s.comp;
        lowCoeffs = s.lowCoeffs;
        highCoeffs = s.highCoeffs;
        valid = s.valid;
        index = s.index;
        return *this;
    }

    /// Test for equality based on name only. 
    bool operator==(const Species& s) const {
        return (s.name == name);
    }

    bool operator!=(const Species& s) const {
        return !(*this == s);
    }

    /// Used to sort lists of species by index number.
    bool operator<(const Species& s) const {
        return (index < s.index);
    }


    string name;                 //!<   Species name
    string id;                   //!<   ID tag from 'date' field in input
    string phase;                //!<   Phase string. Usually "G", "L", or "S".
    double tlow;                 //!<   Min temperature for thermo data fit
    double tmid;                 //!<   Mid temperature for thermo data fit
    double thigh;                //!<   Max temperature for thermo data fit

    /// list of Constituent objects defining elemental composition
    vector<Constituent> elements;

    /// map from element symbols to atom numbers
    mutable map<string, double> comp;

    /// polynomial coefficients for the lower temperature range
    vector_fp lowCoeffs;

    /// polynomial coefficients for the upper temperature range
    vector_fp highCoeffs;

    /// flag set by the validation routines
    int valid;

    /// position in the list of species in the input file
    int index;
};

/// A list of Species
typedef vector<Species>      speciesList;

/// A map from species names to Species objects
typedef map<string, Species> speciesTable;

}

#endif





