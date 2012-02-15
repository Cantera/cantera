/**
 *  @file Species.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CKR_SPECIES_H
#define CKR_SPECIES_H

#include "ckr_defs.h"
#include "Constituent.h"
#include <map>
#include <vector>

namespace ckr
{

/**
 *   Holds species data read in from entries in the THERMO section of
 *   a chemkin or nasa9 fortran formatted input file.
 */
class Species
{
public:

    /// Construct an empty Species object
    Species();

    //! Copy constructor
    Species(const Species& s);

    /// Destructor
    ~Species();

    /// Assignment operator
    Species& operator=(const Species& s);

    /// Test for equality based on name only.
    bool operator==(const Species& s) const;

    bool operator!=(const Species& s) const;

    /// Used to sort lists of species by index number.
    bool operator<(const Species& s) const;

    //! Type of thermodynamic representation
    /*!
     *  0 This is a 2 region NASA polynomial representation
     *
     *  1 This is a multiple temperature region NASA9 polynomial
     *    representation.
     */
    int thermoFormatType;

    //! Species Name
    std::string name;
    std::string id;              //!<   ID tag from 'date' field in input
    std::string phase;           //!<   Phase string. Usually "G", "L", or "S".
    double tlow;                 //!<   Min temperature for thermo data fit
    double tmid;                 //!<   Mid temperature for thermo data fit
    double thigh;                //!<   Max temperature for thermo data fit

    /// list of Constituent objects defining elemental composition
    std::vector<Constituent> elements;

    /// map from element symbols to atom numbers
    mutable std::map<std::string, double> comp;

    /// polynomial coefficients for the lower temperature range
    vector_fp lowCoeffs;

    /// polynomial coefficients for the upper temperature range
    vector_fp highCoeffs;

    //! Number of temperature regions
    int nTempRegions;

    std::vector<vector_fp*> region_coeffs;
    vector_fp minTemps;
    vector_fp maxTemps;

    /// flag set by the validation routines
    int valid;

    /// position in the list of species in the input file
    int index;

    std::string m_commentsRef;

private:
    //! Delete private data
    void delR();
};

//! Shorthand for a list of Species
typedef std::vector<Species> speciesList;

//! A map from species names to Species objects
typedef std::map<std::string, Species> speciesTable;
}

#endif





