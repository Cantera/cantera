/**
 *  @file Species.cpp
 *
 */

// Copyright 2001 California Institute of Technology


#include "Species.h"
#include <iostream>
#include <stdio.h>

namespace ckr
{

// Construct an empty Species object
Species::Species() :
    thermoFormatType(0),
    name("<empty>"),
    id("<none>"),
    phase(""),
    tlow(0.0),
    tmid(0.0),
    thigh(0.0),
    nTempRegions(2),
    valid(0),
    index(-1)
{
}

// Destructor
Species::~Species()
{
    delR();
}

void Species::delR()
{
    for (size_t i = 0; i < region_coeffs.size(); i++) {
        if (region_coeffs[i]) {
            delete region_coeffs[i];
            region_coeffs[i] = 0;
        }
    }
}

//! Copy constructor
Species::Species(const Species& s)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = s;
}

// Assignment operator
Species& Species::operator=(const Species& s)
{
    if (&s == this) {
        return *this;
    }
    thermoFormatType = s.thermoFormatType;
    name = s.name;
    id = s.id;
    phase = s.phase;
    tlow = s.tlow;
    tmid = s.tmid;
    thigh = s.thigh;
    nTempRegions = s.nTempRegions;
    elements = s.elements;
    comp = s.comp;
    lowCoeffs = s.lowCoeffs;
    highCoeffs = s.highCoeffs;
    delR();
    for (size_t i = 0; i < s.region_coeffs.size(); i++) {
        region_coeffs.push_back(new vector_fp(*(s.region_coeffs[i])));
    }
    minTemps = s.minTemps;
    maxTemps = s.maxTemps;
    m_commentsRef = s.m_commentsRef;
    valid = s.valid;
    index = s.index;
    return *this;
}

// Test for equality based on name only.
bool Species::operator==(const Species& s) const
{
    return (s.name == name);
}

bool Species::operator!=(const Species& s) const
{
    return !(*this == s);
}

// Used to sort lists of species by index number.
bool Species::operator<(const Species& s) const
{
    return (index < s.index);
}





}
