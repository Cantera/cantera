/**
 *  @file atomicWeightDB.cpp
 *
 * internal database of default atomic weights
 *
 */

// Copyright 2001  California Institute of Technology

#include <map>
#include <string>
#include <iostream>
using namespace std;

namespace ckr
{

static double _weights[] = {
    1.00797,    4.0026,   6.939,    9.01220,   10.811,            // H - B
    12.01115,    14.0067,  15.9994,  18.9984,    20.183,  22.9898, // C - Na
    24.312,    26.9815,  28.086,  30.9738,    32.064,            // Mg - S
    35.453,    39.948,   39.102,  40.08,      44.956,          // Cl - Sc
    47.9,     50.942,  51.996,   54.938,     55.847,           // Ti - Fe
    58.9332,   58.71,    63.54,   65.37,      69.72,            // Co - Ga
    72.59,     74.9216,  78.96,    79.9009,     83.8,             // Ge - Kr
    85.47,   87.62,    88.905,  91.22,     92.9064,          // Rb - Nb
    95.94,     98,      101.07,   102.906,    106.42,            // Mo - Pd
    107.868,  112.41,   114.82,   118.71,     121.75,            // Ag - Sb
    127.6,    126.905,  131.29,   132.905,    137.33,            // Te - Ba
    138.906,  140.12,   140.908,  144.24,     145,               // La - Pm
    150.36,   151.96,   157.25,   158.925,    162.5,             // Sm - Dy
    164.93,   167.26,   168.934,  173.04,     174.967,           // Ho - Lu
    178.49,   180.948,  183.85,   186.207,    190.2,             // Hf - Os
    192.22,   195.08,   196.967,  200.59,     204.383,           // Ir - Tl
    207.2,    208.98,   209,      210,        222,               // Pb - Rn
    223,      226.025,  227.028,  232.038,    231.036,           // Fr - Pa
    238.029,  237.048,  244,      243,        247,               // U - Cm
    247,      251,      252,      257,        258,               // Bk - Md
    259,      269,      2.0141,   5.45e-4,    -1.0               // No - E
};



static char _symbols[][3] = {
    "H", "He", "Li", "Be", "B",
    "C", "N", "O", "F", "Ne", "Na",
    "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb",
    "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I", "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Ei", "Fm", "Md",
    "No", "Lw", "D", "E", "!"
};



/**
 *  Get table of atomic weights from the internal database.
 *  @param  weights  atomic symbol -> atomic weight map
 *
 *  Example usage:
 *  @code
 *  #include "atomicWeightDB.h"
 *  ...
 *  map<string, double> atw;
 *  getDefaultAtomicWeights(atw);
 *  double copperAtomicWeight = atw["Cu"];
 *  ...
 *  @endcode
 *  Note that if the atomic weight is requested for an unknown
 *  element symbol, the value zero will be returned.
 */
//void getDefaultAtomicWeights(ct::ctmap_sd& weights) {
void getDefaultAtomicWeights(map<string, double>& weights)
{

    // erase existing entries, if any
    //weights.clear();
    const int MAX_NUM = 200;
    int n;
    for (n = 0; n < MAX_NUM; n++) {
        if (_symbols[n][0] == '!') {
            break;
        }
        weights[_symbols[n]] = _weights[n];
    }
}

void writeKnownElements(ostream& s, string fmt)
{

    const int MAX_NUM = 200;
    int n;
    if (fmt == "CK") {
        for (n = 0; n < MAX_NUM; n++) {
            if (_symbols[n][0] == '!') {
                break;
            }
            s << "  " << string(_symbols[n]) << "/" << _weights[n] << "/" << endl;
        }
    } else if (fmt == "XML") {
        s << "<known_elements>" << endl;
        for (n = 0; n < MAX_NUM; n++) {
            if (_symbols[n][0] == '!') {
                break;
            }
            s << "  <element>" << _symbols[n] << "<wt>"
              << _weights[n] << "</wt></element>" << endl;
        }
        s << "</known_elements>" << endl;
    }
}
}

