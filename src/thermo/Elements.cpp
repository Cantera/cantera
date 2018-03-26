/**
 *  @file Elements.cpp
 *  This file contains a database of atomic weights.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Elements.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

/*! Database for atomic weights
 * Values are taken from the 1989 Standard Atomic Weights, CRC
 *
 * units = kg / kg-mol (or equivalently gm / gm-mol)
 *
 * This structure was picked because it's simple, compact, and extensible.
 */
struct atomicWeightData {
    string symbol; //!< Element symbol, first letter capitalized
    string fullName; //!< Element full name, first letter lowercase
    double atomicWeight; //!< Element atomic weight in kg / kg-mol
};

/*! Database for named isotopic weights
 * Values are taken from the 1989 Standard Atomic Weights, CRC
 *
 * units = kg / kg-mol (or equivalently gm / gm-mol)
 *
 * This structure was picked because it's simple, compact, and extensible.
 */
struct isotopeWeightData {
    string symbol; //!< Isotope symbol, first letter capitalized
    string fullName; //!< Isotope full name, first letter lowercase
    double atomicWeight; //!< Isotope atomic weight in kg / kg-mol
    int atomicNumber; //!< Isotope atomic number
};

/*!
 * @var static struct atomicWeightData atomicWeightTable[]
 * \brief atomicWeightTable is a vector containing the atomic weights database.
 *
 * atomicWeightTable[] is a static function with scope limited to this file.
 * It can only be referenced via the functions in this file.
 *
 * The size of the table is given by the initial instantiation.
 */
static struct atomicWeightData atomicWeightTable[] = {
    {"H",  "hydrogen",       1.00794},
    {"He", "helium",         4.002602},
    {"Li", "lithium",        6.941  },
    {"Be", "beryllium",      9.012182},
    {"B",  "boron",         10.811  },
    {"C",  "carbon",        12.011  },
    {"N",  "nitrogen",      14.00674},
    {"O",  "oxygen",        15.9994 },
    {"F",  "fluorine",      18.9984032},
    {"Ne", "neon",          20.1797 },
    {"Na", "sodium",        22.98977},
    {"Mg", "magnesium",     24.3050 },
    {"Al", "aluminum",      26.98154},
    {"Si", "silicon",       28.0855 },
    {"P",  "phosphorus",    30.97376},
    {"S",  "sulfur",        32.066  },
    {"Cl", "chlorine",      35.4527 },
    {"Ar", "argon",         39.948  },
    {"K",  "potassium",     39.0983 },
    {"Ca", "calcium",       40.078  },
    {"Sc", "scandium",      44.95591},
    {"Ti", "titanium",      47.88   },
    {"V",  "vanadium",      50.9415 },
    {"Cr", "chromium",      51.9961 },
    {"Mn", "manganese",     54.9381 },
    {"Fe", "iron",          55.847  },
    {"Co", "cobalt",        58.9332 },
    {"Ni", "nickel",        58.69   },
    {"Cu", "copper",        63.546  },
    {"Zn", "zinc",          65.39   },
    {"Ga", "gallium",       69.723  },
    {"Ge", "germanium",     72.61   },
    {"As", "arsenic",       74.92159},
    {"Se", "selenium",      78.96   },
    {"Br", "bromine",       79.904  },
    {"Kr", "krypton",       83.80   },
    {"Rb", "rubidium",      85.4678 },
    {"Sr", "strontium",     87.62   },
    {"Y",  "yttrium",       88.90585},
    {"Zr", "zirconium",     91.224  },
    {"Nb", "nobelium",      92.90638},
    {"Mo", "molybdenum",    95.94   },
    {"Tc", "technetium",    97.9072 },
    {"Ru", "ruthenium",    101.07   },
    {"Rh", "rhodium",      102.9055 },
    {"Pd", "palladium",    106.42   },
    {"Ag", "silver",       107.8682 },
    {"Cd", "cadmium",      112.411  },
    {"In", "indium",       114.82   },
    {"Sn", "tin",          118.710  },
    {"Sb", "antimony",     121.75   },
    {"Te", "tellurium",    127.6    },
    {"I",  "iodine",       126.90447},
    {"Xe", "xenon",        131.29   },
    {"Cs", "cesium",       132.90543},
    {"Ba", "barium",       137.327  },
    {"La", "lanthanum",    138.9055 },
    {"Ce", "cerium",       140.115  },
    {"Pr", "praseodymium", 140.90765},
    {"Nd", "neodymium",    144.24   },
    {"Pm", "promethium",   144.9127 },
    {"Sm", "samarium",     150.36   },
    {"Eu", "europium",     151.965  },
    {"Gd", "gadolinium",   157.25   },
    {"Tb", "terbium",      158.92534},
    {"Dy", "dysprosium",   162.50   },
    {"Ho", "holmium",      164.93032},
    {"Er", "erbium",       167.26   },
    {"Tm", "thulium",      168.93421},
    {"Yb", "ytterbium",    173.04   },
    {"Lu", "lutetium",     174.967  },
    {"Hf", "hafnium",      178.49   },
    {"Ta", "tantalum",     180.9479 },
    {"W",  "tungsten",     183.85   },
    {"Re", "rhenium",      186.207  },
    {"Os", "osmium",       190.2    },
    {"Ir", "iridium",      192.22   },
    {"Pt", "platinum",     195.08   },
    {"Au", "gold",         196.96654},
    {"Hg", "mercury",      200.59   },
    {"Tl", "thallium",     204.3833 },
    {"Pb", "lead",         207.2    },
    {"Bi", "bismuth",      208.98037},
    {"Po", "polonium",     208.9824 },
    {"At", "astatine",     209.9871 },
    {"Rn", "radon",        222.0176 },
    {"Fr", "francium",     223.0197 },
    {"Ra", "radium",       226.0254 },
    {"Ac", "actinium",     227.0279 },
    {"Th", "thorium",      232.0381 },
    {"Pa", "protactinium", 231.03588},
    {"U",  "uranium",      238.0508 },
    {"Np", "neptunium",    237.0482 },
    {"Pu", "plutonium",    244.0482 },
};

/*!
 * @var static struct isotopeWeightData isotopeWeightTable[]
 * \brief isotopeWeightTable is a vector containing the atomic weights database.
 *
 * isotopeWeightTable[] is a static function with scope limited to this file.
 * It can only be referenced via the functions in this file.
 *
 * The size of the table is given by the initial instantiation.
 */
static struct isotopeWeightData isotopeWeightTable[] = {
    {"D",  "deuterium", 2.0, 1},
    {"Tr", "tritium",   3.0, 1},
    {"E", "electron", 0.000545, 0},
};

double getElementWeight(const std::string& ename)
{
    int numElements = numElementsDefined();
    int numIsotopes = numIsotopesDefined();
    string symbol = trimCopy(ename);
    string name = toLowerCopy(symbol);
    for (int i = 0; i < numElements; i++) {
        if (symbol == atomicWeightTable[i].symbol) {
            return atomicWeightTable[i].atomicWeight;
        } else if (name == atomicWeightTable[i].fullName) {
            return atomicWeightTable[i].atomicWeight;
        }
    }
    for (int i = 0; i < numIsotopes; i++) {
        if (symbol == isotopeWeightTable[i].symbol) {
            return isotopeWeightTable[i].atomicWeight;
        } else if (name == isotopeWeightTable[i].fullName) {
            return isotopeWeightTable[i].atomicWeight;
        }
    }
    throw CanteraError("getElementWeight", "element not found: " + ename);
}

double getElementWeight(int atomicNumber)
{
    int num = numElementsDefined();
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementWeight", "atomicWeightTable", atomicNumber, num);
    }
    return atomicWeightTable[atomicNumber - 1].atomicWeight;
}

string getElementSymbol(const std::string& ename)
{
    int numElements = numElementsDefined();
    int numIsotopes = numIsotopesDefined();
    string name = toLowerCopy(trimCopy(ename));
    for (int i = 0; i < numElements; i++) {
        if (name == atomicWeightTable[i].fullName) {
            return atomicWeightTable[i].symbol;
        }
    }
    for (int i = 0; i < numIsotopes; i++) {
        if (name == isotopeWeightTable[i].fullName) {
            return isotopeWeightTable[i].symbol;
        }
    }
    throw CanteraError("getElementSymbol", "element not found: " + ename);
}

string getElementSymbol(int atomicNumber)
{
    int num = numElementsDefined();
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementSymbol", "atomicWeightTable", atomicNumber,
                         num);
    }
    return atomicWeightTable[atomicNumber - 1].symbol;
}

string getElementName(const std::string& ename)
{
    int numElements = numElementsDefined();
    int numIsotopes = numIsotopesDefined();
    string symbol = trimCopy(ename);
    for (int i = 0; i < numElements; i++) {
        if (symbol == atomicWeightTable[i].symbol) {
            return atomicWeightTable[i].fullName;
        }
    }
    for (int i = 0; i < numIsotopes; i++) {
        if (symbol == isotopeWeightTable[i].symbol) {
            return isotopeWeightTable[i].fullName;
        }
    }
    throw CanteraError("getElementName", "element not found: " + ename);
}

string getElementName(int atomicNumber)
{
    int num = numElementsDefined();
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementName", "atomicWeightTable", atomicNumber,
                         num);
    }
    return atomicWeightTable[atomicNumber - 1].fullName;
}

int getAtomicNumber(const std::string& ename)
{
    int numElements = numElementsDefined();
    int numIsotopes = numIsotopesDefined();
    string symbol = trimCopy(ename);
    string name = toLowerCopy(symbol);
    for (int i = 0; i < numElements; i++) {
        if (symbol == atomicWeightTable[i].symbol) {
            return i+1;
        } else if (name == atomicWeightTable[i].fullName) {
            return i+1;
        }
    }
    for (int i = 0; i < numIsotopes; i++) {
        if (symbol == isotopeWeightTable[i].symbol) {
            return isotopeWeightTable[i].atomicNumber;
        } else if (name == isotopeWeightTable[i].fullName) {
            return isotopeWeightTable[i].atomicNumber;
        }
    }
    throw CanteraError("getAtomicNumber", "element not found: " + ename);
}

int numElementsDefined()
{
    return sizeof(atomicWeightTable) / sizeof(struct atomicWeightData);
}

int numIsotopesDefined()
{
    return sizeof(isotopeWeightTable) / sizeof(struct isotopeWeightData);
}

}
