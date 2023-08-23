/**
 *  @file Elements.cpp
 *  This file contains a database of atomic weights.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Elements.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

/**
 * Database for atomic weights.
 * Values are used from CIAAW. Atomic weights of the elements 2017
 * when a single value is given. Available online at
 * http://www.ciaaw.org/atomic-weights.htm
 *
 * When a range of values is given in the CIAAW table, the "conventional
 * atomic weight" from the IUPAC Periodic Table is used. Available
 * online at https://iupac.org/wp-content/uploads/2018/12/IUPAC_Periodic_Table-01Dec18.pdf
 *
 * If no value is given in either source, it is because no stable isotopes of
 * that element are known and the atomic weight of that element is listed here
 * as -1.0
 *
 * units = kg / kg-mol (or equivalently gm / gm-mol)
 *
 * This structure was picked because it's simple, compact, and extensible.
 */
struct atomicWeightData {
    string symbol; //!< Element symbol, first letter capitalized
    string fullName; //!< Element full name, first letter lowercase
    double atomicWeight; //!< Element atomic weight in kg / kg-mol, if known. -1 if no stable isotope
};

/**
 * Database for named isotopic weights.
 * Values are used from Kim, et al. @cite kim2019.
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

/**
 * @var static vector<atomicWeightData> atomicWeightTable
 * @brief atomicWeightTable is a vector containing the atomic weights database.
 *
 * atomicWeightTable is a static variable with scope limited to this file.
 * It can only be referenced via the functions in this file.
 *
 * The size of the table is given by the initial instantiation.
 */
static vector<atomicWeightData> atomicWeightTable {
    {"H",  "hydrogen",        1.008},
    {"He", "helium",          4.002602},
    {"Li", "lithium",         6.94},
    {"Be", "beryllium",       9.0121831},
    {"B",  "boron",          10.81},
    {"C",  "carbon",         12.011},
    {"N",  "nitrogen",       14.007},
    {"O",  "oxygen",         15.999},
    {"F",  "fluorine",       18.998403163},
    {"Ne", "neon",           20.1797},
    {"Na", "sodium",         22.98976928},
    {"Mg", "magnesium",      24.305},
    {"Al", "aluminum",       26.9815384},
    {"Si", "silicon",        28.085},
    {"P",  "phosphorus",     30.973761998},
    {"S",  "sulfur",         32.06},
    {"Cl", "chlorine",       35.45},
    {"Ar", "argon",          39.95},
    {"K",  "potassium",      39.0983},
    {"Ca", "calcium",        40.078},
    {"Sc", "scandium",       44.955908},
    {"Ti", "titanium",       47.867},
    {"V",  "vanadium",       50.9415},
    {"Cr", "chromium",       51.9961},
    {"Mn", "manganese",      54.938043},
    {"Fe", "iron",           55.845},
    {"Co", "cobalt",         58.933194},
    {"Ni", "nickel",         58.6934},
    {"Cu", "copper",         63.546},
    {"Zn", "zinc",           65.38},
    {"Ga", "gallium",        69.723},
    {"Ge", "germanium",      72.630},
    {"As", "arsenic",        74.921595},
    {"Se", "selenium",       78.971},
    {"Br", "bromine",        79.904},
    {"Kr", "krypton",        83.798},
    {"Rb", "rubidium",       85.4678},
    {"Sr", "strontium",      87.62},
    {"Y",  "yttrium",        88.90584},
    {"Zr", "zirconium",      91.224},
    {"Nb", "nobelium",       92.90637},
    {"Mo", "molybdenum",     95.95},
    {"Tc", "technetium",     -1.0},
    {"Ru", "ruthenium",     101.07},
    {"Rh", "rhodium",       102.90549},
    {"Pd", "palladium",     106.42},
    {"Ag", "silver",        107.8682},
    {"Cd", "cadmium",       112.414},
    {"In", "indium",        114.818},
    {"Sn", "tin",           118.710},
    {"Sb", "antimony",      121.760},
    {"Te", "tellurium",     127.60 },
    {"I",  "iodine",        126.90447},
    {"Xe", "xenon",         131.293},
    {"Cs", "cesium",        132.90545196},
    {"Ba", "barium",        137.327},
    {"La", "lanthanum",     138.90547},
    {"Ce", "cerium",        140.116},
    {"Pr", "praseodymium",  140.90766},
    {"Nd", "neodymium",     144.242},
    {"Pm", "promethium",     -1.0},
    {"Sm", "samarium",      150.36},
    {"Eu", "europium",      151.964},
    {"Gd", "gadolinium",    157.25},
    {"Tb", "terbium",       158.925354},
    {"Dy", "dysprosium",    162.500},
    {"Ho", "holmium",       164.930328},
    {"Er", "erbium",        167.259},
    {"Tm", "thulium",       168.934218},
    {"Yb", "ytterbium",     173.045},
    {"Lu", "lutetium",      174.9668},
    {"Hf", "hafnium",       178.49},
    {"Ta", "tantalum",      180.94788},
    {"W",  "tungsten",      183.84},
    {"Re", "rhenium",       186.207},
    {"Os", "osmium",        190.23 },
    {"Ir", "iridium",       192.217},
    {"Pt", "platinum",      195.084},
    {"Au", "gold",          196.966570},
    {"Hg", "mercury",       200.592},
    {"Tl", "thallium",      204.38},
    {"Pb", "lead",          207.2 },
    {"Bi", "bismuth",       208.98040},
    {"Po", "polonium",       -1.0},
    {"At", "astatine",       -1.0},
    {"Rn", "radon",          -1.0},
    {"Fr", "francium",       -1.0},
    {"Ra", "radium",         -1.0},
    {"Ac", "actinium",       -1.0},
    {"Th", "thorium",       232.0377},
    {"Pa", "protactinium",  231.03588},
    {"U",  "uranium",       238.02891},
    {"Np", "neptunium",      -1.0},
    {"Pu", "plutonium",      -1.0},
    {"Am", "americium",      -1.0},
    {"Cm", "curium",         -1.0},
    {"Bk", "berkelium",      -1.0},
    {"Cf", "californium",    -1.0},
    {"Es", "einsteinium",    -1.0},
    {"Fm", "fermium",        -1.0},
    {"Md", "mendelevium",    -1.0},
    {"No", "nobelium",       -1.0},
    {"Lr", "lawrencium",     -1.0},
    {"Rf", "rutherfordium",  -1.0},
    {"Db", "dubnium",        -1.0},
    {"Sg", "seaborgium",     -1.0},
    {"Bh", "bohrium",        -1.0},
    {"Hs", "hassium",        -1.0},
    {"Mt", "meitnerium",     -1.0},
    {"Ds", "darmstadtium",   -1.0},
    {"Rg", "roentgenium",    -1.0},
    {"Cn", "copernicium",    -1.0},
    {"Nh", "nihonium",       -1.0},
    {"Gl", "flerovium",      -1.0},
    {"Mc", "moscovium",      -1.0},
    {"Lv", "livermorium",    -1.0},
    {"Ts", "tennessine",     -1.0},
    {"Og", "oganesson",      -1.0},
};

/**
 * @var static vector<isotopeWeightData> isotopeWeightTable
 * @brief isotopeWeightTable is a vector containing the atomic weights database.
 *
 * isotopeWeightTable is a static variable with scope limited to this file.
 * It can only be referenced via the functions in this file.
 *
 * The size of the table is given by the initial instantiation.
 */
static vector<isotopeWeightData> isotopeWeightTable {
    // M. Wang et al. The AME2016 atomic mass evaluation. Chinese Physics C.
    // doi:10.1088/1674-1137/41/3/030003.
    {"D",  "deuterium", 2.0141017781, 1},
    {"Tr", "tritium", 3.0160492820, 1},
    {"E", "electron", ElectronMass * Avogadro, 0},
};

// This is implemented as a separate function from elementSymbols() because this pattern
// allows elementSymbols() to return a const reference to the data.
vector<string> elementVectorsFromSymbols() {
    vector<string> values;
    for (const auto& atom : atomicWeightTable) {
        values.push_back(atom.symbol);
    }
    return values;
}

const vector<string>& elementSymbols() {
    const static vector<string> values = elementVectorsFromSymbols();
    return values;
}

// This is implemented as a separate function from elementNames() because this pattern
// allows elementNames() to return a const reference to the data.
vector<string> elementVectorsFromNames() {
    vector<string> values;
    for (const auto& atom : atomicWeightTable) {
        values.push_back(atom.fullName);
    }
    return values;
}

const vector<string>& elementNames() {
    const static vector<string> values = elementVectorsFromNames();
    return values;
}

map<string, double> mapAtomicWeights() {
    map<string, double> symMap;

    for (auto const& atom : atomicWeightTable) {
        symMap.emplace(atom.symbol, atom.atomicWeight);
        symMap.emplace(atom.fullName, atom.atomicWeight);
    }
    for (auto const& isotope : isotopeWeightTable) {
        symMap.emplace(isotope.symbol, isotope.atomicWeight);
        symMap.emplace(isotope.fullName, isotope.atomicWeight);
    }
    return symMap;
}

const map<string, double>& elementWeights() {
    const static map<string, double> symMap = mapAtomicWeights();
    return symMap;
}

double getElementWeight(const string& ename)
{
    const auto& elementMap = elementWeights();
    double elementWeight = 0.0;
    string symbol = trimCopy(ename);
    auto search = elementMap.find(symbol);
    if (search != elementMap.end()) {
        elementWeight = search->second;
    } else {
        string name = toLowerCopy(symbol);
        search = elementMap.find(name);
        if (search != elementMap.end()) {
            elementWeight = search->second;
        }
    }
    if (elementWeight > 0.0) {
        return elementWeight;
    } else if (elementWeight < 0.0) {
        throw CanteraError("getElementWeight",
            "element '{}' has no stable isotopes", ename);
    }
    throw CanteraError("getElementWeight", "element not found: " + ename);
}

double getElementWeight(int atomicNumber)
{
    int num = static_cast<int>(numElementsDefined());
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementWeight", "atomicWeightTable", atomicNumber, num);
    }
    double elementWeight = atomicWeightTable[atomicNumber - 1].atomicWeight;
    if (elementWeight < 0.0) {
        throw CanteraError("getElementWeight",
            "element '{}' has no stable isotopes", getElementName(atomicNumber));
    }
    return elementWeight;
}

string getElementSymbol(const string& ename)
{
    string name = toLowerCopy(trimCopy(ename));
    for (const auto& atom : atomicWeightTable) {
        if (name == atom.fullName) {
            return atom.symbol;
        }
    }
    for (const auto& atom : isotopeWeightTable) {
        if (name == atom.fullName) {
            return atom.symbol;
        }
    }
    throw CanteraError("getElementSymbol", "element not found: " + ename);
}

string getElementSymbol(int atomicNumber)
{
    int num = static_cast<int>(numElementsDefined());
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementSymbol", "atomicWeightTable", atomicNumber, num);
    }
    return atomicWeightTable[atomicNumber - 1].symbol;
}

string getElementName(const string& ename)
{
    string symbol = trimCopy(ename);
    for (const auto& atom : atomicWeightTable) {
        if (symbol == atom.symbol) {
            return atom.fullName;
        }
    }
    for (const auto& atom : isotopeWeightTable) {
        if (symbol == atom.symbol) {
            return atom.fullName;
        }
    }
    throw CanteraError("getElementName", "element not found: " + ename);
}

string getElementName(int atomicNumber)
{
    int num = static_cast<int>(numElementsDefined());
    if (atomicNumber > num || atomicNumber < 1) {
        throw IndexError("getElementName", "atomicWeightTable", atomicNumber, num);
    }
    return atomicWeightTable[atomicNumber - 1].fullName;
}

int getAtomicNumber(const string& ename)
{
    size_t numElements = numElementsDefined();
    size_t numIsotopes = numIsotopesDefined();
    string symbol = trimCopy(ename);
    string name = toLowerCopy(symbol);
    for (size_t i = 0; i < numElements; i++) {
        if (symbol == atomicWeightTable[i].symbol) {
            return static_cast<int>(i) + 1;
        } else if (name == atomicWeightTable[i].fullName) {
            return static_cast<int>(i) + 1;
        }
    }
    for (size_t i = 0; i < numIsotopes; i++) {
        if (symbol == isotopeWeightTable[i].symbol) {
            return isotopeWeightTable[i].atomicNumber;
        } else if (name == isotopeWeightTable[i].fullName) {
            return isotopeWeightTable[i].atomicNumber;
        }
    }
    throw CanteraError("getAtomicNumber", "element not found: " + ename);
}

size_t numElementsDefined()
{
    return atomicWeightTable.size();
}

size_t numIsotopesDefined()
{
    return isotopeWeightTable.size();
}

}
