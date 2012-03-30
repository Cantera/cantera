/**
 *  @file Elements.cpp
 *  This file contains a database of atomic weights.
 */
//  Copyright 2003  California Institute of Technology

#include "cantera/thermo/Elements.h"
#include "cantera/base/ctexceptions.h"

using namespace std;

#include <cstdlib>

namespace Cantera
{

/*! Database for atomic molecular weights
 *  Values are taken from the 1989 Standard Atomic Weights, CRC
 *
 *  awTable[] is a static function with scope limited to this file.
 *  It can only be referenced via the LookupWtElements() function.
 *
 *  units = kg / kg-mol (or equivalently gm / gm-mol)
 *
 * This structure was picked because it's simple, compact, and extensible.
 */
struct awData {
    char name[4];         //!< Null Terminated name, First letter capitalized
    double atomicWeight;  //!< atomic weight in kg / kg-mol
};

/*!
 *  @var static struct awData aWTable[]
 *  \brief aWTable is a vector containing the atomic weights database.
 *
 *  The size of the table is given by the initial instantiation.
 */
static struct awData aWTable[] = {
    {"H",    1.00794},
    {"D",    2.0    },
    {"Tr",   3.0    },
    {"He",   4.002602},
    {"Li",   6.941  },
    {"Be",   9.012182},
    {"B",   10.811  },
    {"C",   12.011  },
    {"N",   14.00674},
    {"O",   15.9994 },
    {"F",   18.9984032},
    {"Ne",  20.1797 },
    {"Na",  22.98977},
    {"Mg",  24.3050 },
    {"Al",  26.98154},
    {"Si",  28.0855 },
    {"P",   30.97376},
    {"S",   32.066  },
    {"Cl",  35.4527 },
    {"Ar",  39.948  },
    {"K",   39.0983 },
    {"Ca",  40.078  },
    {"Sc",  44.95591},
    {"Ti",  47.88   },
    {"V",   50.9415 },
    {"Cr",  51.9961 },
    {"Mn",  54.9381 },
    {"Fe",  55.847  },
    {"Co",  58.9332 },
    {"Ni",  58.69   },
    {"Cu",  63.546  },
    {"Zn",  65.39   },
    {"Ga",  69.723  },
    {"Ge",  72.61   },
    {"As",  74.92159},
    {"Se",  78.96   },
    {"Br",  79.904  },
    {"Kr",  83.80   },
    {"Rb",  85.4678 },
    {"Sr",  87.62   },
    {"Y",   88.90585},
    {"Zr",  91.224  },
    {"Nb",  92.90638},
    {"Mo",  95.94   },
    {"Tc",  97.9072 },
    {"Ru", 101.07   },
    {"Rh", 102.9055 },
    {"Pd", 106.42   },
    {"Ag", 107.8682 },
    {"Cd", 112.411  },
    {"In", 114.82   },
    {"Sn", 118.710  },
    {"Sb", 121.75   },
    {"Te", 127.6    },
    {"I",  126.90447},
    {"Xe", 131.29   },
    {"Cs", 132.90543},
    {"Ba", 137.327  },
    {"La", 138.9055 },
    {"Ce", 140.115  },
    {"Pr", 140.90765},
    {"Nd", 144.24   },
    {"Pm", 144.9127 },
    {"Sm", 150.36   },
    {"Eu", 151.965  },
    {"Gd", 157.25   },
    {"Tb", 158.92534},
    {"Dy", 162.50   },
    {"Ho", 164.93032},
    {"Er", 167.26   },
    {"Tm", 168.93421},
    {"Yb", 173.04   },
    {"Lu", 174.967  },
    {"Hf", 178.49   },
    {"Ta", 180.9479 },
    {"W",  183.85   },
    {"Re", 186.207  },
    {"Os", 190.2    },
    {"Ir", 192.22   },
    {"Pt", 195.08   },
    {"Au", 196.96654},
    {"Hg", 200.59   },
    {"Ti", 204.3833 },
    {"Pb", 207.2    },
    {"Bi", 208.98037},
    {"Po", 208.9824 },
    {"At", 209.9871 },
    {"Rn", 222.0176 },
    {"Fr", 223.0197 },
    {"Ra", 226.0254 },
    {"Ac", 227.0279 },
    {"Th", 232.0381 },
    {"Pa", 231.03588},
    {"U",  238.0508 },
    {"Np", 237.0482 },
    {"Pu", 244.0482 }
};

doublereal LookupWtElements(const std::string& ename)
{
    int num = sizeof(aWTable) / sizeof(struct awData);
    string s3 = ename.substr(0,3);
    for (int i = 0; i < num; i++) {
        if (s3 == aWTable[i].name) {
            return (aWTable[i].atomicWeight);
        }
    }
    throw CanteraError("LookupWtElements", "element not found");
    return -1.0;
}

}
