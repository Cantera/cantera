/**
 * @file import.h
 *
 * functions to import objects from files.
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_IMPORTERS_H
#define CT_IMPORTERS_H

#include "importCK.h"
//#include "importXML.h"
#include "importCTML.h"

namespace Cantera {

    /**
     * import the specifications for a phase from a file, including
     * elements, species, and reactions. The file formats currently
     * supported are
     *
     * - CKML "Chemical Kinetics Markup Language." This XML-based
     * markup language has been developed by M. Aivazis and R. Muller
     * at Caltech and is used in their Fuego software package.
     *
     * - CK This is the name we give to the widely-used format
     * developed by Kee, Miller, and Rupley for the Chemkin-II
     * software package.
     *
     */
    bool importFromFile(Thermo* th, Kinetics* k, map<string, string>& params);

}

#endif
