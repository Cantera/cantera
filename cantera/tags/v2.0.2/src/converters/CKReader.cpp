/**
 *  @file CKReader.cpp
 *
 */

// Copyright 2001  California Institute of Technology

#include <fstream>
#include <string>
using namespace std;

#include "CKParser.h"
#include "CKReader.h"
#include "thermoFunctions.h"

#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#include "writelog.h"
#include <cstdio>
#include "ckr_defs.h"

//#include "cantera/base/global.h"
//#define APP Cantera::Application

namespace ckr
{

/**
 *  read and optionally validate an input file in Chemkin format.
 *  @param inputFile   path to the input file
 *  @param thermoDatabase  path to the species database file
 *  @param logfile  path to the file where log messages should be written
 *  @return true if no errors were encountered, false otherwise
 */
bool CKReader::read(const std::string& inputFile, const std::string& thermoDatabase,
                    const std::string& logfile)
{

    clock_t t0, t1;

    t0 = clock();

    ifstream ckinfile(inputFile.c_str());
    ofstream log(logfile.c_str());
    try {
        // construct a parser for the input file
        CKParser parser(&ckinfile, inputFile, &log);
        parser.verbose = verbose;
        parser.debug = debug;

        // write header information to the log file
        struct tm* newtime;
        time_t aclock;
        time(&aclock);                    /* Get time in seconds */
        newtime = localtime(&aclock);     /* Convert time to struct tm form */


        log << "CKReader version 1.0" << endl
            << "http://www.cantera.org" << endl << endl
            << asctime(newtime) << endl
            << setw(20) << "input file: "
            << setw(30) << inputFile << endl;

        if (thermoDatabase != "")
            log << setw(20) << "species database: "
                << setw(30) << thermoDatabase << endl;

        if (!validate)
            log << endl << "***************  Warning  ***************" << endl
                <<         "     mechanism validation disabled"        << endl
                <<         "*****************************************" << endl;

        if (debug) {
            log << "*** DEBUG MODE ***" << endl;
        } else {
            log << "debugging disabled." << endl;
        }

        //-----------   process ELEMENT section  ----------------------

        bool elok = parser.readElementSection(elements);
        int nel = static_cast<int>(elements.size());
        vector<string> elementSymbols;
        for (int j = 0; j < nel; j++) {
            elementSymbols.push_back(elements[j].name);
        }

        if (verbose) {
            log.flags(ios::showpoint);
            log.precision(6);
            log.width(0);

            log << endl << newTask("reading elements") << endl;

            // write summary to log file
            for (int i = 0; i < nel; i++) {
                log << i+1 << ".  " << pad(elements[i].name,2) << "  ";
                double wt = elements[i].atomicWeight;
                if (wt == 0.0) {
                    log << "<error>";
                } else {
                    log << wt;
                }
                if (!elements[i].weightFromDB) {
                    log << " (specified)";
                }
                if (elements[i].comment != "") {
                    log << "      ! " << elements[i].comment;
                }
                log << endl;
            }
        }
        log << "\nread " << nel << " elements." << endl;

        if (!elok) {
            log << "\nerrors were encountered." << endl;
            return false;
        }
        if (nel == 0) {
            return false;
        }


        //------------   process SPECIES section  ------------------------

        vector<string> speciesSymbols;
        bool spok = parser.readSpeciesSection(species);
        int nsp = static_cast<int>(species.size());

        if (verbose) {
            log << newTask("reading species") << endl;
        }

        for (int i = 0; i < nsp; i++) {
            Species& s = species[i];
            if (verbose) {
                log << i+1 << ".  " << s.name << endl;
            }
            speciesSymbols.push_back(s.name);
        }
        log << "\nread " << nsp << " species." << endl;

        if (!spok) {
            log << "\nerrors were encountered." << endl;
            return false;
        }
        if (nsp == 0) {
            return false;
        }


        //-------------   process THERMO section  -------------------------

        if (verbose) {
            log << newTask("looking up species definitions") << endl;
        }

        // if a thermo database is specified, get the default Tmin, Tmid, Tmax
        vector_fp temp;

        if (thermoDatabase != "") {

            if (verbose) log << "reading default temperature ranges from "
                                 << thermoDatabase  << endl;

            ifstream thermofile(thermoDatabase.c_str());
            CKParser thermReader(&thermofile, thermoDatabase, &log);
            thermReader.verbose = verbose;
            thermReader.debug = debug;
            int dbflag = HasTempRange;
            vector<string> dummy;
            thermReader.readThermoSection(dummy, speciesData, temp, dbflag, log);
        }


        bool hasthermo = parser.advanceToKeyword("THERM","REAC");

        int k, optionFlag = 0;
        int undefined = static_cast<int>(species.size());
        string nm;
        vector<string> undef;
        bool allsp = (speciesSymbols[0] == "<ALL>");
        if (hasthermo &&
                parser.readThermoSection(speciesSymbols,
                                         speciesData, temp, optionFlag, log)) {
            if (allsp) {
                nsp = static_cast<int>(speciesData.size());
                for (k = 0; k < nsp; k++) {
                    Species s;
                    s.name = speciesSymbols[k];
                    species.push_back(s);
                }
            }
            undefined = 0;
            for (k = 0; k < nsp; k++) {
                nm = species[k].name;
                species[k] = speciesData[species[k].name];
                if (species[k].name == "<empty>") {
                    undefined++;
                    undef.push_back(nm);
                    species[k].name = nm;
                }
            }
            int localdefs = nsp - undefined;
            if (localdefs > 0 && verbose) log << "found definitions for "
                                                  << localdefs
                                                  << " of "
                                                  << nsp
                                                  << " species in the input file. "
                                                  << endl;
        } else {
            undef = speciesSymbols;
            if (verbose) {
                log << "no THERMO section in input file." << endl;
            }
        }

        if (undefined > 0 && thermoDatabase != ""
                && optionFlag != NoThermoDatabase) {

            if (verbose) log << "searching external database "
                                 << thermoDatabase << " for species definitions..."
                                 << endl;

            ifstream thermofile(thermoDatabase.c_str());
            CKParser thermoReader(&thermofile, thermoDatabase, &log);
            thermoReader.verbose = verbose;
            thermoReader.debug = debug;
            int dbflag = HasTempRange;
            thermoReader.readThermoSection(undef, speciesData, temp, dbflag, log);
            undefined = 0;
            if (allsp) {
                species.clear();
                nsp = static_cast<int>(speciesData.size());
                for (k = 0; k < nsp; k++) {
                    Species s;
                    s.name = undef[k];
                    species.push_back(s);
                }
            }
            for (int k = 0; k < nsp; k++) {
                if (species[k].valid == 0) {
                    nm = species[k].name;
                    species[k] = speciesData[species[k].name];
                    if (species[k].name == "<empty>") {
                        undefined++;
                        species[k].name = nm;
                    }
                }
            }
        }

        if (validate && !validateSpecies(log)) {
            //Cantera::setError("read","error in species");
            return false;
        }


        //-------------   process REACTIONS section  -------------------------

        if (verbose) {
            log << newTask("reading reactions") << endl;
        }

        ckinfile.close();
        ifstream ckinfile2(inputFile.c_str());

        // construct a new parser for the input file
        CKParser parser2(&ckinfile2, inputFile, &log);
        parser2.verbose = verbose;
        parser2.debug = debug;

        parser2.readReactionSection(speciesSymbols, elementSymbols, reactions, units);
        log << "\nread " << static_cast<int>(reactions.size())
            << " reactions." << endl;

        bool rxnok = true;
        if (validate) {
            rxnok = rxnok && validateReactions(log);
        }
        bool writeok = true;
        if (verbose || validate) {
            writeok = writeReactions(log);
        }
        rxnok = rxnok && writeok;
        if (!rxnok) {
            return false;
        }

        log << "\nSuccess... ";
        t1 = clock();
        log << "elapsed CPU time = "
            << double(t1 - t0)/CLOCKS_PER_SEC
            << " s" << endl;
        if (!validate) {
            log << "*** no validation performed ***" << endl;
        }
    }

    catch (CK_Exception& e) {
        log << e.errorMessage() << endl;
        //Cantera::setError("CKReader::read",e.errorMessage());
        return false;
    } catch (std::exception& e) {
        log << "an exception was raised in CKReader:\n";
        log << e.what() << std::endl;
        return false;
    }

    return true;
}


/// print a summary of all reactions to the log file
bool CKReader::writeReactions(std::ostream& log)
{

    bool ok = true;
    //    int ns = species.size();
    int nrxns = static_cast<int>(reactions.size());
    log.flags(ios::unitbuf);
    log.precision(6);

    log << endl;
    for (int n = 0; n < nrxns; n++) {
        Reaction& r = reactions[n];

        log << "reaction " << r.number << endl;
        log << "   ";
        printReactionEquation(log, r);
        log << endl;

        // rate coefficient
        if (r.isFalloffRxn) {

            log << "   high P rate coeff: ";
            ok = ok && writeRateCoeff(r.kf, log) ;

            log << "   low P rate coeff: ";
            ok = ok && writeRateCoeff(r.kf_aux, log);

            ok = ok && writeFalloff(r.falloffType, r.falloffParameters, log);

        } else {
            log << "   rate coeff: ";
            ok = ok && writeRateCoeff(r.kf, log);
        }
        if (r.isReversible && r.krev.A > 0) {
            log << "   reverse rate coeff: ";
            ok = ok && writeRateCoeff(r.krev, log);
        }
        int ne = static_cast<int>(r.e3b.size());

        if (ne > 0) {
            vector<string> enhSpecies;
            getMapKeys(r.e3b, enhSpecies);
            log << "   enhanced collision efficiencies:" << endl;
            log << "       ";
            for (int nn = 0; nn < ne; nn++) {
                log << enhSpecies[nn] << " "
                    << r.e3b[enhSpecies[nn]];
                if (nn < ne-1) {
                    log << ",  ";
                }
            }
            log << endl;
        }
        if (r.isDuplicate) log
                    << "   declared duplicate reaction. See reaction "
                    << r.duplicate << "." << endl;
        log << endl;
    }
    return ok;
}



/// validate the species
bool CKReader::validateSpecies(std::ostream& log)
{
    size_t nel = elements.size();
    size_t nsp = species.size();
    double tol;

    log << newTask("validating species");

    // check for undeclared elements
    vector<string> esyms;

    log << "   checking that all species have been defined... ";
    for (size_t k = 0; k < nsp; k++) {
        Species& s = species[k];
        if (s.valid == 0) {
            log << endl << "   species " << s.name << " undefined ";
            s.valid = -1;
        }
    }
    if (valid(species)) {
        log << "OK" << endl;
    } else {
        log << endl;
        return false;
    }

    log << "   checking that all species elements have been declared... ";
    for (size_t k = 0; k < nsp; k++) {
        Species& s = species[k];

        getMapKeys(s.comp, esyms);
        size_t nm = esyms.size();

        for (size_t m = 0; m < nm; m++) {
            size_t j;
            for (j = 0; j < nel; j++) {
                if (esyms[m] == elements[j].name) {
                    break;
                }
            }
            if (j == nel) {
                log << endl << "   species "
                    << s.name << ": undeclared element "
                    << esyms[m];
                s.valid = -1;
            }
        }
    }
    if (valid(species)) {
        log << "OK" << endl;
    } else {
        log << endl;
        return false;
    }

    log << "   checking consistency of species thermo data... ";
    tol = 0.01;
    if (checkThermo(log, species, tol)) {
        log << "OK" << endl;
    } else {
        log << endl;
        return false;
    }
    return true;

}


/// validate the reactions
bool CKReader::validateReactions(std::ostream& log)
{

    bool ok = true;
    //    int ns = species.size();
    int nrxns = static_cast<int>(reactions.size());

    vector<int> unbal;
    log << "checking that all reactions balance...";
    if (checkBalance(log, speciesData, reactions, unbal)) {
        log << " OK" << endl;
    } else {
        int nu = static_cast<int>(unbal.size());
        for (int iu = 0; iu < nu; iu++) {
            log << "   error... reaction " << unbal[iu]
                << " does not balance" << endl;
        }
        ok = false;
    }

    log << "checking for duplicate reactions...";

    for (int nn = 0; nn < nrxns; nn++) {
        Reaction& r1 = reactions[nn];
        for (int mm = nn + 1; mm < nrxns; mm++) {
            Reaction& r2 = reactions[mm];
            if (r1 == r2) {
                r1.duplicate = mm + 1;
                r2.duplicate = nn + 1;
                if (!r1.isDuplicate || !r2.isDuplicate) {
                    log << endl << "   error... undeclared duplicate reactions: "
                        << nn + 1 << " and " << mm + 1;
                    ok = false;
                } else {
                    log << endl << "   declared duplicate reactions: "
                        << nn + 1
                        << " and " << mm + 1;
                }
            }
        }
    }
    if (ok) {
        log << "...OK" << endl;
    }

    return ok;
}


/**
 *  Check the thermodynamic property parameterizations for all species.
 *  The following are verified:
 *  - The heat capacity is positive throughout the full temperature range;
 *  - The entropy at Tmin is positive;
 *  - The heat capacity, entropy, and enthalpy evaluated at Tmid using
 *    both the high and low polynomial coefficients are the same to within
 *    relative error tol
 *  - The heat capacity at Tmax is not greater than the equipartition limit
 *    for the number of atoms in the molecule
 */
bool checkThermo(std::ostream& log, speciesList& sp, double tol)
{
    const double dt = 0.0001;
    double t, cp0, h0, s0, cp1, h1, s1;
    int nsp = static_cast<int>(sp.size());
    const int n_points = 20;

    int k;
    bool ok = true;
    for (k = 0; k < nsp; k++) {
        Species& s = sp[k];

        if (s.valid <= 0) {
            ok = false;
            log << endl << "species " << s.name
                << " contains an error." << endl;
        }
        if (!ok) {
            return false;
        }
    }

    log << endl << "   Checking that cp/R is positive... ";

    for (k = 0; k < nsp; k++) {
        Species& s = sp[k];

        // check that cp is positive

        cp0 = 0.0;
        for (int j = 0; j < n_points; j++) {
            t = j*(s.thigh - s.tlow)/(n_points - 1) + s.tlow;

            cp0 = cp(t, s);
            if (cp0 < 0.0) {
                log << endl << "   error... Cp/R < 0 at T = " << t
                    << " for species " << s.name << endl;
                s.valid = -1;
                ok = false;
            }
        }
    }
    if (ok) {
        log << "ok" << endl;
    } else {
        return ok;
    }


    // check that S(tlow) > 0

    log << "   Checking that the species entropies are positive... ";
    for (k = 0; k < nsp; k++) {
        Species& s = sp[k];
        if (entropy(s.tlow, s) <= 0.0) {
            log << endl << "   error... negative entropy for species "
                << s.name << endl;
            s.valid = -1;
            ok = false;
        }
    }
    if (ok) {
        log << "ok" << endl;
    } else {
        return ok;
    }

    log << "   Checking that properties are continuous at the midpoint temperature... ";
    for (k = 0; k < nsp; k++) {
        Species& s = sp[k];

        // check continuity at Tmid
        t = s.tmid - dt;
        cp0 = cp(t, s);
        h0 = enthalpy(t, s) + cp0*dt;
        s0 = entropy(t, s) + dt*cp0/t;

        t = s.tmid + dt;
        cp1 = cp(t, s);
        h1 = enthalpy(t, s) - cp1*dt;
        s1 = entropy(t, s) - cp1*dt/t;

        if (absval((cp0 - cp1)/cp0) > tol) {
            log << endl << "Warning... species " << s.name
                << ": discontinuity in Cp at Tmid = "
                << s.tmid << endl;
            log << "Cp/R (low, high) = " << cp0
                << ", " << cp1 << endl;
            ok = false;
        }

        if (absval((h0 - h1)/h0) > tol) {
            log << endl << "Warning... species " << s.name
                << ": discontinuity in enthalpy at Tmid = "
                << s.tmid << endl;
            log << "h/R (low, high) = "
                << h0 << ", " << h1 << endl;
            ok = false;
        }

        if (absval((s0 - s1)/s0) > tol) {
            log << endl << "Warning... species " << s.name
                << ": discontinuity in entropy at Tmid = "
                << s.tmid << endl;
            log << "s/R (low, high) = "
                << s0 << ", " << s1 << endl;
            ok = false;
        }
    }
    if (ok) {
        log << "ok \n\n\n";
    } else {
        log << "\n\n\n";
    }

    log << "   Checking that cp is less that the high-temperature\n"
        << "   limiting value predicted by equipartition of energy.\n";
    log << "   Note that this limit does not account for the electronic\n"
        << "   contribution to cp, and so may be violated in some cases."
        << endl << endl;


    for (k = 0; k < nsp; k++) {
        Species& s = sp[k];

        // check that cp at Tmax is less than the equipartion value
        // This does not include any possible electronic contribution

        cp0 = cp(s.thigh, s);
        int nel = static_cast<int>(s.elements.size());
        int i;
        double na = 0.0;
        for (i = 0; i < nel; i++)
            if (s.elements[i].name != "E") {
                na += s.elements[i].number;
            }
        int natoms = int(floor(na));
        double cpmax;
        switch (natoms) {
        case 1:
            cpmax = 2.5;
            break;
        case 2:
            cpmax = 4.5;
            break;
        default:
            cpmax = 3.0*natoms - 3.0;
        }

        if (cp0 > cpmax) {
            double over = 100.0*(cp0 - cpmax)/cpmax;
            log << endl << "Warning... species " << s.name
                << ": cp(Tmax) greater than equipartition value \nby "
                <<  over << " percent.";
            if ((natoms > 2) && (cp0 - cpmax < 0.5))
                log << endl << "      (if molecule is linear, cp is ok)"
                    << endl;
        }
    }

    return valid(sp);
}



/**
 * Check that all reactions balance.
 * @param speciesData species property dataset used to look up
 * elemental compositions.
 * @param r list of reactions to check
 * @param unbalanced list of integers specifying reaction numbers of
 * unbalanced reactions.
 * @return true if all reactions balance
 * @todo use reaction number stored in reaction object
 */
bool checkBalance(std::ostream& f, speciesTable& speciesData,
                  reactionList& r, std::vector<int>& unbalanced, double tolerance)
{
    int nrxn = static_cast<int>(r.size());
    string rname, pname;
    vector<string> elementSyms;
    unsigned int m;

    unbalanced.clear();
    map<string, double> atoms;

    for (int i = 0; i < nrxn; i++) {
        atoms.clear();
        int nr = static_cast<int>(r[i].reactants.size());
        int np = static_cast<int>(r[i].products.size());
        int j;
        double stoichCoeff;
        for (j = 0; j < nr; j++) {
            rname = r[i].reactants[j].name;
            stoichCoeff = r[i].reactants[j].number;
            vector<Constituent>& elements = speciesData[rname].elements;
            for (m = 0; m < elements.size(); m++) {
                atoms[elements[m].name] -= stoichCoeff * elements[m].number;
            }
        }
        for (j = 0; j < np; j++) {
            pname = r[i].products[j].name;
            stoichCoeff = r[i].products[j].number;
            vector<Constituent>& elements = speciesData[pname].elements;
            for (m = 0; m < elements.size(); m++) {
                atoms[elements[m].name] += stoichCoeff * elements[m].number;
            }
        }
        double atms;
        getMapKeys(atoms, elementSyms);
        for (m = 0; m < elementSyms.size(); m++) {
            atms = atoms[elementSyms[m]];
            if (fabs(atms) > tolerance) {
                //if (atoms[elementSyms[m]] != 0.0) {
                //                    cout << "Reaction " << i+1 << " has an unbalanced element: "
                //   << elementSyms[m] << "  "
                //   << atoms[elementSyms[m]] << endl;
                unbalanced.push_back(i+1);
                break;
            }
        }
    }
    return (unbalanced.empty());
}

}


