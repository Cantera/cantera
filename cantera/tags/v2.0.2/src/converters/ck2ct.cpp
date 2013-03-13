/**
 * @file ck2ct.cpp
 *   Convert CK-format reaction mechanism files to Cantera input format.
 */
#include <iostream>
#include <string>
#include <ctype.h>

#include "cantera/base/config.h"

#ifdef HAS_SSTREAM
#include <sstream>
#endif
using namespace std;

#include "CKReader.h"
#include "Reaction.h"
#include "writelog.h"

#include "ck2ct.h"
#include "cantera/base/stringUtils.h"
#include <time.h>
#include "cantera/base/ct_defs.h"
#include "cantera/base/ctml.h"


using namespace Cantera;

namespace pip
{

struct trdata {
    int geom;
    doublereal welldepth, diam, dipole, polar, rot;
};

static map<string, trdata> _trmap;
static bool _with_transport = false;

static void getTransportData(string trfile)
{

    _with_transport = true;
    ifstream s(trfile.c_str());
    if (!s) throw CanteraError("getTransportData",
                                   "could not open transport database "+trfile);

    /*
     *  The first thing we will do is to read the entire transport
     *  database and place its contents into a map structure,
     *  indexed by the name of the species.
     */
    string rest;
    while (! s.eof()) {
        /*
         * HKM Note: the USE_STRINGSTREAM block works for files
         * with comments in them. The other block gets hung up
         * somehow. Should probably default to the USE_STRINGSTREAM
         * option.
         */
#ifdef HAS_SSTREAM
        /*
         * Read a line from the file
             *
             *  SOLARIS 10 NOTES: Optimized version of solaris seg faults
             *         without the '\n' argument for some reason, probably
             *         an internal solaris stl bug.
         */
        getline(s, rest,'\n');
        /*
         * In the transport database, we allow comment lines that
         * consist of '#' and '!' as the first character in the
         * in the line. We also don't bother to parse short lines that
         * can't possibly have enough data in them to comprise a
         * properly formatted record.
         */
        if (rest.size() > 5 && rest[0] != '#' && rest[0] != '!') {
            /*
             * copy the string into a stringstream and parse the line
             * into the trdata object
             */
            std::istringstream ioline(rest);
            trdata t;
            string nm;
            ioline >> nm >> t.geom >> t.welldepth >> t.diam
                   >> t.dipole >> t.polar >> t.rot;
            /*
             * Add the trdata object into the map database by making a
             * copy of it, and index it by the species name.
             */
            if (nm != "") {
                _trmap[nm] = t; // t.name] = t;
            }
        }
#else
        trdata t;
        string nm;
        s >> nm;
        if (nm[0] != '!' && !s.eof()) {
            s >> t.geom >> t.welldepth >> t.diam
              >> t.dipole >> t.polar >> t.rot;

            // get the rest of the line, in case there are comments
            getline(s, rest, '\n');
            if (nm != "") {
                _trmap[nm] = t; // t.name] = t;
            }
        }
#endif
    }
}


// add a NASA polynomial parameterization
static void addNASA(FILE* f,
                    const vector_fp& low, const vector_fp& high,
                    doublereal minx, doublereal midx,
                    doublereal maxx)
{

    fprintf(f,"    thermo = (\n");
    fprintf(f,"       NASA( [%8.2f, %8.2f], ", minx, midx);
    fprintf(f,"[%17.9E, %17.9E, \n", low[0], low[1]);
    fprintf(f,"              %17.9E, %17.9E, %17.9E,\n",
            low[2], low[3], low[4]);
    fprintf(f,"              %17.9E, %17.9E] ),\n", low[5], low[6]);
    fprintf(f,"       NASA( [%8.2f, %8.2f], ", midx, maxx);
    fprintf(f,"[%17.9E, %17.9E, \n", high[0], high[1]);
    fprintf(f,"              %17.9E, %17.9E, %17.9E,\n",
            high[2], high[3], high[4]);
    fprintf(f,"              %17.9E, %17.9E] )\n", high[5], high[6]);
    fprintf(f,"             )");
}



//! Add a NASA polynomial parameterization to the cti file
/*!
 * This little tidbit of code writes out the polynomials to the cti file
 */
static void addNASA9(FILE* f,
                     const std::vector<vector_fp*> &region_coeffs,
                     const vector_fp& minTemps, const vector_fp& maxTemps)
{
    size_t nReg = region_coeffs.size();
    if (minTemps.size() != nReg) {
        throw CanteraError("addNASA9", "incompat");
    }
    if (maxTemps.size() != nReg) {
        throw CanteraError("addNASA9", "incompat");
    }

    fprintf(f,"    thermo = (\n");
    for (size_t i = 0; i < nReg; i++) {
        double minT = minTemps[i];
        double maxT = maxTemps[i];
        const vector_fp& coeffs = *(region_coeffs[i]);
        if ((int) coeffs.size() != 9) {
            throw CanteraError("addNASA9", "incompat");
        }
        fprintf(f,"       NASA9( [%8.2f, %8.2f], ", minT, maxT);
        fprintf(f,"[%17.9E, %17.9E, %17.9E,\n", coeffs[0],
                coeffs[1], coeffs[2]);
        fprintf(f,"              %17.9E, %17.9E, %17.9E,\n",
                coeffs[3], coeffs[4], coeffs[5]);
        fprintf(f,"              %17.9E, %17.9E, %17.9E] )",
                coeffs[6], coeffs[7], coeffs[8]);
        if (i < nReg - 1) {
            fprintf(f,",\n");
        } else {
            fprintf(f,"\n");
        }
    }
    fprintf(f,"             )");
}


static void addTransportParams(FILE* f, string name)
{

    trdata td;
    if (_with_transport && _trmap.find(name) != _trmap.end()) {
        td = _trmap[name];
    } else {
        throw CanteraError("addTransportParams",
                           "no transport data for species "+name);
    }

    fprintf(f,",\n    transport = gas_transport(\n");
    int geom = td.geom;
    switch (geom) {
    case 0:
        fprintf(f,"                     geom = \"atom\",\n");
        break;
    case 1:
        fprintf(f,"                     geom = \"linear\",\n");
        break;
    case 2:
        fprintf(f,"                     geom = \"nonlinear\",\n");
        break;
    default:
        throw CanteraError("addTransportParams",
                           "Unrecognized geometry flag for species " + name);
    }
    fprintf(f,"                     diam = %g,\n",td.diam);
    fprintf(f,"                     well_depth = %g",td.welldepth);
    if (td.polar != 0.0) {
        fprintf(f,",\n                     polar = %g",td.polar);
    }
    if (td.dipole != 0.0) {
        fprintf(f,",\n                     dipole = %g",td.dipole);
    }
    if (td.rot != 0.0) {
        fprintf(f,",\n                     rot_relax = %g",td.rot);
    }
    fprintf(f,")");
}


static void addFalloff(FILE* f, string type,
                       const vector_fp& params)
{
    if (type == "Troe") {
        fprintf(f, "%s", (",\n         falloff = Troe(A = " +
                          fp2str(params[0]) + ", T3 = " +
                          fp2str(params[1]) + ", T1 = " +
                          fp2str(params[2])).c_str());
        if (params.size() >= 4) {
            fprintf(f, "%s", (", T2 = " + fp2str(params[3])).c_str());
        }
        fprintf(f, ")");
    } else if (type == "SRI") {
        fprintf(f, "%s", (",\n         falloff = SRI(A = " +
                          fp2str(params[0]) + ", B = " +
                          fp2str(params[1]) + ", C = " +
                          fp2str(params[2])).c_str());
        if (params.size() >= 5) {
            fprintf(f, "%s", (", D = " + fp2str(params[3]) +
                              ", E = " + fp2str(params[4])).c_str());
        }
        fprintf(f,  ")");
    }
}

/**
 * Write out a species cti block to the output file.
 *
 */
static void addSpecies(FILE* f, string idtag, const ckr::Species& sp)
{
    string spname = sp.name;
    if (spname.size() == 0) {
        throw CanteraError("addSpecies",
                           "Species name is empty");
    }
    fprintf(f,"\nspecies(name = \"%s\",\n",spname.c_str());
    int nel = static_cast<int>(sp.elements.size());
    int m, num;
    string nm, str="";
    for (m = 0; m < nel; m++) {
        /*
         * Copy the element name into the string, nm. Lower case the
         * second letter, if needed.
         */
        nm = sp.elements[m].name;
        nm[0] = (char) toupper(nm[0]);
        if (nm.size() == 2) {
            nm[1] = (char) tolower(nm[1]);
        }
        /*
         *  Obtain the current number of atoms in the species.
         *  Linearize the number (HKM question? can we employ real values here
         *  instead?)
         */
        num = int(sp.elements[m].number);
        /*
         * Add the name and number to end of the string, str
         */
        str += " "+nm+":"+int2str(num)+" ";

    }

    fprintf(f,"    atoms = \"%s\",\n", str.c_str());
    // Add the NASA block according to the thermoFormatType value
    if (sp.thermoFormatType == 0) {
        if (sp.lowCoeffs.size() == 0) {
            throw CanteraError("addSpecies",
                               "Low Nasa Thermo Polynomial was not found");
        }
        if (sp.highCoeffs.size() == 0) {
            throw CanteraError("addSpecies",
                               "High Nasa Thermo Polynomial was not found");
        }
        if (sp.tlow >= sp.thigh) {
            throw CanteraError("addSpecies",
                               "Low temp limit is greater or equal to high temp limit");
        }
        addNASA(f, sp.lowCoeffs, sp.highCoeffs,
                sp.tlow, sp.tmid, sp.thigh);
    } else if (sp.thermoFormatType == 1) {
        // This new typs is a multiregion 9 coefficient formulation
        addNASA9(f, sp.region_coeffs, sp.minTemps, sp.maxTemps);
    } else {
        throw CanteraError("addSpecies", "Unknown thermoFormatType");
    }

    if (_with_transport) {
        addTransportParams(f, sp.name);
    }
    if (sp.id != "" || sp.m_commentsRef != "") {
        fprintf(f,",\n    note = \"");
        if (sp.id != "") {
            fprintf(f, "%s", sp.id.c_str());
        }
        if (sp.m_commentsRef != "") {
            fprintf(f, " %s", sp.m_commentsRef.c_str());
        }
        fprintf(f, "\"");
    }
    fprintf(f,"\n       )\n");
}


static void addReaction(FILE* f, string idtag, int i,
                        const ckr::Reaction& rxn,
                        const ckr::ReactionUnits& runits, doublereal version)
{

    fprintf(f, "%s", ("\n#  Reaction " + int2str(i+1) + "\n").c_str());
    int nc = static_cast<int>(rxn.comment.size());
    vector<string> options;

    for (int nn = 0; nn < nc; nn++)
        if (rxn.comment[nn] != "") fprintf(f,  "# %s \n",
                                               rxn.comment[nn].c_str());

    string eqn = ckr::reactionEquation(rxn);

    if (rxn.isThreeBodyRxn) {
        fprintf(f,  "three_body_reaction( \"%s\",", eqn.c_str());
    } else if (rxn.isFalloffRxn) {
        fprintf(f,  "falloff_reaction( \"%s\",", eqn.c_str());
    } else {
        fprintf(f,  "reaction(  \"%s\",", eqn.c_str());
    }

    if (rxn.isFalloffRxn) {

        if (rxn.kf.type == ckr::Arrhenius) {
            fprintf(f,"\n         kf = [%10.5E, %g, %g]", rxn.kf.A, rxn.kf.n, rxn.kf.E);
        }
        if (rxn.kf_aux.type == ckr::Arrhenius) {
            fprintf(f,",\n         kf0   = [%10.5E, %g, %g]", rxn.kf_aux.A, rxn.kf_aux.n, rxn.kf_aux.E);
        }
        if (rxn.falloffType == ckr::Lindemann) {
            addFalloff(f, "Lindemann",rxn.falloffParameters);
        } else if (rxn.falloffType == ckr::Troe) {
            addFalloff(f, "Troe",rxn.falloffParameters);
        } else if (rxn.falloffType == ckr::SRI) {
            addFalloff(f, "SRI",rxn.falloffParameters);
        } else {
            throw CanteraError("addReaction","unknown falloff type");
        }
    } else {
        if (rxn.kf.type == ckr::Arrhenius) {
            fprintf(f,"  [%10.5E, %g, %g]", rxn.kf.A, rxn.kf.n, rxn.kf.E);
        } else {
            throw CanteraError("addReaction",
                               "unknown kf_type to reaction: " + int2str(rxn.kf.type));
        }
    }

    // reaction orders
    int nord = static_cast<int>(rxn.fwdOrder.size());
    if (nord > 0) {
        map<string, double>::const_iterator b = rxn.fwdOrder.begin(),
                                            e = rxn.fwdOrder.end();
        string estr = "";
        for (; b != e; ++b) {
            estr += " "+b->first+":"+fp2str(b->second)+" ";
        }
        fprintf(f,  ",\n         order = \"%s\"", estr.c_str());
    }

    int ne = static_cast<int>(rxn.e3b.size());
    if (rxn.thirdBody != "<none>") {
        if (rxn.thirdBody != "M") {
            ;
        } else if (ne > 0.0) {
            map<string, double>::const_iterator b = rxn.e3b.begin(),
                                                e = rxn.e3b.end();
            string estr = "";
            for (; b != e; ++b) {
                estr += " "+b->first+":"+fp2str(b->second)+" ";
            }
            fprintf(f,  ",\n         efficiencies = \"%s\"", estr.c_str());
        }
    }
    if (rxn.kf.A <= 0.0) {
        options.push_back("negative_A");
    }
    if (rxn.isDuplicate) {
        options.push_back("duplicate");
    }
    size_t nopt = options.size();
    if (nopt > 0) {
        fprintf(f,  ",\n         options = [");
        for (size_t n = 0; n < nopt; n++) {
            fprintf(f,  "\"%s\"", options[n].c_str());
            if (n < nopt-1) {
                fprintf(f,  ", ");
            }
        }
        fprintf(f,  "]");
    }
    fprintf(f,  ")\n");
}

void writeline(FILE* f)
{
    fprintf(f,  "#-------------------------------------------------------------------------------\n");
}

/*!
 * This routine is the main routine.
 *
 * @param r reference to a ckreader object that has already
 * read a chemkin formatted mechanism. This is the input to the routine.
 *
 * @param root Reference to the root node of an XML description of the
 *             mechanism. The node will be filled up with the description
 *             of the mechanism. This is the output to the routine.
 */
void ck2ct(FILE* f, string idtag, ckr::CKReader& r, bool hastransport)
{

    popError();
    doublereal version = 1.0;

    fprintf(f,  "units(length = \"cm\", time = \"s\", quantity = \"mol\", ");
    string e_unit = " ";
    int eunit = r.units.ActEnergy;
    if (eunit == ckr::Cal_per_Mole) {
        e_unit = "cal/mol";
    } else if (eunit == ckr::Kcal_per_Mole) {
        e_unit = "kcal/mol";
    } else if (eunit == ckr::Joules_per_Mole) {
        e_unit = "J/mol";
    } else if (eunit == ckr::Kjoules_per_Mole) {
        e_unit = "kJ/mol";
    } else if (eunit == ckr::Kelvin) {
        e_unit = "K";
    } else if (eunit == ckr::Electron_Volts) {
        e_unit = "eV";
    }
    fprintf(f,  "act_energy = \"%s\")\n\n", e_unit.c_str());


    fprintf(f,"\nideal_gas(name = \"%s\",\n",idtag.c_str());

    string enames;
    int nel = static_cast<int>(r.elements.size());
    int i;
    map<string, string> emap;
    string elnm;
    for (i = 0; i < nel; i++) {
        elnm = r.elements[i].name;
        elnm[0] = (char) toupper(elnm[0]);
        if (elnm.size() == 2) {
            elnm[1] = (char) tolower(elnm[1]);
        }
        emap[r.elements[i].name] = elnm;
        enames += " "+elnm+" ";
    }
    fprintf(f,"      elements = \"%s\",\n",enames.c_str());

    string spnames = "";
    int nsp = static_cast<int>(r.species.size());
    for (i = 0; i < nsp; i++) {
        spnames += " "+r.species[i].name+" ";
        if ((i+1) % 10 == 0) {
            spnames += "\n                  ";
        }
    }
    fprintf(f,"      species = \"\"\"%s\"\"\",\n", spnames.c_str());
    fprintf(f,"      reactions = \"all\",\n");
    if (hastransport) {
        fprintf(f,"      transport = \"Mix\",\n");
    }
    fprintf(f,"      initial_state = state(temperature = 300.0,\n");
    fprintf(f,"                        pressure = OneAtm)");
    fprintf(f,  "    )\n");

    fprintf(f,  "\n\n\n");
    writeline(f);
    fprintf(f,  "#  Species data \n");
    writeline(f);

    for (i = 0; i < nsp; i++) {
        addSpecies(f, idtag, r.species[i]);
    }

    fprintf(f,  "\n\n\n");
    writeline(f);
    fprintf(f,  "#  Reaction data \n");
    writeline(f);


    int nrxns = static_cast<int>(r.reactions.size());

    int irxn = 0;
    string idktag = idtag;
    for (i = 0; i < nrxns; i++) {

        // if krev.A is non-zero, then the reverse rate coefficient is
        // being explicitly specified rather than being computed from
        // thermochemistry. In this case, convert the reaction into
        // two irreversible reactions.

        if (r.reactions[i].krev.A != 0.0) {
            fprintf(f,  "\n# [CK Reaction (+%d)]\n",i+1);
            addReaction(f, idktag, irxn,
                        ckr::forwardReaction(r.reactions[i]), r.units, version);
            irxn++;
            fprintf(f,  "# [CK Reaction (-%d)]\n",i+1);
            addReaction(f, idktag, irxn,
                        ckr::reverseReaction(r.reactions[i]), r.units, version);
            irxn++;
        }

        // Otherwise, just add the whole reaction, which may or may
        // not be reversible.
        else {
            if (i != irxn) {
                fprintf(f,  "\n# [CK Reaction (%d)]\n",i+1);
            }
            addReaction(f, idktag, irxn, r.reactions[i],
                        r.units, version);
            irxn++;
        }
    }
}

/*
    static int fixtext(string infile, string outfile) {
        ifstream fin(infile.c_str());
        ofstream fout(outfile.c_str());
        if (!fout) {
            throw CanteraError("fixtext","could not open "+outfile+" for writing.");
        }
        char ch;
        char last_eol = ' ';
        const char char10 = char(10);
        const char char13 = char(13);
        string line;
        while (1 > 0) {
            line = "";
            while (1 > 0) {
                fin.get(ch);
                if (fin.eof()) break;
                if (ch == char13 || (ch == char10
                        && (last_eol != char13)))  {
                    last_eol = ch;
                    break;
                }
                if (isprint(ch)) line += ch;
            }
            fout << line << endl;
            if (fin.eof()) break;
        }
        fin.close();
        fout.close();
        return 0;
    }
*/

int convert_ck(const char* in_file, const char* db_file,
               const char* tr_file, const char* id_tag, bool debug, bool validate)
{
    ckr::CKReader r;

    r.validate = validate;
    r.debug = debug;
    //int i=1;

    string infile = string(in_file);
    string dbfile = string(db_file);
    string trfile = string(tr_file);
    //string outfile = string(out_file);
    string idtag = string(id_tag);
    string logfile;
    if (dbfile == "-") {
        dbfile = "";
    }
    if (trfile == "-") {
        trfile = "";
    }

    string::size_type idot = infile.rfind('.');
    string ctifile, ext;
    if (idot != string::npos) {
        ext = infile.substr(idot, infile.size());
        ctifile = infile.substr(0,idot)+".cti";
    } else {
        ctifile = infile+".cti";
    }

    FILE* f = fopen(ctifile.c_str(),"w");

    struct tm* newtime;
    time_t aclock;
    ::time(&aclock);                /* Get time in seconds */
    newtime = localtime(&aclock);   /* Convert time to struct tm form */

    try {

        //string tmpinfile = tmpDir()+
        //fixtext(infile, tmpinfile);

        //string tmpdbfile = "";
        //string tmptrfile = "";
        //if (dbfile != "") {
        //    tmpdbfile = tmpDir()+"/.tmp_"+dbfile;
        //    fixtext(dbfile, tmpdbfile);
        //}
        //if (trfile != "") {
        //    tmptrfile = tmpDir()+"/.tmp_"+trfile;
        //    fixtext(trfile, tmptrfile);
        //}

        logfile = "ck2cti.log";
        if (!r.read(infile, dbfile, logfile)) {
            throw CanteraError("convert_ck",
                               "error encountered in input file " + string(infile)
                               + "\nsee file ck2cti.log for more information.\n");
        }

        fprintf(f,  "#\n");
        fprintf(f,  "# Generated from file %s\n# by ck2cti on %s#\n",
                infile.c_str(), asctime(newtime));
        if (trfile != "") {
            fprintf(f,  "# Transport data from file %s.\n\n",
                    trfile.c_str());
            getTransportData(trfile);
        }
        bool hastransport = (trfile != "");
        ck2ct(f, idtag, r, hastransport);
    } catch (CanteraError& err) {
        err.save();
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}
}


