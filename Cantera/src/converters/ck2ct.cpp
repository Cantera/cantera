/**
 * @file ck2ct.cpp
 *
 * Convert CK-format reaction mechanism files to Cantera input format.
 *
 */
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <string>

// on cygwin, #include <sstream> doesn't work
#if defined(__CYGWIN__)
#undef USE_STRINGSTREAM
#else
#define USE_STRINGSTREAM
#endif

#ifdef USE_STRINGSTREAM
#include <sstream>
#endif
using namespace std;

#include "CKReader.h"
#include "Reaction.h"
#include "writelog.h"

#include "ck2ct.h"
#include <time.h>
#include "../ct_defs.h"
#include "ctml.h"
using namespace Cantera;

namespace pip {

    struct trdata {
        int geom;
        doublereal welldepth, diam, dipole, polar, rot;
    };

    static map<string, trdata> _trmap;
    static bool _with_transport = false;

    static void getTransportData(string trfile) {

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
#ifdef USE_STRINGSTREAM
	  /*
	   * Read a line from the file
	   */
	  getline(s, rest);
	  /*
	   * In the transport database, we allow comment lines that
	   * consist of '#' and '!' as the first character in the
	   * in the line. We also don't bother to parse short lines that
	   * can't possibly have enough data in them to comprise a
	   * properly formatted record.
	   */
	  if (rest[0] != '#' && rest[0] != '!' && rest.size() > 5) {
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
	    getline(s, rest);
	    if (nm != "") { 
	      _trmap[nm] = t; // t.name] = t;
	    }
	  }
#endif
	}
    }


    // add a NASA polynomial parameterization
    static void addNASA(  
        const vector_fp& low, const vector_fp& high,  
        doublereal minx, doublereal midx, 
        doublereal maxx) {

        printf("    thermo = (\n");
        printf("       NASA( [%8.2f, %8.2f], ", minx, midx);
        printf("[%17.9E, %17.9E, \n", low[0], low[1]);
        printf("              %17.9E, %17.9E, %17.9E,\n", low[2], low[3], low[4]);
        printf("              %17.9E, %17.9E] ),\n", low[5], low[6]);
        printf("       NASA( [%8.2f, %8.2f], ", midx, maxx);
        printf("[%17.9E, %17.9E, \n", high[0], high[1]);
        printf("              %17.9E, %17.9E, %17.9E,\n", high[2], high[3], high[4]);
        printf("              %17.9E, %17.9E] )\n", high[5], high[6]);
        printf("             )");
    }


    static void addTransportParams(string name) {

        trdata td;
        if (_with_transport && _trmap.find(name) != _trmap.end()) {
            td = _trmap[name];
        }
        else {
            throw CanteraError("addTransportParams",
                "no transport data for species "+name);
        }

        printf(",\n    transport = gas_transport(\n");
        int geom = td.geom;
        switch (geom) {
        case 0: printf("                     geom = \"atom\",\n"); break;
        case 1: printf("                     geom = \"linear\",\n"); break;
        case 2: printf("                     geom = \"nonlinear\",\n"); break;
        }
        printf("                     diam = %8.2f,\n",td.diam);
        printf("                     well_depth = %8.2f",td.welldepth);
        if (td.polar != 0.0)
            printf(",\n                     polar = %8.2f",td.polar);
        if (td.dipole != 0.0)
            printf(",\n                     dipole = %8.2f",td.dipole);
        if (td.rot != 0.0)
            printf(",\n                     rot_relax = %8.2f",td.rot);
        printf(")");
    }


    static void addFalloff(string type, 
        const vector_fp& params) {
        if (type == "Troe") {
            cout << ",\n         falloff = Troe(A = " 
                 << fp2str(params[0]) << ", T3 = "
                 << fp2str(params[1]) << ", T1 = "
                 << fp2str(params[2]);
            if (params.size() >= 4) {
                cout << ", T2 = " << fp2str(params[3]);
            }
            cout << ")";
        }
        else if (type == "SRI") {
            cout << ",\n         falloff = SRI(A = " 
                 << fp2str(params[0]) << ", B = "
                 << fp2str(params[1]) << ", C = "
                 << fp2str(params[2]);
            if (params.size() >= 5) {
                cout << ", D = " << fp2str(params[3])
                     << ", E = " << fp2str(params[4]);
            }
            cout << ")";
        }
    }

    /**
     * addSpecies():
     *
     */
    static void addSpecies(string idtag, const ckr::Species& sp) {
        string spname = sp.name;
        printf("\nspecies(name = \"%s\",\n",spname.c_str());
        int nel = sp.elements.size();
        int m, num;
        string nm, str="";
        doublereal charge = 0.0;
        for (m = 0; m < nel; m++) {
	  /*
	   * Copy the element name into the string, nm. Lower case the
	   * second letter, if needed.
	   */
	  nm = sp.elements[m].name;
          nm[0] = toupper(nm[0]);
	  if (nm.size() == 2) nm[1] = tolower(nm[1]);
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

          /* if the species contains the special element E (electron),
           * then set the charge.
           */
          if (nm == "E") charge = -sp.elements[m].number;
        }

        printf("    atoms = \"%s\",\n", str.c_str());
        addNASA(sp.lowCoeffs, sp.highCoeffs, 
            sp.tlow, sp.tmid, sp.thigh);

        if (_with_transport)
            addTransportParams(sp.name);
        if (sp.id != "") printf(",\n    note = \"%s\"", sp.id.c_str());
        printf("\n       )\n");
    }


    static void addReaction(string idtag, int i,
        const ckr::Reaction& rxn,
        const ckr::ReactionUnits& runits, doublereal version) {

        cout << "\n#  Reaction " << i+1 << endl;
        int nc = rxn.comment.size();
        for (int nn = 0; nn < nc; nn++) 
            if (rxn.comment[nn] != "") cout << "# " << rxn.comment[nn] << endl;
        string eqn = ckr::reactionEquation(rxn);

        if (rxn.isThreeBodyRxn) 
            cout << "three_body_reaction( \"" << eqn << "\", ";
        else if (rxn.isFalloffRxn)
            cout << "falloff_reaction( \"" << eqn << "\", ";
        else
            cout << "reaction( \"" << eqn << "\", ";

        if (rxn.isFalloffRxn) {

            if (rxn.kf.type == ckr::Arrhenius) {
                printf("\n         kf = [%10.5E, %g, %g]", rxn.kf.A, rxn.kf.n, rxn.kf.E);
            }
            if (rxn.kf_aux.type == ckr::Arrhenius) {
                printf(",\n         kf0   = [%10.5E, %g, %g]", rxn.kf_aux.A, rxn.kf_aux.n, rxn.kf_aux.E);
            }
            if (rxn.falloffType == ckr::Lindemann)
                addFalloff("Lindemann",rxn.falloffParameters);
            else if (rxn.falloffType == ckr::Troe)
                addFalloff("Troe",rxn.falloffParameters);
            else if (rxn.falloffType == ckr::SRI)
                addFalloff("SRI",rxn.falloffParameters);
            else 
                throw CanteraError("addReaction","unknown falloff type");
        }
        else {
            if (rxn.kf.type == ckr::Arrhenius) {
                printf("  [%10.5E, %g, %g]", rxn.kf.A, rxn.kf.n, rxn.kf.E);
            }
        }


        int ne = rxn.e3b.size();
        if (rxn.thirdBody != "<none>") {
            if (rxn.thirdBody != "M") {
                ;
            }
            else if (ne > 0.0) {
                map<string, double>::const_iterator b = rxn.e3b.begin(), 
                                                    e = rxn.e3b.end();
                string estr = "";
                for (; b != e; ++b) {
                    estr += " "+b->first+":"+fp2str(b->second)+" ";
                }
                cout << ",\n         efficiencies = \"" << estr << "\"";
            }
        }
        if (rxn.isDuplicate) {
            cout << ",\n         options = \'duplicate\'";
        }
    cout << ")" << endl;
    }

    void writeline() {
        cout << "#-------------------------------------------------------------------------------" << endl;
    }

    /*!
     * This routine is the main routine. It
     * 
     * @param r reference to a ckreader object that has already read a chemkin formatted
     *          mechanism. This is the input to the routine.
     * @param root Reference to the root node of an XML description of the
     *             mechanism. The node will be filled up with the description
     *             of the mechanism. This is the output to the routine.
     */
    void ck2ct(string idtag, ckr::CKReader& r, bool hastransport) {

        popError();
        doublereal version = 1.0;

        cout << "units(length = \"cm\", time = \"s\", quantity = \"mol\", ";
        string e_unit;
        int eunit = r.units.ActEnergy;
        if (eunit == ckr::Cal_per_Mole)
            e_unit = "cal/mol";
        else if (eunit == ckr::Kcal_per_Mole)
            e_unit = "kcal/mol";
        else if (eunit == ckr::Joules_per_Mole)
            e_unit = "J/mol";
        else if (eunit == ckr::Kjoules_per_Mole) 
            e_unit = "kJ/mol";
        else if (eunit == ckr::Kelvin) 
            e_unit = "K";
        else if (eunit == ckr::Electron_Volts)
            e_unit = "eV";
        cout << "act_energy = " << "\"" << e_unit << "\")\n\n";


        printf("\nideal_gas(name = \"%s\",\n",idtag.c_str());

        string enames;
        int nel = r.elements.size();
        int i;
        map<string, string> emap;
        string elnm;
        for (i = 0; i < nel; i++) {
            elnm = r.elements[i].name;
            elnm[0] = toupper(elnm[0]);
            if (elnm.size() == 2) elnm[1] = tolower(elnm[1]);
            emap[r.elements[i].name] = elnm;
            enames += " "+elnm+" ";
            //addElement(earray, idtag, r.elements[i]);
        }
        printf("      elements = \"%s\",\n",enames.c_str());

        string spnames = "";
        int nsp = r.species.size();
        for (i = 0; i < nsp; i++) {
            spnames += " "+r.species[i].name+" ";
            if ((i+1) % 10 == 0) spnames += "\n                  ";
        }
        printf("      species = \"\"\"%s\"\"\",\n", spnames.c_str());
        printf("      reactions = \"all\",\n");
        if (hastransport) {
            printf("      transport = \"Mix\",\n");
        } 
        printf("      initial_state = state(temperature = 300.0,\n");
        printf("                        pressure = OneAtm)");
        cout << "    )" << endl;

        cout << "\n\n\n";
        writeline();
        cout << "#  Species data \n";
        writeline();

        for (i = 0; i < nsp; i++) {
            addSpecies(idtag, r.species[i]);
        }

        cout << "\n\n\n";
        writeline();
        cout << "#  Reaction data \n";
        writeline();


        int nrxns = r.reactions.size();

         int irxn = 0;
         string idktag = idtag;
        for (i = 0; i < nrxns; i++) {

            // if krev.A is non-zero, then the reverse rate coefficient is
            // being explicitly specified rather than being computed from 
            // thermochemistry. In this case, convert the reaction into 
            // two irreversible reactions.

            if (r.reactions[i].krev.A != 0.0) {
                addReaction(idktag, irxn, 
                    ckr::forwardReaction(r.reactions[i]), r.units, version);
                irxn++;
                addReaction(idktag, irxn, 
                    ckr::reverseReaction(r.reactions[i]), r.units, version);
                irxn++;
            }

            // Otherwise, just add the whole reaction, which may or may
            // not be reversible.
            else { 
                addReaction(idktag, irxn, r.reactions[i], 
                    r.units, version);
                irxn++;
            }
        }
//         incl.addAttribute("min",1);
//         incl.addAttribute("max", irxn);
    }



    int convert_ck(const char* in_file, const char* db_file,
        const char* tr_file, const char* id_tag) {
        ckr::CKReader r;

        r.validate = true;
        //int i=1;

        string infile = string(in_file);
        string dbfile = string(db_file);
        string trfile = string(tr_file);
        //string outfile = string(out_file);
        string idtag = string(id_tag);
        string logfile;
        if (dbfile == "-") dbfile = "";
        if (trfile == "-") trfile = "";

        struct tm *newtime;
        time_t aclock;
        ::time( &aclock );              /* Get time in seconds */
        newtime = localtime( &aclock ); /* Convert time to struct tm form */
 
        try {

            logfile = "ck2cti.log";
            if (!r.read(infile, dbfile, logfile)) {
                throw CanteraError("convert_ck",
                    "error encountered in input file " + string(infile) 
                    + "\nsee file ck2cti.log for more information.\n");
            }

            cout << "#" << endl;
            cout << "# Generated from file " 
                 << infile << "\n# by ck2cti on " << asctime(newtime) << "#" << endl;
            if (trfile != "") {
                cout << "# Transport data from file "+trfile+".\n" << endl;
                getTransportData(trfile);
            }
            bool hastransport = (trfile != "");
            ck2ct(idtag, r, hastransport);
        }
        catch (CanteraError) {
            return -1;
        }
        return 0;
    }
}


