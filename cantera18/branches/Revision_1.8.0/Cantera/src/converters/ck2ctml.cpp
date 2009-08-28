/**
 * @file ck2ctml.cpp
 *
 * Convert CK-format reaction mechanism files to CTML format.
 *
 */
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <string>
#ifdef USE_STRINGSTREAM
#include <sstream>
#endif
using namespace std;

#include "CKReader.h"
#include "Reaction.h"
#include "writelog.h"


#include "xml.h"
#include "ck2ctml.h"

using namespace Cantera;

namespace ctml {

    struct trdata {
        trdata() {name = "-";}
        string name;
        int geom;
        doublereal welldepth, diam, dipole, polar, rot;
    };

    // add a NASA polynomial parameterization
    static void addNASA(XML_Node& node,  
        const vector_fp& low, const vector_fp& high,  
        doublereal minx, doublereal midx, 
        doublereal maxx) {
        XML_Node& f = node.addChild("NASA");
        if (minx != -999.0) f.addAttribute("Tmin",minx);
        if (midx != -999.0) f.addAttribute("Tmax",midx);
        addFloatArray(f,"coeffs",low.size(),&low[0]);
        XML_Node& fh = node.addChild("NASA");
        if (midx != -999.0) fh.addAttribute("Tmin",midx);
        if (maxx != -999.0) fh.addAttribute("Tmax",maxx);
        addFloatArray(fh,"coeffs",high.size(),&high[0]);
    }

    /*
     * HKM -
     * Took out because it was defined as static but not used
     * in this file
     * 
     *  static void addShomate(XML_Node& node, 
     *                         const vector_fp& low, const vector_fp& high,  
     *                         doublereal minx, doublereal midx, 
     *                         doublereal maxx) {
     *    XML_Node& f = node.addChild("Shomate");
     *    if (minx != -999.0) f.addAttribute("Tmin",minx);
     *    if (maxx != -999.0) f.addAttribute("Tmid",midx);
     *    if (maxx != -999.0) f.addAttribute("Tmax",maxx);
     *    addFloatArray(f,"low",low.size(),low.begin());
     *    addFloatArray(f,"high",high.size(),high.begin());
     * }
     */

    static void addArrhenius(XML_Node& node,  
        doublereal A, doublereal b, doublereal E, int order, 
        string id, string E_units) {

#ifdef OLD_VERSION
        // versions prior to 1.4.1
        string abe = fp2str(A)+" "+fp2str(b)+" "+fp2str(E);
        XML_Node& r = node.addChild("Arrhenius",abe);
        r.addAttribute("order",order);
#else
        // version 1.4.1
        XML_Node& rn = node.addChild("Arrhenius");
        if (id != "") rn.addAttribute("name",id);
        string units;
        if (order == 1) units = "/s";
        else if (order == 2) units = "cm3/mol/s";
        else if (order == 3) units = "cm6/mol2/s";
        else throw CanteraError("addArrhenius",
            "unsupported rxn order: "+int2str(order));
        addFloat(rn, "A", A,  units);
        addFloat(rn, "b", b);
        addFloat(rn, "E", E, E_units);
#endif
    }

    /*
     * HKM -
     * Took out because it was defined as static but not used
     * in this file
     * 
     * static void addRateCoeff(XML_Node& node, string title, string type,  
     *                          string direction, int order, doublereal A,
     *                          doublereal b, doublereal E, 
     *                          string unitsys, string E_units) {
     *    XML_Node& r = node.addChild("rateCoeff");
     *    r.addAttribute("id",title);
     *    if (type == "Arrhenius")
     *        addArrhenius(r, A, b, E, order, unitsys, E_units);
     *    r.addAttribute("direction",direction);
     *  }
    */

    static void addFalloff(XML_Node& node, string type, 
        const vector_fp& params) {
        int np = params.size();
        string p;
        for (int n = 0; n < np; n++) p += " "+fp2str(params[n])+" ";
        XML_Node& f = node.addChild("falloff",p);
        f.addAttribute("type",type);
    }

    /*
     * HKM -
     * Took out because it was defined as static but not used
     * in this file
     * 
     * static void addTroeFalloff(XML_Node& node, 
     *                            const vector_fp& params) {
     *   XML_Node& f = node.addChild("falloff");
     *   f.addAttribute("type","Troe");
     *   addFloat(f,"A",params[0]);
     *   addFloat(f,"T***",params[1]);
     *   addFloat(f,"T*",params[2]);
     *   if (params.size() == 4)
     *       addFloat(f,"T**",params[3]);
     *  }
     */
    /*
     * HKM -
     * Took out because it was defined as static but not used
     * in this file
     * 
     * static void addSRIFalloff(XML_Node& node, 
     *                         const vector_fp& params) {
     *   XML_Node& f = node.addChild("falloff");
     *   f.addAttribute("type","SRI");
     *   addFloat(f,"A",params[0]);
     *   addFloat(f,"B",params[1]);
     *   addFloat(f,"C",params[2]);
     *   if (params.size() == 5) {
     *       addFloat(f,"D",params[3]);
     *       addFloat(f,"E",params[4]);
     *   }
     *  }
     */
    /*
     * HKM -
     * Took out because it was defined as static but not used
     * in this file
     * 
     * static void addElement(XML_Node& node, string idtag, 
     *                        const ckr::Element& el) {
     *   XML_Node& elx = node.addChild("element");
     *   elx.addAttribute("id",idtag+"_e_"+el.name);
     *   string elnm = el.name;
     *   elnm[0] = toupper(elnm[0]);
     *   if (el.name.size() == 2) elnm[1] = tolower(elnm[1]);
     *   elx.addAttribute("name",elnm);
     *   elx.addAttribute("atomicWt",el.atomicWeight);
     *   //addFloat(elx, "mass", el.atomicWeight, "amu");
     *   //addString(elx, "elementSymbol", el.name);
     * }
     */
    
    /**
     * addSpecies():
     *
     *  Add a species element to an XML description of a mechanism
     *
     *  Input
     * --------
     *  node :  Reference to the current parent node, where the
     *          species description is going to be placed.
     */
    static void addSpecies(XML_Node& node, string idtag, const ckr::Species& sp) {
        string spname = sp.name;
        node.addComment(spname);
        XML_Node& spx = node.addChild("species");
        spx.addAttribute("id",idtag+"_s_"+sp.name);
        spx.addAttribute("name",sp.name);
        //spx.addAttribute("phaseType",sp.phase);
        if (stripws(sp.id) != "") spx.addChild("note",sp.id);
        //addString(spx,"comment","ideal gas");
        int nel = static_cast<int>(sp.elements.size());
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

	/*
	 * Add the child element, atomArray, to the species xml node. 
	 */
	spx.addChild("atomArray", str);
        addFloat(spx,"charge",charge);
        XML_Node& cp1 = spx.addChild("thermo");
        addNASA(cp1, sp.lowCoeffs, sp.highCoeffs, 
		sp.tlow, sp.tmid, sp.thigh); 
    }

    static void addReaction(XML_Node& node, string idtag, int i,
        const ckr::Reaction& rxn,
        const ckr::ReactionUnits& runits, doublereal version) {

        node.addComment(idtag+" reaction "+int2str(i+1));
        string eqn = ckr::reactionEquation(rxn);
        int eqlen = static_cast<int>(eqn.size());
        for (int nn = 0; nn < eqlen; nn++) {
            if (eqn[nn] == '<') eqn[nn] = '[';
            if (eqn[nn] == '>') eqn[nn] = ']';
        }
        //        node.addComment(eqn);
        XML_Node& r = node.addChild("reaction");
        r.addAttribute("id",idtag+"_rxn_"+int2str(i+1));
        r.addChild("equation", eqn);

        string e_unit;
        int eunit = runits.ActEnergy;
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

        int nr = static_cast<int>(rxn.reactants.size());
        int np = static_cast<int>(rxn.products.size());

        int n;
        doublereal order = 0.0;

        string rct="", prd="";
        for (n = 0; n < nr; n++) {
            //XML_Node& reac = r.addChild("reactant");
            //reac.addAttribute("name",rxn.reactants[n].name);
            rct += (" " + rxn.reactants[n].name + ":"
		    + int2str((int) rxn.reactants[n].number) + " ");
            //reac.addAttribute("phase",idtag);
            doublereal nn = rxn.reactants[n].number;
            order += nn;
        }
        r.addChild("reactants",rct);
            //if (nn != 1) addInteger(reac, "number", int(nn));
        
        //                reac.addAttribute("phase",idtag);
        for (n = 0; n < np; n++) {
            prd += (" " + rxn.products[n].name + ":"
		    + int2str((int) rxn.products[n].number) + " ");
            //doublereal nn = rxn.products[n].number;
            //prod.addAttribute("name",rxn.products[n].name);
            //if (nn != 1) addInteger(prod, "number", int(nn));
        }         
        r.addChild("products",prd);
            
        if (rxn.isChemActRxn || rxn.isThreeBodyRxn) order += 1.0;
        

        XML_Node& kf = r.addChild("rateCoeff");
        //kf.addAttribute("units","mol,cm,s");
        //kf.addAttribute("Eunits",e_unit);

        //kf.addAttribute("id",r["id"]+"_kf");
        if (rxn.kf.type == ckr::Arrhenius)
            addArrhenius(kf, rxn.kf.A, rxn.kf.n, rxn.kf.E, 
                int(order), "", e_unit);

        if (rxn.isFalloffRxn) {
            addArrhenius(kf, rxn.kf_aux.A, rxn.kf_aux.n, rxn.kf_aux.E, 
                int(order+1), "k0", e_unit);
            
            if (rxn.falloffType == ckr::Lindemann)
                addFalloff(kf,"Lindemann",rxn.falloffParameters);
            else if (rxn.falloffType == ckr::Troe)
                addFalloff(kf,"Troe",rxn.falloffParameters);
            else if (rxn.falloffType == ckr::SRI)
                addFalloff(kf,"SRI",rxn.falloffParameters);
            else 
                throw CanteraError("addReaction","unknown falloff type");
        }

        int ne = static_cast<int>(rxn.e3b.size());
        if (rxn.thirdBody != "<none>") {
            if (rxn.thirdBody != "M") {
                XML_Node& e3b = kf.addChild("efficiencies",rxn.thirdBody+":1.0");
                e3b.addAttribute("default",0.0);
                //addFloat(e3b, rxn.thirdBody, 1.0);
            }
            else if (ne > 0.0) {
                map<string, double>::const_iterator b = rxn.e3b.begin(), 
                                                    e = rxn.e3b.end();
                string estr = "";
                for (; b != e; ++b) {
                    estr += " "+b->first+":"+fp2str(b->second)+" ";
                    //addFloat(e3b,b->first,b->second);
                }
                estr += "\n";
                XML_Node& e3b = kf.addChild("efficiencies",estr);
                e3b.addAttribute("default",1.0);
            }
        }

        if (rxn.isReversible)
            r.addAttribute("reversible","yes");
        else
            r.addAttribute("reversible","no");
        //addBool(r,"reversible",rxn.isReversible);        

        if (rxn.duplicate != 0)
            r.addChild("duplicate","idtag_rxn_"+int2str(rxn.duplicate));
        //            addString(r, "duplicate","reaction_"+int2str(rxn.duplicate));
    
        if (rxn.isFalloffRxn) 
            r.addAttribute("type","falloff");
        if (rxn.isChemActRxn) 
            r.addAttribute("type","chemAct");
        if (rxn.isThreeBodyRxn) 
            r.addAttribute("type","threeBody");
    }

    static void addState(XML_Node& phase, string x) {
        XML_Node& state = phase.addChild("state");
        XML_Node& temp = state.addChild("temperature",300.0);
        temp.addAttribute("units","K");
        XML_Node& pres = state.addChild("pressure",1.0);
        pres.addAttribute("units","atm");
        state.addChild("moleFractions",x+":1.0");
    }

    /*!
     *  addTransport() searches a transport database, pointed to by
     *  the istream, s,
     *  for the species names listed in node. It then adds the
     *  pertinent transport properties to the XML_Node database
     *  pointed to by node.
     */
    static void addTransport(istream& s, XML_Node& node) {
        /*
         *  The first thing we will do is to read the entire transport
         *  database and place its contents into a map structure,
         *  indexed by the name of the species.
         */
        map<string, trdata> indx;
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
	    ioline >> t.name >> t.geom >> t.welldepth >> t.diam 
		   >> t.dipole >> t.polar >> t.rot;
	    /*
	     * Add the trdata object into the map database by making a
	     * copy of it, and index it by the species name.
	     */
	    if (t.name != "") { 
	      indx[t.name] = t;
	    }
	  }
#else
	  trdata t;
	  s >> t.name;
          if (t.name[0] != '!' && !s.eof()) {
              s >> t.geom >> t.welldepth >> t.diam 
                >> t.dipole >> t.polar >> t.rot;

              // get the rest of the line, in case there are comments
              getline(s, rest);
              if (t.name != "") { 
                  indx[t.name] = t;
              }
          }
#endif
	}

        vector<XML_Node*> sp;
        node.getChildren("species",sp);
        int ns = static_cast<int>(sp.size());
        for (int n = 0; n < ns; n++) {
            XML_Node& spx = *sp[n];
            string nm = spx["name"];
            trdata t = indx[nm];
            if (t.name == "-") {
                throw CanteraError("addTransport",
                    "no transport data for species "+nm);
            }
            XML_Node& tr = spx.addChild("transport");
        
            switch (t.geom) {
            case 0:
                addString(tr,"geometry","atom"); break;
            case 1:
                addString(tr,"geometry","linear"); break;
            case 2:
                addString(tr,"geometry","nonlinear"); break;
            default: ;
            }
            //if (t.welldepth != 0.0)
                addFloat(tr,"LJ_welldepth",t.welldepth,"Kelvin");
                //if (t.diam != 0.0)
                addFloat(tr,"LJ_diameter",t.diam,"A");
                //if (t.dipole != 0.0)
                addFloat(tr,"dipoleMoment",t.dipole,"Debye");
                //if (t.polar != 0.0)
                addFloat(tr,"polarizability",t.polar,"A^3");
                //if (t.rot != 0.0)
                addFloat(tr,"rotRelax",t.rot);
        }
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
    void ck2ctml(string idtag, ckr::CKReader& r, XML_Node& root) {

        popError();
        doublereal version = 1.3;

        XML_Node& ph = root.addChild("phase");
        ph.addAttribute("id",idtag);

        addState(ph, r.species[0].name);

        XML_Node& eos = ph.addChild("thermo");
        eos.addAttribute("model","IdealGas");

        string enames;
        int nel = static_cast<int>(r.elements.size());
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
        XML_Node& earray = ph.addChild("elementArray",enames);
        earray.addAttribute("datasrc","elements.xml");

        string spnames = "";
        root.addComment("species data");
        XML_Node& spdata = root.addChild("speciesData");
        spdata.addAttribute("id",idtag+"_species_data");
        int nsp = static_cast<int>(r.species.size());
        for (i = 0; i < nsp; i++) {
            spnames += " "+r.species[i].name+" ";
            if ((i+1) % 10 == 0) spnames += "\n";
            addSpecies(spdata, idtag, r.species[i]);
        }
        XML_Node& sparray = ph.addChild("speciesArray",spnames+"\n");
        sparray.addAttribute("datasrc","#"+idtag+"_species_data");
        //sparray.addAttribute("database",idtag+"_species_data");

        XML_Node& rxns = ph.addChild("reactionArray");
        rxns.addAttribute("datasrc","#"+idtag+"_rxn_data");
        //rxns.addAttribute("datasrc",idtag+"_rxn_data");
        XML_Node& incl = rxns.addChild("include");
        XML_Node& ktype = ph.addChild("kinetics");
        ktype.addAttribute("model","GasKinetics");

        root.addComment("reaction data");
        XML_Node& kin = root.addChild("reactionData");
        kin.addAttribute("id",idtag+"_rxn_data");

        int nrxns = static_cast<int>(r.reactions.size());
        incl.addAttribute("prefix",idtag+"_rxn_");

        int irxn = 0;
        string idktag = idtag;
        for (i = 0; i < nrxns; i++) {

            // if krev.A is non-zero, then the reverse rate coefficient is
            // being explicitly specified rather than being computed from 
            // thermochemistry. In this case, convert the reaction into 
            // two irreversible reactions.

            if (r.reactions[i].krev.A != 0.0) {
                addReaction(kin, idktag, irxn, 
                    ckr::forwardReaction(r.reactions[i]), r.units, version);
                irxn++;
                addReaction(kin, idktag, irxn, 
                    ckr::reverseReaction(r.reactions[i]), r.units, version);
                irxn++;
            }

            // Otherwise, just add the whole reaction, which may or may
            // not be reversible.
            else {
                addReaction(kin, idktag, irxn, r.reactions[i], 
                    r.units, version);
                irxn++;
            }
        }
        incl.addAttribute("min",1);
        incl.addAttribute("max", irxn);
    }


    int convert_ck(const char * const in_file, const char * const db_file,
        const char * const tr_file, const char * const out_file, const char* const id_tag) {
        ckr::CKReader r;
        r.validate = true;
        //int i=1;
 
        string infile = string(in_file);
        string dbfile = string(db_file);
        string trfile = string(tr_file);
        string outfile = string(out_file);
        string idtag = string(id_tag);
        string logfile;
        if (dbfile == "-") dbfile = "";
        if (trfile == "-") trfile = "";
 
        try {

            logfile = "ck2ctml.log";
            if (!r.read(infile, dbfile, logfile)) {
                throw CanteraError("convert_ck",
                    "error encountered in input file " + string(infile) 
                    + "\nsee file ck2ctml.log for more information.\n");
            }

            XML_Node root("ctml");
            root["version"] = CTML_Version;
            root.addComment("generated from "+infile+" by ck2ctml.");
            if (trfile != "") 
            root.addComment("transport data from "+trfile+".");

            ck2ctml(idtag, r, root);
    
            if (trfile != "") {
                ifstream ftr(trfile.c_str());
                if (!ftr) {
                    throw CanteraError("convert_ck",
                        "could not open transport database.");
                }
                XML_Node* sparray = root.findByName("speciesData");
                if (sparray) {
                    addTransport(ftr, *sparray);
                }
                else
                    throw CanteraError("convert_ck", "could not find speciesData");
                ftr.close();
            }

            if (outfile == "") {
                root.writeHeader(cout);
                root.write(cout);
            }
            else {
                ofstream ff(outfile.c_str());
                root.writeHeader(ff);
                root.write(ff);
                ff.close(); 
                cout << "CTML file " << outfile << " written." << endl;
            }  
        }
        catch (CanteraError) {
            return -1;
        }
        
        return 0;
    }

}


