/*
 *  @file importCTML.cpp
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2002  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "importCTML.h"
#include "mix_defs.h"
#include <time.h>

//   STL includes
#include <map>
#include <string>
#include <vector>
using namespace std;

//   Cantera includes
#include "speciesThermoTypes.h"
#include "ThermoPhase.h"
#include "ThermoFactory.h"
#include "SpeciesThermoFactory.h"
#include "KineticsFactory.h"
#include "reaction_defs.h"
#include "ReactionData.h"
#include "global.h"
#include "stringUtils.h"

#include "xml.h"
#include "ctml.h"
using namespace ctml;

#include <stdio.h>

namespace Cantera {

    typedef vector<XML_Node*> nodeset_t;
    typedef XML_Node node_t;

    /// Number of reactant molecules
    static int nReacMolecules(ReactionData& r) {
        return accumulate(r.rstoich.begin(), r.rstoich.end(), 0);
    }

    const doublereal DefaultPref = 1.01325e5;   // one atm

    /**
     * Get an XML tree from a file. If successful, a pointer to the
     * document root is returned. If not, a null pointer is returned.
     */
    XML_Node* get_XML(string file) {
        string inname = findInputFile(file);
        if (inname == "") return 0;
        ifstream fin(inname.c_str());

        XML_Node* rootPtr = new XML_Node;
        try {
            rootPtr->build(fin);
        }
        catch (...) {
            return 0;
        }
        fin.close();
        return rootPtr;
    }        


    /** 
     * Install a NASA polynomial thermodynamic property
     * parameterization for species k.
     */
    void installNasaThermo(SpeciesThermo& sp, int k, XML_Node& f) {
        doublereal tmin, tmid, tmax;
        tmin = fpValue(f["Tmin"]);
        tmid = fpValue(f["Tmid"]);
        tmax = fpValue(f["Tmax"]);

        vector<XML_Node*> fa;
        f.getChildren("floatArray",fa);
        vector_fp c0, c1;
        getFloatArray(*fa[0], c0, false);
        getFloatArray(*fa[1], c1, false);
        array_fp c(15);
        c[0] = tmid;
        doublereal p0 = OneAtm;
        if ((*fa[0])["title"] == "low") {
            c[1] = c0[5];
            c[2] = c0[6];
            copy(c0.begin(), c0.begin()+5, c.begin() + 3);
            c[8] = c1[5];
            c[9] = c1[6];
            copy(c1.begin(), c1.begin()+5, c.begin() + 10);
        }
        else {
            c[1] = c1[5];
            c[2] = c1[6];
            copy(c1.begin(), c1.begin()+5, c.begin() + 3);
            c[8] = c0[5];
            c[9] = c0[6];
            copy(c0.begin(), c0.begin()+5, c.begin() + 10);
        }
        sp.install(k, NASA, c.begin(), tmin, tmax, p0);
    }

    /** 
     * Install a Shomate polynomial thermodynamic property
     * parameterization for species k.
     */
    void installShomateThermo(SpeciesThermo& sp, int k, XML_Node& f) {
        doublereal tmin, tmid, tmax;
        tmin = fpValue(f["Tmin"]);
        tmid = fpValue(f["Tmid"]);
        tmax = fpValue(f["Tmax"]);

        vector<XML_Node*> fa;
        f.getChildren("floatArray",fa);
        vector_fp c0, c1;
        getFloatArray(*fa[0], c0, false);
        getFloatArray(*fa[1], c1, false);
        array_fp c(15);
        c[0] = tmid;
        doublereal p0 = OneAtm;
        if ((*fa[0])["title"] == "low") {
            copy(c0.begin(), c0.end(), c.begin() + 1);
            copy(c1.begin(), c1.end(), c.begin() + 8);
        }
        else {
            copy(c1.begin(), c1.end(), c.begin() + 1);
            copy(c0.begin(), c0.end(), c.begin() + 8);
        }
        sp.install(k, SHOMATE, c.begin(), tmin, tmax, p0);
    }

    /** 
     * Install a Shomate polynomial thermodynamic property
     * parameterization for species k.
     */
    void installSimpleThermo(SpeciesThermo& sp, int k, XML_Node& f) {
        doublereal tmin, tmax;
        tmin = fpValue(f["Tmin"]);
        tmax = fpValue(f["Tmax"]);
        if (tmax == 0.0) tmax = 1.0e30;

        //map<string, doublereal> fp;
        //getFloats(f, fp);
        vector_fp c(4);
        c[0] = getFloat(f, "t0", "-");
        c[1] = getFloat(f, "h0", "-");
        c[2] = getFloat(f, "s0", "-");
        c[3] = getFloat(f, "cp0", "-");
        doublereal p0 = OneAtm;
        sp.install(k, SIMPLE, c.begin(), tmin, tmax, p0);
    }

    bool installSpecies(int k, XML_Node& s, thermo_t& p, SpeciesThermo& spthermo, int rule) {

            // get the composition of the species
            XML_Node& a = s.child("atomArray");
            map<string,string> comp;
            getMap(a, comp);

            // check that all elements in the species
            // exist in 'p'
            map<string,string>::const_iterator _b = comp.begin();
            for (; _b != comp.end(); ++_b) {
                if (p.elementIndex(_b->first) < 0) {
                    if (rule == 0) 
                        throw CanteraError("installSpecies", 
                            "species " + s["name"] + 
                            " contains undeclared element " + _b->first);
                    else
                        return false;
                }
            }

            int m, nel = p.nElements();
            vector_fp ecomp(nel, 0.0);            
            for (m = 0; m < nel; m++) {
                ecomp[m] = atoi(comp[p.elementName(m)].c_str());
            }

	    /*
	     * Define a map and get all of the floats in the
	     * current XML species block
	     */
            map<string, double> fd;
            getFloats(s, fd);

	    /*
	     * Set a default for the size parameter to one
	     */
            if (fd["size"] == 0.0) fd["size"] = 1.0;
            p.addUniqueSpecies(s["name"], ecomp.begin(),
				 fd["charge"], fd["size"]);

            // get thermo
            XML_Node& thermo = s.child("thermo");
            vector<XML_Node*> tp = thermo.children();
            int nc = tp.size();
            if (nc == 1) {
                XML_Node& f = *tp[0];
                if (f.name() == "NASA") {
                    installNasaThermo(spthermo, k, f);
                }
                else if (f.name() == "Shomate") {
                    installShomateThermo(spthermo, k, f);
                }
                else if (f.name() == "const_cp") {
                    installSimpleThermo(spthermo, k, f);
                }
                else 
                    throw CanteraError("importCTML",
                        "Unsupported species thermo parameterization"
                        " for species "+s["name"]);
            }
            else 
                throw CanteraError("importCTML",
                    "Multiple thermo parameterizations given for "
                    "species "+s["name"]);

            return true;
        }


    /**
     * Get the reactants or products of a reaction.
     */
    bool getReagents(XML_Node& rxn, kinetics_t& kin, int rp,
        string default_phase, 
        vector_int& spnum, vector_int& stoich, vector_fp& order,
        int rule) {

        string rptype;
        if (rp == 1) rptype = "reactants";
        else rptype = "products";
        XML_Node& rg = rxn.child(rptype);
        vector<string> key, val;
        getPairs(rg, key, val);

        int ns = key.size();
        int stch, isp;
        doublereal ord;
        string ph, sp;
        for (int n = 0; n < ns; n++) {
            sp = key[n];
            ph = ""; //snode["phase"];
            //if (ph == "") ph = default_phase;
            isp = kin.kineticsSpeciesIndex(sp,"<any>");
            if (isp < 0) {
                if (rule == 1) 
                    return false;
                else {
                    throw CanteraError("getReagents",
                        "undeclared reactant or product species "+sp);
                    return false;
                }
            }
            spnum.push_back(isp);
            stch = atoi(val[n].c_str());
            stoich.push_back(stch);
            ord = doublereal(stch);
            order.push_back(ord);
        }
        return true;
    }


    void getArrhenius(XML_Node& node, int& order, doublereal& A, doublereal& b, 
        doublereal& E) {
        
        // get rxn order to do unit conversion for pre-exponential
        order = intValue(node["order"]);
        
        //nodeset_t c = node.children();

        vector<string> abe;
        getStringArray(node, abe);
        A = fpValue(abe[0]);
        b = fpValue(abe[1]);
        E = fpValue(abe[2]);

        string u = (*node.parent())["units"];
        string eu = (*node.parent())["Eunits"];

        doublereal cmult = 1.0;
        if (u != "") {
            if (u == "mol,cm,s") 
                cmult = 1.0e-6 / CtMoles_per_mole; 
            else if (u == "molec,cm,s") 
                cmult = 1.0e-6*Avogadro;
            else 
                throw CanteraError("getArrhenius","unknown units for A");
        }
        A *= pow(cmult, order - 1);

        doublereal gasConstant = 1.0;
        if (eu != "") {
            if (eu == "cal/mol") 
                gasConstant = 1.987;
            else if (eu == "kcal/mol") 
                gasConstant = 1.987e-3;
            else if (eu == "J/mol") 
                gasConstant = 8.314;
            else if (eu == "kJ/mol") 
                gasConstant = 8.314e-3;
            else if (eu == "K") 
                gasConstant = 1.0;
            else if (eu == "eV") 
                gasConstant = 1.0/11600.0;
            else 
                throw CanteraError("getArrhenius",
                    "unknown units for activation energy: "+eu);
        }
        E /= gasConstant;
    }                


    /**
     * Get falloff parameters for a reaction.
     */
    void getFalloff(node_t& f, ReactionData& rdata) {
        string type = f["type"];
        vector<string> p;
        getStringArray(f,p);
        vector_fp c;
        int np = p.size();
        for (int n = 0; n < np; n++) {
            c.push_back(fpValue(p[n]));
        }
        if (type == "Troe") {
            if (np == 4) rdata.falloffType = TROE4_FALLOFF;
            else rdata.falloffType = TROE3_FALLOFF;
        }
        else if (type == "SRI") {
            if (np == 5) rdata.falloffType = SRI5_FALLOFF;
            else rdata.falloffType = SRI3_FALLOFF;
        }
        rdata.falloffParameters = c;
    }

    /**
     * Get the enhanced collision efficiencies. It is assumed that the
     * reaction mechanism is homogeneous, so that all species belong
     * to phase(0) of 'kin'.
     */
    void getEfficiencies(node_t& eff, kinetics_t& kin, ReactionData& rdata) {

        // set the default collision efficiency
        rdata.default_3b_eff = fpValue(eff["default"]);

        vector<string> key, val;
        getPairs(eff, key, val);
        int ne = key.size();
        //map<string, doublereal> e;
        //getFloats(eff, e, false); 
        //map<string, doublereal>::const_iterator bb = e.begin(), ee = e.end();
        string nm;
        string phse = kin.thermo(0).id();
        int n, k;
        for (n = 0; n < ne; n++) { // ; bb != ee; ++bb) {
            nm = key[n];// bb->first;
            k = kin.kineticsSpeciesIndex(nm, phse);
            rdata.thirdBodyEfficiencies[k] = fpValue(val[n]); // bb->second;
        }
    }


    /**
     * Get the rate coefficient for a reaction.
     */
    void getRateCoefficient(node_t& kf, kinetics_t& kin, ReactionData& rdata) {

        int nc = kf.nChildren();
        const nodeset_t& kf_children = kf.children();
        vector_fp clow(3,0.0), chigh(3,0.0);
        int nr = nReacMolecules(rdata);
        for (int m = 0; m < nc; m++) {
            node_t& c = *kf_children[m];
            string nm = c.name();
            int order=0;

            if (nm == "Arrhenius") {
                vector_fp coeff(3);
                getArrhenius(c, order, coeff[0], coeff[1], coeff[2]);
                if (order == 0) order = nr;
                if (order == nr || rdata.reactionType == THREE_BODY_RXN) 
                    chigh = coeff;
                else if (order == nr + 1) clow = coeff;
                else {
                    cerr << "\n\n\n" << endl;
                    kf.write(cerr);
                    throw CanteraError("importCTML",
                        "wrong Arrhenius coeff order");
                }
            }
            else if (nm == "falloff") {
                getFalloff(c, rdata);
            }
            else if (nm == "efficiencies") {
                getEfficiencies(c, kin, rdata);
            }
        }

        if (rdata.reactionType == CHEMACT_RXN) 
            rdata.rateCoeffParameters = clow;        
        else
            rdata.rateCoeffParameters = chigh;

        if (rdata.reactionType == FALLOFF_RXN)
            rdata.auxRateCoeffParameters = clow;
        else if (rdata.reactionType == CHEMACT_RXN) 
            rdata.auxRateCoeffParameters = chigh;        

    }


    /**
     * Create a new ThermoPhase object.
     */
    ThermoPhase* newPhase(XML_Node& xmlphase) {
        XML_Node& th = xmlphase.child("thermo");
        string model = th["model"];
        ThermoPhase* t = newThermoPhase(model);
        importPhase(xmlphase, t);
        return t;
    }


    /**
     * Set the thermodynamic state.
     */
    static void setState(XML_Node& phase, ThermoPhase* th) {
        if (!phase.hasChild("state")) return;
        XML_Node state = phase.child("state");
        doublereal t, p, rho;
        string comp = getString(state,"moleFractions");
        if (comp != "") 
            th->setMoleFractionsByName(comp);
        else {
            comp = getString(state,"massFractions");
            if (comp != "") 
                th->setMassFractionsByName(comp);
        }
        if (state.hasChild("temperature")) {
            t = getFloat(state, "temperature", "temperature");
            th->setTemperature(t);
        }
        if (state.hasChild("pressure")) {
            p = getFloat(state, "pressure", "pressure");
            th->setPressure(p);
        }
        if (state.hasChild("density")) {
            rho = getFloat(state, "density", "density");
            th->setDensity(rho);
        }
    }
        
    /**
     * Import a phase specification.
     */
    bool importPhase(XML_Node& phase, ThermoPhase* th) {

        if (phase.name() != "phase") 
            throw CanteraError("importPhase",
                "Current XML_Node is not a phase element.");

        th->setID(phase.id());                // set the phase id

        // Number of spatial dimensions. Defaults to 3 (bulk phase)
        if (phase.hasAttrib("dim")) {
            int idim = intValue(phase["dim"]);
            if (idim < 1 || idim > 3)
                throw CanteraError("importPhase",
                    "unphysical number of dimensions: "+phase["dim"]);
            th->setNDim(idim);
        }
        else
            th->setNDim(3);     // default


         // equation of state
        if (phase.hasChild("thermo")) {
            XML_Node& eos = phase.child("thermo");
            if (eos["model"] == "Incompressible") {
                if (th->eosType() == cIncompressible) {
                    //map<string, doublereal> d;
                    doublereal rho = getFloat(eos, "density", "-");
                    //doublereal rho = d["density"];
                    th->setParameters(1, &rho);
                }
                else {
                    throw CanteraError("importCTML",
                        "wrong equation of state type");
                }
            }
            else if (eos["model"] == "Surface") {
                if (th->eosType() == cSurf) {
                    map<string, doublereal> d;
                    //getFloats(eos, d);
                    doublereal n = fpValue(eos("site_density"));
                    th->setParameters(1, &n);
                }
                else {
                    throw CanteraError("importCTML",
                        "wrong equation of state type");
                }
            }
        }


        /*************************************************
         *  Add the elements.
         ************************************************/


        // get the declared element names
        XML_Node& elements = phase.child("elementArray");
        vector<string> enames;
        getStringArray(elements, enames);
        
        // // element database defaults to elements.xml
        string element_database; // = "elements.xml";
        if (elements.hasAttrib("datasrc")) 
            element_database = elements["datasrc"];
        XML_Node* db = find_XML(element_database,&phase.root(),"","",
            "elementData");

        int nel = enames.size();
        int i;
        string enm;
        for (i = 0; i < nel; i++) {
            XML_Node* e = db->findByAttr("name",enames[i]);
            if (e) {
                th->addUniqueElement(*e);
            }
            else {
                throw CanteraError("importPhase","no data for element "
                    +enames[i]);
            }
        }
        delete db;
        db = 0;


        /***************************************************************
         * Add the species. First get the speciesArray element, then 
         * the species database.
         ***************************************************************/

        XML_Node& species = phase.child("speciesArray");

        int sprule = 0;
        if (species.hasChild("skip")) {
            XML_Node& sk = species.child("skip");
            string eskip = sk["element"];
            if (eskip == "undeclared") {
                sprule = 1;
            }
        }
            
        db = find_XML(species["datasrc"], &phase.root(), species["idRef"],
            "","speciesData");


        /*******************************************************
         * Set the species thermo manager.
         * Function 'newSpeciesThermoMgr' looks at the species
         * in the database to see what thermodynamic property
         * parameterizations are used, and selects a class
         * that can handle the parameterizations found.
         ******************************************************/

        delete &th->speciesThermo();
        SpeciesThermo* spth = newSpeciesThermoMgr(db);
        th->setSpeciesThermo(spth);
        SpeciesThermo& spthermo = th->speciesThermo();


        /*
         * Get the array of species name strings.
         */                          
        vector<string> spnames;
        getStringArray(species, spnames);
        int nsp = spnames.size();

        map<string,bool> declared;
        string name;
        int k = 0;
        for (i = 0; i < nsp; i++) {
            name = spnames[i];
            
            // Check that every species is only declared once
            if (declared[name]) {
                throw CanteraError("importPhase",
                    "duplicate species: "+name);
            }
            declared[name] = true;

            /*
             * Find the species in the database by name.
             */
            XML_Node* s = db->findByAttr("name",spnames[i]);
            if (s) {
                if (installSpecies(k, *s, *th, spthermo, sprule)) 
                    ++k;
            }
            else {
                throw CanteraError("importPhase","no data for species "
                    +name);
            }
        }
        th->freezeSpecies();
        th->initThermo();
        setState(phase, th);

        th->saveSpeciesData(db);
        return true;
    }        



    bool installReaction(int i, XML_Node& r, Kinetics* k, 
        string default_phase, int rule) {

        Kinetics& kin = *k;
        ReactionData rdata;
        rdata.reactionType = ELEMENTARY_RXN;
        vector_int reac, prod;
        string eqn, type;
        int nn, eqlen;
        vector_fp dummy;

        eqn = r("equation");
        eqlen = eqn.size();
        for (nn = 0; nn < eqlen; nn++) {
            if (eqn[nn] == '[') eqn[nn] = '<';
            if (eqn[nn] == ']') eqn[nn] = '>';
        }

        bool ok;
        // get the reactants
        ok = getReagents(r, kin, 1, default_phase, rdata.reactants, 
            rdata.rstoich, rdata.order, rule);

        // get the products
        ok = ok && getReagents(r, kin, -1, default_phase, rdata.products, 
            rdata.pstoich, dummy, rule);
        if (!ok) {
            cout << "skipping " << eqn << endl;
            return false;
        }

        rdata.equation = eqn;
        rdata.reversible = false;
        rdata.number = i;
        rdata.rxn_number = i;

        string typ = r["type"];
        if (typ == "falloff") {
            rdata.reactionType = FALLOFF_RXN;
            rdata.falloffType = SIMPLE_FALLOFF;
        }
        else if (typ == "chemAct") {
            rdata.reactionType = CHEMACT_RXN;
            rdata.falloffType = SIMPLE_FALLOFF;
        }
        else if (typ == "threeBody") {
            rdata.reactionType = THREE_BODY_RXN;
        }
        else if (typ != "")
            throw CanteraError("installReaction",
                "Unknown reaction type: " + typ);
            
        string isrev = r["reversible"];
        if (isrev == "yes" || isrev == "true")
            rdata.reversible = true;

        getRateCoefficient(r.child("rateCoeff"), kin, rdata);
        kin.addReaction(rdata);
        return true;
    }


    bool installReactionArrays(XML_Node& p, Kinetics& kin, 
        string default_phase) {
        vector<XML_Node*> rarrays;
        int itot = 0;
        p.getChildren("reactionArray",rarrays);
        int na = rarrays.size();
        if (na == 0) return false;
        for (int n = 0; n < na; n++) {
            XML_Node& rxns = *rarrays[n];
            XML_Node* rdata = find_XML(rxns["datasrc"],&rxns.root(),
                "","","reactionData");

            int rxnrule = 0;
            if (rxns.hasChild("skip")) {
                XML_Node& sk = rxns.child("skip");
                string sskip = sk["species"];
                if (sskip == "undeclared") {
                    rxnrule = 1;
                }
            }

            vector<XML_Node*> incl;
            rxns.getChildren("include",incl);
            int ninc = incl.size();
            for (int nii = 0; nii < ninc; nii++) {
                int nrxns = 0;
                XML_Node& ii = *incl[nii];
                vector<string> rxn_ids;
                string pref = ii["prefix"];
                int imin = atoi(ii["min"].c_str());
                int imax = atoi(ii["max"].c_str());
                if (imin != 0 && imax != 0) {
                    nrxns = imax - imin + 1;
                    for (int nn=0; nn<nrxns; nn++) {
                        rxn_ids.push_back(pref+int2str(imin+nn));
                    }
                }
                
                int i;
                for (i = 0; i < nrxns; i++) {
                    XML_Node* r = rdata->findID(rxn_ids[i],1);
                    if (r) {
                        if (installReaction(itot, *r, &kin, 
                                default_phase, rxnrule)) ++itot;
                    }
                }
            }
        }
        kin.finalize();
        return true;
    }


    /**
     * Import a reaction mechanism for a phase or an interface.
     */
    bool importKinetics(XML_Node& phase, vector<ThermoPhase*> th, 
        Kinetics* k) {

        Kinetics& kin = *k;

        // This phase will be the default one
        string default_phase = phase["id"];

        // if other phases are involved in the reaction mechanism,
        // they must be listed in a 'phaseArray' child
        // element. Homogeneous mechanisms do not need to include a
        // phaseArray element.

        vector<string> phase_ids;
        if (phase.hasChild("phaseArray")) {
            XML_Node& pa = phase.child("phaseArray");
            getStringArray(pa, phase_ids);
        }
        phase_ids.push_back(default_phase);
            
        int np = phase_ids.size();

        int nt = th.size();

        // for each referenced phase, attempt to find its id among those
        // phases specified. 
        bool phase_ok;

        string phase_id;
        for (int n = 0; n < np; n++) {
            phase_id = phase_ids[n];
            phase_ok = false;

            // loop over the supplied 'ThermoPhase' objects representing
            // phases, to find an object with the same id.
            for (int m = 0; m < nt; m++) {
                if (th[m]->id() == phase_id) {
                    phase_ok = true;

                    // if no phase with this id has been added to 
                    //the kinetics manager yet, then add this one
                    if (kin.phaseIndex(phase_id) < 0) {
                        kin.addPhase(*th[m]);
                    }
                }
            }
            if (!phase_ok) {
                throw CanteraError("importKinetics",
                    "phase "+phase_id+" not found.");
            }
        }

        // allocates arrays, etc. Must be called after the phases have 
        // been added to 'kin', so that the number of species in each
        // phase is known.
        kin.init();

        // Install the reactions.
        return installReactionArrays(phase, kin, default_phase);
     }

    /**
     * Build a single-phase solution.
     */
    bool buildSolutionFromXML(XML_Node& root, string id, string nm, 
        ThermoPhase* th, Kinetics* k) {
        XML_Node* x;
        x = find_XML("", &root, id, "", nm);
        if (!x) return false;

        importPhase(*x, th);

        vector<ThermoPhase*> phases(1);
        phases[0] = th;
        importKinetics(*x, phases, k);

        return true;
    }
}    
 
