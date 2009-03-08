/**
 *  @file importCTML.cpp
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to Cantera objects. The purpose
 *     of these routines is to intialize the Cantera objects with data
 *     from the ctml tree structures.
 */

/* $Author: dggoodwin $
 * $Revision: 1.48 $
 * $Date: 2006/07/11 16:07:46 $
 */

// Copyright 2002  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "importCTML.h"
#include "mix_defs.h"
#include <time.h>

//   Cantera includes
#include "speciesThermoTypes.h"
#include "ThermoPhase.h"
#include "SurfPhase.h"
#include "EdgePhase.h"
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


// these are all used to check for duplicate reactions
class rxninfo {
public:
    vector< map<int, doublereal> > rdata;
    vector<string> eqn;
    vector<int> dup, nr, typ;
    vector<bool> rev;
};

rxninfo* _rxns = 0;
#define _reactiondata _rxns->rdata
#define _eqn _rxns->eqn
#define _dup _rxns->dup
#define _nr _rxns->nr
#define _typ _rxns->typ
#define _rev _rxns->rev

namespace Cantera {

    /*
     * First we define a couple of typedefs that will
     * be used throught this file
     */
    typedef const vector<XML_Node*>    nodeset_t;
    typedef XML_Node                   node_t;

    const doublereal DefaultPref = 1.01325e5;   // one atm


    /// split a string at a '#' sign. Used to separate a file name
    /// from an id string.
    static void split(const string& src, string& file, string& id) { 
	string::size_type ipound = src.find('#');
        if (ipound != string::npos) {
            id = src.substr(ipound+1,src.size());
            file = src.substr(0,ipound);
        }
        else {
            id = "";
            file = src;
        }
    }


    /**
     * This routine will locate an XML node in either the input
     * XML tree or in another input file specified by the file
     * part of the file_ID string. Searches are based on the
     * ID attribute of the XML element only.
     *
     * @param file_ID This is a concatenation of two strings seperated
     *                by the "#" character. The string before the
     *                pound character is the file name of an xml
     *                file to carry out the search. The string after
     *                the # character is the ID attribute 
     *                of the xml element to search for. 
     *                The string is interpreted as a file string if
     *                no # character is in the string.
     *
     * @param root    If the file string is empty, searches for the
     *                xml element with matching ID attribute are
     *                carried out from this XML node.
     */
    XML_Node* get_XML_Node(const string& file_ID, XML_Node* root) {
        string fname, idstr;
        XML_Node *db, *doc;
        split(file_ID, fname, idstr);
        if (fname == "") {
            if (!root) throw CanteraError("get_XML_Node",
                "no file name given. file_ID = "+file_ID);
            db = root->findID(idstr, 3);
        } 
        else {
            doc = get_XML_File(fname);
            if (!doc) throw CanteraError("get_XML_Node", 
                "get_XML_File failed trying to open "+fname);
            db = doc->findID(idstr, 3);
        }
        if (!db) {
            throw CanteraError("get_XML_Node", 
                "id tag '"+idstr+"' not found.");
        }
        return db;
    }

   /**
     * This routine will locate an XML node in either the input
     * XML tree or in another input file specified by the file
     * part of the file_ID string. Searches are based on the
     * XML element name and the ID attribute of the XML element.
     * An exact match of both is usually required. However, the
     * ID attribute may be set to "", in which case the first
     * xml element with the correct element name will be returned.
     *
     * @param nameTarget This is the XML element name to look for.
     *                   
     * @param file_ID This is a concatenation of two strings seperated
     *                by the "#" character. The string before the
     *                pound character is the file name of an xml
     *                file to carry out the search. The string after
     *                the # character is the ID attribute 
     *                of the xml element to search for. 
     *                The string is interpreted as a file string if
     *                no # character is in the string.
     *
     * @param root    If the file string is empty, searches for the
     *                xml element with matching ID attribute are
     *                carried out from this XML node.
     */
    XML_Node* get_XML_NameID(const string& nameTarget,
			     const string& file_ID, 
			     XML_Node* root) {
        string fname, idTarget;
        XML_Node *db, *doc;
        split(file_ID, fname, idTarget);
        if (fname == "") {
            if (!root) return 0;
            db = root->findNameID(nameTarget, idTarget);
        } else {
            doc = get_XML_File(fname);
            if (!doc) return 0;
            db = doc->findNameID(nameTarget, idTarget);
        }
        return db;
    }


    /**
     * Install a species into a ThermoPhase object, which defines
     * the phase thermodynamics and speciation.
     *
     *  This routine first gathers the information from the Species XML
     *  tree and calls addUniqueSpecies() to add it to the
     *  ThermoPhase object, p.
     *  This information consists of:
     *         ecomp[] = element composition of species.
     *         chgr    = electric charge of species
     *         name    = string name of species
     *         sz      = size of the species 
     *                 (option double used a lot in thermo)
     *
     *  Then, the routine processes the "thermo" XML element and
     *  calls underlying utility routines to read the XML elements
     *  containing the thermodynamic information for the reference
     *  state of the species. Failures or lack of information trigger
     *  an "UnknownSpeciesThermoModel" exception being thrown.
     */
    bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
        SpeciesThermo& spthermo, int rule, SpeciesThermoFactory* factory) {

	// get the composition of the species
	const XML_Node& a = s.child("atomArray");
	map<string,string> comp;
	getMap(a, comp);

	// check that all elements in the species
	// exist in 'p'. If rule != 0, quietly skip 
        // this species if some elements are undeclared;
        // otherwise, throw an exception
	map<string,string>::const_iterator _b = comp.begin();
	for (; _b != comp.end(); ++_b) {
            if (p.elementIndex(_b->first) < 0) {
                if (rule == 0) {
                    throw CanteraError("installSpecies", 
                        "Species " + s["name"] + 
                        " contains undeclared element " + _b->first);
                }
                else
                    return false;
            }
	}

        // construct a vector of atom numbers for each 
        // element in phase p. Elements not declared in the
        // species (i.e., not in map comp) will have zero
        // entries in the vector.
	int m, nel = p.nElements();
	vector_fp ecomp(nel, 0.0);            
	for (m = 0; m < nel; m++) {
            ecomp[m] = atoi(comp[p.elementName(m)].c_str());
	}


        // get the species charge, if any. Note that the charge need
        // not be explicitly specified if special element 'E'
        // (electron) is one of the elements.
	doublereal chrg = 0.0;
	if (s.hasChild("charge")) chrg = getFloat(s, "charge");

        // get the species size, if any. (This is used by surface
        // phases to represent how many sites a species occupies.)
	doublereal sz = 1.0;
	if (s.hasChild("size")) sz = getFloat(s, "size");

        // add the species to phase p.
	p.addUniqueSpecies(s["name"], &ecomp[0], chrg, sz);

        // install the thermo parameterization for this species into
        // the species thermo manager for phase p.
        factory->installThermoForSpecies(k, s, spthermo);
        
	return true;
    }



    /**
     * Check a reaction to see if the elements balance.
     */
    void checkRxnElementBalance(Kinetics& kin, 
        const ReactionData &rdata, doublereal errorTolerance) {
	int index, klocal, n, kp, kr, m, nel;
        double kstoich;

        map<string, double> bal, balr, balp;
        bal.clear();
        balp.clear();
        balr.clear();

	int np = rdata.products.size();

        // iterate over the products
	for (index = 0; index < np; index++) {
            kp = rdata.products[index];     // index of the product in 'kin'
            n = kin.speciesPhaseIndex(kp);  // phase this product belongs to
            klocal = kp - kin.kineticsSpeciesIndex(0,n); // index within this phase
            kstoich = rdata.pstoich[index]; // product stoichiometric coeff
            const ThermoPhase& ph = kin.speciesPhase(kp); 
            nel = ph.nElements();
            for (m = 0; m < nel; m++) {
                bal[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
                balp[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
            }
	}
	int nr = rdata.reactants.size();
	for (index = 0; index < nr; index++) {
            kr = rdata.reactants[index];
            n = kin.speciesPhaseIndex(kr);
            //klocal = kr - kin.start(n);
            klocal = kr - kin.kineticsSpeciesIndex(0,n);
            kstoich = rdata.rstoich[index];
            const ThermoPhase& ph = kin.speciesPhase(kr);
            nel = ph.nElements();
            for (m = 0; m < nel; m++) {
                bal[ph.elementName(m)] -= kstoich*ph.nAtoms(klocal,m);
                balr[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
            }
	}

        map<string, double>::iterator b = bal.begin();
        string msg = "\n\tElement    Reactants    Products";
        bool ok = true;
        doublereal err;
        for (; b != bal.end(); ++b) {
            err = fabs(b->second/(balr[b->first] + balp[b->first]));
            if (err > errorTolerance) {
                ok = false;
                msg += "\n\t"+b->first+"           "+ fp2str(balr[b->first])
                       +"           "+ fp2str(balp[b->first]);
            }
        }
        if (!ok) {
            msg = "The following reaction is unbalanced:\n\t"
                  + rdata.equation + "\n" + msg + "\n";
            throw CanteraError("checkRxnElementBalance",msg);
        }
    }



    /**
     * Get the reactants or products of a reaction. The information
     * is returned in the spnum, stoich, and order vectors. The
     * length of the vectors is the number of different types of
     * reactants or products found for the reaction.
     *
     * Input
     * --------
     *  rxn -> xml node pointing to the reaction element
     *         in the xml tree.
     *  kin -> Reference to the kinetics object to install
     *         the information into.
     *  rp = 1 -> Go get the reactants for a reaction
     *      -1 -> Go get the products for a reaction
     *  default_phase = String name for the default phase
     *          to loop up species in.
     *  Output
     * -----------
     *  spnum = vector of species numbers found. 
     *          Length is number of reactants or products. 
     *  stoich = stoichiometric coefficient of the reactant or product
     *          Length is number of reactants or products. 
     *  order = Order of the reactant and product in the reaction
     *          rate expression
     *  rule = If we fail to find a species, we will throw an error
     *         if rule != 1. If rule = 1, we simply return false, 
     *         allowing the calling routine to skip this reaction
     *         and continue.
     */
    static bool getReagents(const XML_Node& rxn, kinetics_t& kin, int rp,
        string default_phase, 
        vector_int& spnum, vector_fp& stoich, vector_fp& order,
        int rule) {

        string rptype;

	/*
	 * The id of reactants and products are kept in child elements
	 * of reaction, named "reactants" and "products". We search
	 * the xml tree for these children based on the value of rp,
	 * and store the xml element pointer here.
	 */
        if (rp == 1) rptype = "reactants";
        else rptype = "products";
        const XML_Node& rg = rxn.child(rptype);

	/*
	 * The species and stoichiometric coefficient for the species
	 * are stored as a colon seperated pair. Get all of these 
	 * pairs in the reactions/products object.
	 */
        vector<string> key, val;
        getPairs(rg, key, val);

        int ns = static_cast<int>(key.size());

	/*
	 * Loop over each of the pairs and process them
	 */
        int isp;
        doublereal ord, stch;
        string ph, sp;
        map<string, int> speciesMap;
        for (int n = 0; n < ns; n++) {
            sp = key[n]; // sp is the string name for species
            ph = "";
	    /*
	     * Search for the species in the kinetics object using the
	     * member function kineticsSpeciesIndex(). We will search
	     * for the species in all phases defined in the kinetics operator. 
	     */
            isp = kin.kineticsSpeciesIndex(sp,"<any>");
            if (isp < 0) {
                if (rule == 1) 
                    return false;
                else {
                    throw CanteraError("getReagents",
                        "Undeclared reactant or product species "+sp);
                    return false;
                }
            }

	    /*
	     * For each reagent, we store the the species number, isp
	     * the stoichiometric coefficient, val[n], and the order
	     * species in the reaction rate expression. We assume mass
	     * action kinetics here, but will modify this below for
	     * specified species.
	     */
            spnum.push_back(isp);
            stch = atof(val[n].c_str());
            stoich.push_back(stch);
            ord = doublereal(stch);
            order.push_back(ord);
            /*
             * Needed to process reaction orders below.
             */
            speciesMap[sp] = order.size();
        }

        /*
         * Check to see if reactant reaction orders have been specified. 
         */
        if (rp == 1 && rxn.hasChild("order")) {
            vector<XML_Node*> ord;
            rxn.getChildren("order",ord);
            int norder = static_cast<int>(ord.size());
            int loc;
            doublereal forder;
            for (int nn = 0; nn < norder; nn++) {
                const XML_Node& oo = *ord[nn];
                string sp = oo["species"];
                loc = speciesMap[sp];
                if (loc == 0) 
                    throw CanteraError("getReagents",
                        "reaction order specified for non-reactant: "
                        +sp);
                forder = fpValue(oo());
                if (forder < 0.0) {
                    throw CanteraError("getReagents",
                        "reaction order must be non-negative");
                }
                // replace the forward stoichiometric coefficient
                // stored above in 'order' with the specified
                // reaction order
                order[loc-1] = forder;
            }
        }
        return true;
    }
    

    /**
     * getArrhenius() parses the xml element called Arrhenius. 
     * The Arrhenius expression is
     * \f[        k =  A T^(b) exp (-E_a / RT). \f]
     */
    static void getArrhenius(const XML_Node& node, int& highlow, 
        doublereal& A, doublereal& b, doublereal& E) {
        
        if (node["name"] == "k0") 
            highlow = 0;
        else highlow = 1;
	/*
	 * We parse the children for the A, b, and E conponents.
	 */
        A = getFloat(node, "A", "-");
        b = getFloat(node, "b");
        E = getFloat(node, "E", "actEnergy");
        E /= GasConstant;
    }                

    /**
     * getStick() processes the XML element called Stick that specifies
     * the sticking coefficient reaction. This routine will 
     * translate the sticking coefficient value into a "normal"
     * rate constant for the surface reaction.
     *
     *  Output
     * -----------
     * Output is the normal Arrhenius expressions for a surface
     * reaction rate constant.
     * 
     *   A - units such that rate of rxn has kmol/m^2/s when
     *       A is multiplied by activity concentrations of 
     *       reactants in the normal manner.
     *   n - unitless
     *   E - Units 1/Kelvin
     */
    static void getStick(const XML_Node& node, Kinetics& kin,
        ReactionData& r, doublereal& A, doublereal& b, doublereal& E) {
        int nr = r.reactants.size();
        int k, klocal, not_surf = 0;
        int np = 0;
        doublereal f = 1.0;
        doublereal order;
	/*
	 * species is the name of the special reactant whose surface
	 * flux rate will be calculated.
	 *      isp = species # in the local phase
	 *      ispKinetics = species # in the kinetics object
	 *      ispPhaseIndex = phase # of the special species
	 */
        string spname = node["species"];
        ThermoPhase& th = kin.speciesPhase(spname);
	int isp = th.speciesIndex(spname);
	int ispKinetics = kin.kineticsSpeciesIndex(spname);
	int ispPhaseIndex = kin.speciesPhaseIndex(ispKinetics);
  
        double ispMW = th.molecularWeights()[isp];
	double sc;

        // loop over the reactants
        for (int n = 0; n < nr; n++) {
            k = r.reactants[n];
            order = r.order[n];    // stoich coeff

            // get the phase species k belongs to
            np = kin.speciesPhaseIndex(k);
            const ThermoPhase& p = kin.thermo(np);

            // get the local index of species k in this phase
            klocal = p.speciesIndex(kin.kineticsSpeciesName(k));

            // if it is a surface species, divide f by the standard
            // concentration for this species, in order to convert
            // from concentration units used in the law of mass action
            // to coverages used in the sticking probability
            // expression
            if (p.eosType() == cSurf || p.eosType() == cEdge) {
	      sc = p.standardConcentration(klocal);
	      f /= pow(sc, order);
            }   
            // Otherwise:
            else {
	      // We only allow one species to be in the phase
	      // containing the special sticking coefficient
	      // species.
	      if (ispPhaseIndex == np) {
		not_surf++;
	      } 
	      // Other bulk phase species on the other side
	      // of ther interface are treated like surface
	      // species.
	      else {
		sc = p.standardConcentration(klocal);
		f /= pow(sc, order);
	      }
	    }
        }
        if (not_surf != 1) {
            throw CanteraError("getStick",
                "reaction probabilities can only be used in "
                "reactions with exactly 1 gas/liquid species.");
        }

        doublereal cbar = sqrt(8.0*GasConstant/(Pi*ispMW));
        A = 0.25 * getFloat(node, "A", "-") * cbar * f;
        b = getFloat(node, "b") + 0.5;
        E = getFloat(node, "E", "actEnergy");
        E /= GasConstant;
    }                

    static void getCoverageDependence(const node_t& node, 
        thermo_t& surfphase, ReactionData& rdata) {
        vector<XML_Node*> cov;
        node.getChildren("coverage", cov);
        int k, nc = static_cast<int>(cov.size());
        doublereal e;
        string spname;
        if (nc > 0) {
            for (int n = 0; n < nc; n++) {
                const XML_Node& cnode = *cov[n];
                spname = cnode["species"];
                k = surfphase.speciesIndex(spname);
                rdata.cov.push_back(doublereal(k));
                rdata.cov.push_back(getFloat(cnode, "a"));
                rdata.cov.push_back(getFloat(cnode, "m"));
                e = getFloat(cnode, "e", "actEnergy");
                rdata.cov.push_back(e/GasConstant);
            }
        }
    }

    /**
     * Get falloff parameters for a reaction.
     */
    static void getFalloff(const node_t& f, ReactionData& rdata) {
        string type = f["type"];
        vector<string> p;
        getStringArray(f,p);
        vector_fp c;
        int np = static_cast<int>(p.size());
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
    static void getEfficiencies(const node_t& eff, kinetics_t& kin, ReactionData& rdata) {

        // set the default collision efficiency
        rdata.default_3b_eff = fpValue(eff["default"]);

        vector<string> key, val;
        getPairs(eff, key, val);
        int ne = static_cast<int>(key.size());
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
     * Extract the rate coefficient for a reaction from the xml node, kf.
     * kf should point to a XML element named "rateCoeff".
     * rdata is the partially filled ReactionData object for the reaction.
     * This function will fill in more fields in the ReactionData object.
     * 
     */
    void getRateCoefficient(const node_t& kf, kinetics_t& kin, 
        ReactionData& rdata, int negA) {

        int nc = kf.nChildren();
        nodeset_t& kf_children = kf.children();
        vector_fp clow(3,0.0), chigh(3,0.0);
        //        int nr = nReacMolecules(rdata);
        for (int m = 0; m < nc; m++) {
            const node_t& c = *kf_children[m];
            string nm = c.name();
            int highlow=0;

            if (nm == "Arrhenius") {
                vector_fp coeff(3);
                if (c["type"] == "stick") {
                    getStick(c, kin, rdata, coeff[0], coeff[1], coeff[2]);
                    chigh = coeff;
                }
                else {
                    getArrhenius(c, highlow, coeff[0], coeff[1], coeff[2]);
                    if (highlow == 1 || rdata.reactionType == THREE_BODY_RXN 
                        || rdata.reactionType == ELEMENTARY_RXN) 
                        chigh = coeff;
                    else clow = coeff;
                }
                if (rdata.reactionType == SURFACE_RXN) {
                    getCoverageDependence(c, 
                        kin.thermo(kin.surfacePhaseIndex()), rdata);
                }

                if (coeff[0] <= 0.0 && negA == 0) {
                    throw CanteraError("getRateCoefficient", 
                        "negative or zero A coefficient for reaction "+int2str(rdata.number));
                }
            }
            else if (nm == "falloff") {
                getFalloff(c, rdata);
            }
            else if (nm == "efficiencies") {
                getEfficiencies(c, kin, rdata);
            }
            else if (nm == "electrochem") {
                rdata.beta = fpValue(c["beta"]);
            }
        }
	/*
	 * Store the coefficients in the ReactionData object for return
	 * from this function.
	 */
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
     * Create a new ThermoPhase object and initializes it according to
     * the XML tree database.  This routine first looks up the
     * identity of the model for the solution thermodynamics in the
     * model attribute of the thermo child of the xml phase
     * node. Then, it does a string lookup on the model to figure out
     * what ThermoPhase derived class is assigned. It creates a new
     * instance of that class, and then calls importPhase() to
     * populate that class with the correct parameters from the XML
     * tree.
     */
    ThermoPhase* newPhase(XML_Node& xmlphase) {
        const XML_Node& th = xmlphase.child("thermo");
        string model = th["model"];
        ThermoPhase* t = newThermoPhase(model);
        importPhase(xmlphase, t);
        return t;
    }

    ThermoPhase* newPhase(string infile, string id) {
        XML_Node* root = get_XML_File(infile); 
        if (id == "-") id = "";
        XML_Node* x = get_XML_Node(string("#")+id, root);
        if (x) 
            return newPhase(*x);
        else
            return 0;
    }

        
    /**
     * Import a phase specification.
     *   Here we read an XML description of the phase.
     *   We import descriptions of the elements that make up the
     *   species in a phase.
     *   We import information about the species, including their
     *   reference state thermodynamic polynomials. We then freeze
     *   the state of the species, and finally call initThermo()
     *   a member function of the ThermoPhase object to "finish"
     *   the description.
     *
     *
     * @param phase This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param th   Pointer to the ThermoPhase object which will
     *             handle the thermodynamics for this phase.
     *             We initialize part of the Thermophase object
     *             here, especially for those objects which are
     *             part of the Cantera Kernel.
     */
    bool importPhase(XML_Node& phase, ThermoPhase* th, 
        SpeciesThermoFactory* spfactory) {

        // Check the the supplied XML node in fact represents a 
        // phase.
        if (phase.name() != "phase") 
            throw CanteraError("importPhase",
                "Current const XML_Node is not a phase element.");

        // if no species thermo factory was supplied,
        // use the default one. 
        if (!spfactory) 
            spfactory = SpeciesThermoFactory::factory();

        // set the id attribute of the phase to the 'id' attribute 
        // in the XML tree.
        th->setID(phase.id());
        th->setName(phase.id());

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


	
        // set equation of state parameters. The parameters are
        // specific to each subclass of ThermoPhase, so this is done
        // by method setParametersFromXML in each subclass.
        if (phase.hasChild("thermo")) {
            const XML_Node& eos = phase.child("thermo");
            th->setParametersFromXML(eos);
        }


        /***************************************************************
         * Add the elements.
         ***************************************************************/
        th->addElementsFromXML(phase);


        /***************************************************************
         * Add the species. 
         *
         * Species definitions may be imported from multiple
         * sources. For each one, a speciesArray element must be
         * present.
         ***************************************************************/
        XML_Node* db = 0;
        vector<XML_Node*> sparrays;
        phase.getChildren("speciesArray", sparrays);
        int jsp, nspa = static_cast<int>(sparrays.size());
        vector<XML_Node*> dbases;
        vector_int sprule(nspa,0);

        // loop over the speciesArray elements
        for (jsp = 0; jsp < nspa; jsp++) {

            const XML_Node& species = *sparrays[jsp];

            // If the speciesArray element has a child element
            //   <skip element="undeclared"> 
            // then set sprule[jsp] to 1, so
            // that any species with an undeclared element will be
            // quietly skipped when importing species.
            if (species.hasChild("skip")) {
                const XML_Node& sk = species.child("skip");
                string eskip = sk["element"];
                if (eskip == "undeclared") {
                    sprule[jsp] = 1;
                }
                string dskip = sk["species"];
                if (dskip == "duplicate") {
                    sprule[jsp] += 10;
                }
            }

            string fname, idstr;

            // get a pointer to the node containing the species
            // definitions for the species declared in this 
            // speciesArray element. This may be in the local file
            // containing the phase element, or may be in another
            // file.            
            db = get_XML_Node(species["datasrc"], &phase.root());

            // add this node to the list of species database nodes.
            dbases.push_back(db);
        }


        // if the phase has a species thermo manager already installed,
        // delete it since we are adding new species.
        delete &th->speciesThermo();

        // create a new species thermo manager.  Function
        // 'newSpeciesThermoMgr' looks at the species in the database
        // to see what thermodynamic property parameterizations are
        // used, and selects a class that can handle the
        // parameterizations found.
        SpeciesThermo* spth = newSpeciesThermoMgr(dbases);

        // install it in the phase object
        th->setSpeciesThermo(spth);
        SpeciesThermo& spthermo = th->speciesThermo();

        // used to check that each species is declared only once
        map<string,bool> declared;

        int i, k = 0;

        // loop over the species arrays
        for (jsp = 0; jsp < nspa; jsp++) {

            const XML_Node& species = *sparrays[jsp]; 
            db = dbases[jsp];

            // Get the array of species name strings.
            vector<string> spnames;
            getStringArray(species, spnames);
            int nsp = static_cast<int>(spnames.size());

            // if 'all' is specified, then add all species 
            // defined in this database to the phase
            if (nsp == 1 && spnames[0] == "all") {
                vector<XML_Node*> allsp;
                db->getChildren("species",allsp);
                nsp = static_cast<int>(allsp.size());
                spnames.resize(nsp);
                for (int nn = 0; nn < nsp; nn++) {
		  spnames[nn] = (*allsp[nn])["name"];
		}
            }
            else if (nsp == 1 && spnames[0] == "unique") {
                vector<XML_Node*> uniquesp;
                db->getChildren("species",uniquesp);
                nsp = static_cast<int>(uniquesp.size());
                spnames.clear();
                spnames.resize(nsp);
                string spnm;
                for (int nn = 0; nn < nsp; nn++) {
		  spnm = (*uniquesp[nn])["name"];
                  if (!declared[spnm]) spnames[nn] = spnm;
		}
            }

            string name;
            bool skip;
            for (i = 0; i < nsp; i++) {
                name = spnames[i];
                skip = false;
                if (name == "") skip = true;
                // Check that every species is only declared once
                if (declared[name]) {
                    if (sprule[jsp] >= 10) 
                        skip = true;
                    else
                        throw CanteraError("importPhase",
                            "duplicate species: "+name);
                }
                if (!skip) {
                    declared[name] = true;

                    // Find the species in the database by name.
                    XML_Node* s = db->findByAttr("name",spnames[i]);
                    if (s) {
                        if (installSpecies(k, *s, *th, spthermo, sprule[jsp], 
                                spfactory)) 
                            ++k;
                    }
                    else {
                        throw CanteraError("importPhase","no data for species "
                            +name);
                    }
                }
            }
        }

        // done adding species. 
        th->freezeSpecies();

        th->saveSpeciesData(db);

        // perform any required subclass-specific initialization.
        string id = "";
        th->initThermoXML(phase, id);

        return true;
    }        



    /**
     * This function returns true if two reactions are duplicates of
     * one another, and false otherwise.  The input arguments are two
     * maps from species number to stoichiometric coefficient, one for
     * each reaction. The reactions are considered duplicates if their
     * stoichiometric coefficients have the same ratio for all
     * species.
     */
    doublereal isDuplicateReaction(map<int, doublereal>& r1, 
                                   map<int, doublereal>& r2) {
        
        map<int, doublereal>::const_iterator b = r1.begin(), e = r1.end();
        int k1 = b->first;
        doublereal ratio = 0.0;
        if (r1[k1] == 0.0 || r2[k1] == 0.0) goto next;
        ratio = r2[k1]/r1[k1];
        ++b;
        for (; b != e; ++b) {
            k1 = b->first;
            if (r1[k1] == 0.0 || r2[k1] == 0.0) goto next;
            if (fabs(r2[k1]/r1[k1] - ratio) > 1.e-8) 
                goto next;
        }
        return ratio;
next:
        ratio = 0.0;
        b = r1.begin();
        k1 = b->first;
        if (r1[k1] == 0.0 || r2[-k1] == 0.0) return 0.0;
        ratio = r2[-k1]/r1[k1];
        ++b;
        for (; b != e; ++b) {
            k1 = b->first;
            if (r1[k1] == 0.0 || r2[-k1] == 0.0) return 0.0;
            if (fabs(r2[-k1]/r1[k1] - ratio) > 1.e-8) 
                return 0.0;
        }
        return ratio;
    }


    /**
     *  Install an individual reaction into a kinetics manager. The
     *  data for the reaction is in the xml_node r. In other words, r
     *  points directly to a ctml element named "reaction". i refers
     *  to the number id of the reaction in the kinetics object.
     * 
     * @param i Reaction number.
     * @param r XML_Node containing reaction data.
     * @param k Kinetics manager to which reaction will be added.
     * @param default_phase ...
     * @param rule Rule for handling reactions with missing species 
     *             (skip or flag as error)
     * @param validate_rxn If true, check that this reaction is not a 
     *   duplicate of one already entered, and check that the reaction 
     *   balances.
     */
    static bool installReaction(int i, const XML_Node& r, Kinetics* k, 
				string default_phase, int rule,
				bool validate_rxn) {

        Kinetics& kin = *k;
	/*
	 *  We use the ReactionData object to store initial values read
	 *  in from the xml data. Then, when we have collected everything
	 *  we add the reaction to the kinetics object, k, at the end
	 * of the routine. (Someday this may be rewritten to skip building
         * the ReactionData object).
	 */
        ReactionData rdata;
        rdata.reactionType = ELEMENTARY_RXN;  // default
        vector_int reac, prod;
        string eqn, type;
        int nn, eqlen;
        vector_fp dummy;

        // check to see if the reaction is specified to be a duplicate
        // of another reaction, or to allow a negative pre-exponential.
        int dup = 0;
        if (r.hasAttrib("duplicate")) dup = 1;
        int negA = 0;
        if (r.hasAttrib("negative_A")) negA = 1;


 	/*
	 * This seemingly simple expression goes and finds the child element,
	 * "equation". Then it treats all of the contents of the "equation"
	 * as a string, and returns it the variable eqn. We post process
	 * the string to convert [ and ] characters into < and >, which 
	 * cannot be stored in an XML file.
	 */
        if (r.hasChild("equation"))
            eqn = r("equation");
        else
            eqn = "<no equation>";

        eqlen = static_cast<int>(eqn.size());
        for (nn = 0; nn < eqlen; nn++) {
            if (eqn[nn] == '[') eqn[nn] = '<';
            if (eqn[nn] == ']') eqn[nn] = '>';
        }


        bool ok;
        // get the reactants
        ok = getReagents(r, kin, 1, default_phase, rdata.reactants, 
            rdata.rstoich, rdata.order, rule);

        /*
	 * Get the products. We store the id of products in rdata.products
	 */
        ok = ok && getReagents(r, kin, -1, default_phase, rdata.products, 
            rdata.pstoich, dummy, rule);

        // if there was a problem getting either the reactants or the products, 
        // then abort.
        if (!ok) return false;

        // check whether the reaction is specified to be
        // reversible. Default is irreversible.
        rdata.reversible = false;
        string isrev = r["reversible"];
        if (isrev == "yes" || isrev == "true")
            rdata.reversible = true;


        string typ = r["type"];

        /*
         * If reaction orders are specified, then this reaction
         * does not follow mass-action kinetics, and is not 
         * an elementary reaction. So check that it is not reversible,
         * since computing the reverse rate from thermochemistry only
         * works for elementary reactions. Set the type to global,
         * so that kinetics managers will know to process the reaction
         * orders.
         */
        if (r.hasChild("order")) {
            if (rdata.reversible == true) 
                throw CanteraError("installReaction",
                    "reaction orders may only be given for "
                    "irreversible reactions");
            rdata.global = true;
        }


	/*
	 * Seaarch the reaction element for the attribute "type".
	 * If found, then branch on the type, to fill in appropriate
	 * fields in rdata. 
	 */

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
        else if (typ == "surface") {
            rdata.reactionType = SURFACE_RXN;
        }
        else if (typ == "edge") {
            rdata.reactionType = EDGE_RXN;
        }
        else if (typ != "")
            throw CanteraError("installReaction",
                "Unknown reaction type: " + typ);

        /*
         * Look for undeclared duplicate reactions.
         */
        if (validate_rxn) {
            doublereal c = 0.0;

            map<int, doublereal> rxnstoich;
            rxnstoich.clear();
            int nr = rdata.reactants.size();
            for (nn = 0; nn < nr; nn++) {
                rxnstoich[-1 - rdata.reactants[nn]] -= rdata.rstoich[nn];
            }
            int np = rdata.products.size();
            for (nn = 0; nn < np; nn++) {
                rxnstoich[rdata.products[nn]+1] += rdata.pstoich[nn];
            }
            int nrxns = static_cast<int>(_reactiondata.size());
            for (nn = 0; nn < nrxns; nn++) {
                if ((int(rdata.reactants.size()) == _nr[nn]) 
                    && (rdata.reactionType == _typ[nn])) {
                    c = isDuplicateReaction(rxnstoich, _reactiondata[nn]);
                    if (c > 0.0 
                        || (c < 0.0 && rdata.reversible)
                        || (c < 0.0 && _rev[nn])) {
                        if ((!dup || !_dup[nn])) {
                            string msg = string("Undeclared duplicate reactions detected: \n")
                                         +"Reaction "+int2str(nn+1)+": "+_eqn[nn]
                                         +"\nReaction "+int2str(i+1)+": "+eqn+"\n";
                            _reactiondata.clear();
                            _eqn.clear();
                            _rev.clear();
                            _nr.clear();
                            _typ.clear();
                            _dup.clear();
                            throw CanteraError("installReaction",msg);
                        }
                    }
                }
            }
            _dup.push_back(dup);
            _rev.push_back(rdata.reversible);
            _eqn.push_back(eqn);
            _nr.push_back(rdata.reactants.size());
            _typ.push_back(rdata.reactionType);
            _reactiondata.push_back(rxnstoich);
        }
        
        rdata.equation = eqn;
        rdata.number = i;
        rdata.rxn_number = i;

	/*
	 * Read the rate coefficient data from the XML file. Trigger an
	 * exception for negative A unless specifically authorized.
	 */
        getRateCoefficient(r.child("rateCoeff"), kin, rdata, negA);

	/*
	 * Check to see that the elements balance in the reaction.
	 * Throw an error if they don't
	 */        
        if (validate_rxn) {
	  checkRxnElementBalance(kin, rdata);
	}

	/*
	 * Ok we have read everything in about the reaction. Add it
	 * to the kinetics object by calling the Kinetics member function,
	 * addReaction()
	 */
        kin.addReaction(rdata);
        return true;
    }



   /**
     *  Take information from the XML tree, p, about reactions
     *  and install them into the kinetics object, kin. 
     *  default_phase is the default phase to assume when
     *  looking up species.
     *
     *  At this point, p usually refers to the phase xml element.
     *  One of the children of this element is reactionArray,
     *  the element which determines where in the xml file to
     *  look up the reaction rate data pertaining to the phase.
     *
     *  On return, if reaction instantiation goes correctly, return true.
     *  If there is a problem, return false.
     */
    bool installReactionArrays(const XML_Node& p, Kinetics& kin, 
        string default_phase, bool check_for_duplicates) {

        if (_rxns == 0) {
            _rxns = new rxninfo;
        }
        _eqn.clear();
        _dup.clear();
        _nr.clear();
        _typ.clear();
        _reactiondata.clear();
        _rev.clear();

        vector<XML_Node*> rarrays;
        int itot = 0;
	/*
	 * Search the children of the phase element for the
	 * xml element named reactionArray. If we can't find it,
	 * then return signaling having not found any reactions.
	 * Apparently, we allow multiple reactionArray elements here
	 * Each one will be processed sequentially, with the
	 * end result being purely additive.
	 */
        p.getChildren("reactionArray",rarrays);
        int na = static_cast<int>(rarrays.size());
        if (na == 0) return false;
        for (int n = 0; n < na; n++) {
	  /*
	   * Go get a reference to the current xml element, 
	   * reactionArray. We will process this element now.
	   */
	  const XML_Node& rxns = *rarrays[n];
	  /*
	   * The reactionArray element has an attribute called,
	   * datasrc. The value of the attribute is the xml
	   * element comprising the top of the
	   * tree of reactions for the phase.
	   * Find this datasrc element starting with the root
	   * of the current xml node.
	   */
          const XML_Node* rdata = get_XML_Node(rxns["datasrc"], &rxns.root());
	  //const XML_Node* rdata = find_XML(rxns["datasrc"],&rxns.root(),
          //			     "","","reactionData");
	  /*
	   * If the reactionArray element has a child element named
	   * "skip", and if the attribute of skip called "species" has
	   * a value of "undeclared", we will set rxnrule = 1.
	   * rxnrule is passed to the routine that parses each individual
	   * reaction. I believe what this means is that the parser will
	   * skip all reactions containing an undefined species without
	   * throwing an error condition.
	   */
	  int rxnrule = 0;
	  if (rxns.hasChild("skip")) {
	    const XML_Node& sk = rxns.child("skip");
	    string sskip = sk["species"];
	    if (sskip == "undeclared") {
	      rxnrule = 1;
	    }
	  }
	  int i, nrxns = 0;
	  /*
	   * Search for child elements called include. We only include
	   * a reaction if it's tagged by one of the include fields.
	   * Or, we include all reactions if there are no include fields.
	   */
	  vector<XML_Node*> incl;
	  rxns.getChildren("include",incl);
	  int ninc = static_cast<int>(incl.size());

	  vector<XML_Node*> allrxns;
	  rdata->getChildren("reaction",allrxns);
	  nrxns = static_cast<int>(allrxns.size());
	  // if no 'include' directive, then include all reactions
	  if (ninc == 0) {
	    for (i = 0; i < nrxns; i++) {
	      const XML_Node* r = allrxns[i];
	      if (r) {
		if (installReaction(itot, *r, &kin, 
                        default_phase, rxnrule, check_for_duplicates)) ++itot;
	      }
	    }
	  }
	  else {
	    for (int nii = 0; nii < ninc; nii++) {
	      const XML_Node& ii = *incl[nii];
	      string imin = ii["min"];
	      string imax = ii["max"];

              string::size_type iwild = string::npos;
              if (imax == imin) {
                  iwild = imin.find("*");
                  if (iwild != string::npos) {
                      imin = imin.substr(0,iwild);
                      imax = imin;
                  }
              }

	      for (i = 0; i < nrxns; i++) {
		const XML_Node* r = allrxns[i];
		string rxid;
		if (r) {
		  rxid = (*r)["id"];
                  if (iwild != string::npos) {
                      rxid = rxid.substr(0,iwild);
                  }
		  /*
		   * To decide whether the reaction is included or not
		   * we do a lexical min max and operation. This 
		   * sometimes has surprising results.
		   */
		  if ((rxid >= imin) && (rxid <= imax)) {
		    if (installReaction(itot, *r, &kin, 
                            default_phase, rxnrule, check_for_duplicates)) ++itot;
		  }
		}
	      }
	    }
	  }
        }

	/*
	 * Finalize the installation of the kinetics, now that we know
	 * the true number of reactions in the mechanism, itot.
	 */
        kin.finalize();
        //writer = 0;
        _eqn.clear();
        _dup.clear();
        _nr.clear();
        _typ.clear();
        _reactiondata.clear();
        delete _rxns;
        _rxns = 0;
        return true;
    }


    /**
     * Import a reaction mechanism for a phase or an interface.
     *
     * @param phase This is an xml node containing a description
     *              of a phase. Within the phase is a XML element
     *              called reactionArray containing the location
     *              of the description of the reactions that make
     *              up the kinetics object. 
     *              Also within the phase is an XML element called
     *              phaseArray containing a listing of other phases
     *              that participate in the kinetics mechanism.
     *
     * @param th    This is a list of ThermoPhase pointers containing
     *              the phases that participate in the kinetics
     *              reactions. All of the phases must have already
     *              been initialized and formed within Cantera.
     *              However, their pointers should not have been
     *              added to the Kinetics object; this addition
     *              is carried out here.
     *
     * @param k     This is a pointer to the kinetics manager class
     *              that will be initialized with a kinetics 
     *              mechanism.
     */
    bool importKinetics(const XML_Node& phase, vector<ThermoPhase*> th, 
        Kinetics* k) {

        if (k == 0) return false;

        Kinetics& kin = *k;

        // This phase will be the default one
        string default_phase = phase["id"];

        bool check_for_duplicates = false;
        if (phase.parent()->hasChild("validate")) {
            const XML_Node& d = phase.parent()->child("validate");
            if (d["reactions"] == "yes") check_for_duplicates = true;
        }

        // if other phases are involved in the reaction mechanism,
        // they must be listed in a 'phaseArray' child
        // element. Homogeneous mechanisms do not need to include a
        // phaseArray element.

        vector<string> phase_ids;
        if (phase.hasChild("phaseArray")) {
            const XML_Node& pa = phase.child("phaseArray");
            getStringArray(pa, phase_ids);
        }
        phase_ids.push_back(default_phase);
            
        int np = static_cast<int>(phase_ids.size());
        int nt = static_cast<int>(th.size());

        // for each referenced phase, attempt to find its id among those
        // phases specified. 
        bool phase_ok;

        string phase_id;
        string msg = "";
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
                msg += " "+th[m]->id();
            }
            if (!phase_ok) {
                throw CanteraError("importKinetics",
                    "phase "+phase_id+" not found. Supplied phases are:"+msg);
            }
        }

        // allocates arrays, etc. Must be called after the phases have 
        // been added to 'kin', so that the number of species in each
        // phase is known.
        kin.init();

        // Install the reactions.
        return installReactionArrays(phase, kin, default_phase, check_for_duplicates);
     }

    /**
     * Build a single-phase ThermoPhase object with associated kinetics
     * mechanism.
     */
    bool buildSolutionFromXML(XML_Node& root, string id, string nm, 
        ThermoPhase* th, Kinetics* k) {
        XML_Node* x;
        try {

            x = get_XML_NameID(nm, string("#")+id, &root);
            //            x = get_XML_Node(string("#")+id, &root); 

            /*
             * Fill in the ThermoPhase object by querying the
             * const XML_Node tree located at x.
             */
            importPhase(*x, th);
            /*
             * Create a vector of ThermoPhase pointers of length 1
             * having the current th ThermoPhase as the entry.
             */
            vector<ThermoPhase*> phases(1);
            phases[0] = th;
            /*
             * Fill in the kinetics object k, by querying the
             * const XML_Node tree located by x. The source terms and
             * eventually the source term vector will be constructed
             * from the list of ThermoPhases in the vector, phases.
             */
            importKinetics(*x, phases, k);

            return true;
        }
        catch (CanteraError) {
            throw CanteraError("buildSolutionFromXML","error encountered");
            return false;
        }
    }

    /**
     * Search an XML tree for species data.
     *
     *   This utility routine will search the XML tree for the species
     *   named by the string, kname. It will return the XML_Node
     *   pointer.
     *   Failures of any kind return the null pointer.
     */
    const XML_Node *speciesXML_Node(string kname,
				  const XML_Node *phaseSpecies) {
      /*
       * First look at the species database.
       *  -> Look for the subelement "stoichIsMods"
       *     in each of the species SS databases.
       */
      if (!phaseSpecies) return ((const XML_Node *) 0);
      string jname;
      vector<XML_Node*> xspecies;
      phaseSpecies->getChildren("species", xspecies);
      int jj = xspecies.size();
      for (int j = 0; j < jj; j++) {
	const XML_Node& sp = *xspecies[j];
          jname = sp["name"];
          if (jname == kname) {
            return &sp;
          }
      }
      return ((const XML_Node *) 0);
    }



}    
 
