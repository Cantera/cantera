
#include "surfKinetics.h"
#include "ctml.h"
using namespace ctml;

namespace Cantera {


    /**
     * Import a surface reaction mechanism
     */
    void importInterfaceData(SurfacePhase* ph, SurfKinetics* kin, 
        string fname, string id) {

        ifstream f(fname.c_str());
        XML_Node root;
        root.build(f);

        XML_Node* srxns = root.findID(id);
        map<string, double> fmap;
        getFloats(*srxns, fmap);
        ph->setSiteDensity(fmap["site_density"]);
        XML_Node& spset = srxns->child("SpeciesArray");
        vector<XML_Node*> sp;
        spset.getChildren("species",sp);
        int nsp = sp.size();
        int k;
        for (k = 0; k < nsp; k++) {
            XML_Node& s = *sp[k];
            ph->addSpecies(s["name"], atof(s["size"].c_str()));
        }

        vector<XML_Node*> rxns;
        srxns->child("ReactionArray").getChildren("reaction",rxns);
        int nrxns = rxns.size();
        int i, n;
        string phase;
        vector_int rindex, order, rstoich, pindex, pstoich;

        // get bulk phase data
        int kk1 = kin->bulkPhase(0)->nSpecies();
        int kk2 = 0;
        if (kin->bulkPhase(1)) kk2 = kin->bulkPhase(1)->nSpecies();
        vector<XML_Node*> bphase;
        srxns->getChildren("phase", bphase);
        int nbulk = bphase.size();
        string s, t;
        vector<string> phase_id(2,"<none>");
        for (int nb = 0; nb < nbulk; nb++) {
            phase_id[nb] = (*bphase[nb])["id"];
        }
        for (i = 0; i < nrxns; i++) {
            XML_Node& rxn = *rxns[i];
            vector<XML_Node*> reac;
            rxn.getChildren("reactant",reac);
            int nr = reac.size();
            int k;
            for (n = 0; n < nr; n++) {
                XML_Node& r = *reac[n];
                rstoich.push_back(atoi(r["stoich"].c_str()));
                order.push_back(atoi(r["order"].c_str()));
                phase = r["phase"];
                if (phase == phase_id[0]) { 
                    k = kin->bulkPhase(0)->speciesIndex(r["name"]);
                }
                else if (phase == phase_id[1]) {
                    k = kin->bulkPhase(1)->speciesIndex(r["name"]) + kk1;
                }
                else { 
                    k = ph->speciesIndex(r["name"]) + kk1 + kk2;
                }
                rindex.push_back(k);
            }

            vector<XML_Node*> prod;
            rxn.getChildren("product",prod);
            int np = prod.size();
            for (n = 0; n < np; n++) {
                XML_Node& p = *prod[n];
                pstoich.push_back(atoi(p["stoich"].c_str()));
                phase = p["phase"];
                if (phase == phase_id[0]) { 
                    k = kin->bulkPhase(0)->speciesIndex(p["name"]);
                }
                else if (phase == phase_id[1]) {
                    k = kin->bulkPhase(1)->speciesIndex(p["name"]) + kk1;
                }
                else { 
                    k = ph->speciesIndex(p["name"]) + kk1 + kk2;
                }
                pindex.push_back(k);
            }

            XML_Node& rate = rxn.child("rate");
            map<string, doublereal> rp;
            getFloats(rate, rp);
            vector_fp kf(3);
            kf[0] = rp["A"];
            kf[1] = rp["n"];
            kf[2] = rp["E"];

            kin->addReaction(rindex, rstoich, order, pindex, pstoich, kf);
        }
    }

}
