/**
 *   Cantera Fortran interface library. This library of functions is designed
 *   to encapsulate Cantera functionality and make it available for
 *   use in languages and applications other than C++. A set of
 *   library functions is provided that are declared "extern C". All
 *   Cantera objects are stored and referenced by integers - no
 *   pointers are passed to or from the calling application.
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// Cantera includes
#include "ChemEquil.h"
#include "KineticsFactory.h"
#include "transport/TransportFactory.h"
#include "ctml.h"
#include "importCTML.h"
#include "converters/ck2ctml.h"
#include "../../clib/src/Cabinet.h"
#include "InterfaceKinetics.h"

#include "flib_defs.h"

inline XML_Node* _xml(integer* n) {
    return Cabinet<XML_Node>::cabinet()->item(*n);
}

inline ThermoPhase* _fph(integer* n) {
    return Cabinet<ThermoPhase>::cabinet()->item(*n);
}

inline Kinetics* _fkin(integer* n) {
    return Cabinet<Kinetics>::cabinet()->item(*n);
}

inline ThermoPhase* _fth(integer* n) {
    return Cabinet<ThermoPhase>::cabinet()->item(*n);
}

inline Transport* _ftrans(integer* n) {
    return Cabinet<Transport>::cabinet()->item(*n);
}

inline string f2string(const char* s, ftnlen n) {
    return string(s, n);
}
 
/**
 * Exported functions.
 */
extern "C" {

    //--------------- Phase ---------------------//

    integer DLL_EXPORT phase_nelements_(integer* n) {
        return _fph(n)->nElements();
    }

    integer DLL_EXPORT phase_nspecies_(integer* n) {
        return _fph(n)->nSpecies();
    }

    doublereal DLL_EXPORT phase_temperature_(integer* n) {
        return _fph(n)->temperature();
    }

    integer DLL_EXPORT phase_settemperature_(integer* n, doublereal* t) {
        _fph(n)->setTemperature(*t);
        return 0;
    }

    doublereal DLL_EXPORT phase_density_(integer* n) {
        return _fph(n)->density();
    }

    integer DLL_EXPORT phase_setdensity_(integer* n, doublereal* rho) {
        _fph(n)->setDensity(*rho);
        return 0;
    }

    doublereal DLL_EXPORT phase_molardensity_(integer* n) {
        return _fph(n)->molarDensity();
    }

    doublereal DLL_EXPORT phase_meanmolecularweight_(integer* n) {
        return _fph(n)->meanMolecularWeight();
    }

    integer DLL_EXPORT phase_elementindex_(integer* n, char* nm, ftnlen lennm) {
        string elnm = f2string(nm, lennm);
        return _fph(n)->elementIndex(elnm);
    }

    integer DLL_EXPORT phase_speciesindex_(integer* n, char* nm, ftnlen lennm) {
        string spnm = f2string(nm, lennm);
        return _fph(n)->speciesIndex(spnm);
    }

    integer DLL_EXPORT phase_getmolefractions_(integer* n, doublereal* x) {
        _fph(n)->getMoleFractions(x);
        return 0;
    }

    doublereal DLL_EXPORT phase_molefraction_(integer* n, integer* k) {
        return _fph(n)->moleFraction(*k);
    }

    integer DLL_EXPORT phase_getmassfractions_(integer* n, doublereal* y) {
        ThermoPhase* p = _fph(n);
        p->getMassFractions(p->nSpecies(), y);
        return 0;
    } 

    doublereal DLL_EXPORT phase_massfraction_(integer* n, integer* k) {
        return _fph(n)->massFraction(*k);
    }

    integer DLL_EXPORT phase_setmolefractions_(integer* n, double* x, integer* norm) {
        ThermoPhase* p = _fph(n);
        if (*norm) p->setMoleFractions(x);
        else p->setMoleFractions_NoNorm(x);
        return 0;
    }

    integer DLL_EXPORT phase_setmolefractionsbyname_(integer* n, char* x, ftnlen lx) {
        try {
            ThermoPhase* p = _fph(n);
            compositionMap xx;
            int nsp = p->nSpecies();
            for (int nn = 0; nn < nsp; nn++) {
                xx[p->speciesName(nn)] = -1;
            }
            parseCompString(f2string(x, lx), xx);
            p->setMoleFractionsByName(xx);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT phase_setmassfractions_(integer* n, doublereal* y, integer* norm) {
        ThermoPhase* p = _fph(n);
        if (*norm) p->setMassFractions(y);
        else p->setMassFractions_NoNorm(y);
        return 0;
    }

    integer DLL_EXPORT phase_setmassfractionsbyname_(integer* n, char* y, ftnlen leny) {
        try {
            ThermoPhase* p = _fph(n);
            compositionMap yy;
            int nsp = p->nSpecies();
            for (int nn = 0; nn < nsp; nn++) {
                yy[p->speciesName(nn)] = -1;
            }
            parseCompString(f2string(y, leny), yy);
            p->setMassFractionsByName(yy);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT phase_getatomicweights_(integer* n, doublereal* atw) {
        ThermoPhase* p = _fph(n);
        const vector_fp& wt = p->atomicWeights();
        copy(wt.begin(), wt.end(), atw);
        return 0;
    }

    integer DLL_EXPORT phase_getmolecularweights_(integer* n, doublereal* mw) {
        ThermoPhase* p = _fph(n);
        const vector_fp& wt = p->molecularWeights();
        copy(wt.begin(), wt.end(), mw);
        return 0;
    }


    integer DLL_EXPORT phase_getspeciesname_(integer* n, integer* k, char* nm, ftnlen lennm) {
        try {
            string spnm = _fph(n)->speciesName(*k);
            int lout = min(lennm,spnm.size());
            copy(spnm.c_str(), spnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) nm[nn] = '\0';            
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT phase_getelementname_(integer* n, integer* m, char* nm, ftnlen lennm) {
        try {
            string elnm = _fph(n)->elementName(*m);
            int lout = min(lennm,elnm.size());
            copy(elnm.c_str(), elnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) nm[nn] = '\0';
            return 0;
        }
        catch (CanteraError) { return -1; }
    }


    doublereal DLL_EXPORT phase_natoms_(integer* n, integer* k, integer* m) {
        try {
            return _fph(n)->nAtoms(*k,*m);
        }
        catch (CanteraError) { return -1; }
    }

        


    //-------------- Thermo --------------------//


    //    integer DLL_EXPORT th_thermoIndex(char* id) {
    //    return thermo_index(id);
    //}

    integer DLL_EXPORT newthermofromxml_(integer* mxml) {
        try {
            XML_Node* x = _xml(mxml);
            thermo_t* th = newPhase(*x);
            return Cabinet<ThermoPhase>::cabinet()->add(th);
        }
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT th_nspecies_(integer* n) {
        return _fth(n)->nSpecies();
    }

    integer DLL_EXPORT th_eostype_(integer* n) {
        return _fth(n)->eosType();
    }

    doublereal DLL_EXPORT th_enthalpy_mole_(integer* n) {
        try {return _fth(n)->enthalpy_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_intenergy_mole_(integer* n) {
        try {return _fth(n)->intEnergy_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_entropy_mole_(integer* n) {
        try {return _fth(n)->entropy_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_gibbs_mole_(integer* n) {
        try {return _fth(n)->gibbs_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_cp_mole_(integer* n) {
        try {return _fth(n)->cp_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_cv_mole_(integer* n) {
        try {return _fth(n)->cv_mole();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_pressure_(integer* n) {
        try {return _fth(n)->pressure();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_enthalpy_mass_(integer* n) {
        try {return _fth(n)->enthalpy_mass();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_intEnergy_mass_(integer* n) {
        try {return _fth(n)->intEnergy_mass();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_entropy_mass_(integer* n) {
        try {return _fth(n)->entropy_mass();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_gibbs_mass_(integer* n) {
        try {return _fth(n)->gibbs_mass();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_cp_mass_(integer* n) {
        try {return _fth(n)->cp_mass();}
        catch (CanteraError) {return DERR;}
    }

    doublereal DLL_EXPORT th_cv_mass_(integer* n) {
        try {return _fth(n)->cv_mass();}
        catch (CanteraError) {return DERR;}
    }

    integer DLL_EXPORT th_chempotentials_(integer* n, doublereal* murt) {
        thermo_t* thrm = _fth(n);
        thrm->getChemPotentials(murt);
        return 0;
    }

    integer DLL_EXPORT th_setpressure_(integer* n, doublereal* p) {
        try {
            _fth(n)->setPressure(*p);
            return 0; 
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT th_set_hp_(integer* n, doublereal* vals) {
        try { 
            _fth(n)->setState_HP(vals[0],vals[1]);
            return 0; 
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT th_set_uv_(integer* n, doublereal* vals) {
        try { 
            _fth(n)->setState_UV(vals[0],vals[1]);
            return 0; 
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT th_set_sv_(integer* n, doublereal* vals) {
        try { _fth(n)->setState_SV(vals[0],vals[1]);
        return 0; }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT th_set_sp_(integer* n, doublereal* vals) {
        try { 
            _fth(n)->setState_SP(vals[0],vals[1]);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT th_equil_(integer* n, integer* XY) {
        try { 
            equilibrate(*_fth(n), *XY); return 0; 
        }
        catch (CanteraError) {return -1;}
    }

    doublereal DLL_EXPORT th_refpressure_(integer* n) {
        return _fth(n)->refPressure();
    }

    doublereal DLL_EXPORT th_mintemp_(integer* n, integer* k) {
        return _fth(n)->minTemp(*k);
    }

    doublereal DLL_EXPORT th_maxtemp_(integer* n, integer* k) {
        return _fth(n)->maxTemp(*k);
    }


    integer DLL_EXPORT th_getenthalpies_rt_(integer* n, doublereal* h_rt) {
        thermo_t* thrm = _fth(n);
        thrm->getEnthalpy_RT(h_rt);
        return 0;
    }

    integer DLL_EXPORT th_getentropies_r_(integer* n, doublereal* s_r) {
        thermo_t* thrm = _fth(n);
        thrm->getEntropy_R(s_r);
        return 0;
    }

    integer DLL_EXPORT th_getcp_r_(integer* n, integer* lenm, doublereal* cp_r) {
        thermo_t* thrm = _fth(n);
        thrm->getCp_R(cp_r);
        return 0;
    }


    //-------------- Kinetics ------------------//

    integer DLL_EXPORT newkineticsfromxml_(integer* mxml, integer* iphase, 
        integer* neighbor1, integer* neighbor2, integer* neighbor3, 
        integer* neighbor4) {
        try {
            XML_Node* x = _xml(mxml);
            vector<thermo_t*> phases;
            phases.push_back(_fth(iphase));
            if (neighbor1 >= 0) {
                phases.push_back(_fth(neighbor1));
                if (neighbor2 >= 0) {
                    phases.push_back(_fth(neighbor2));
                    if (neighbor3 >= 0) {
                        phases.push_back(_fth(neighbor3));
                        if (neighbor4 >= 0) {
                            phases.push_back(_fth(neighbor4));
                        }
                    }
                }
            }
            Kinetics* kin = newKineticsMgr(*x, phases);
            if (kin)
                return Cabinet<Kinetics>::cabinet()->add(kin);
            else
                return 0;
        }
        catch (CanteraError) { return -1; }
    }

//     integer DLL_EXPORT installRxnArrays_(integer* pxml, integer* ikin, 
//         char* default_phase) {
//         try {
//             XML_Node* p = _xml(pxml);
//             kinetics_t* k = kin(ikin);
//             string defphase = string(default_phase);
//             installReactionArrays(*p, *k, defphase);
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//     }

    //-------------------------------------
    integer DLL_EXPORT kin_type_(integer* n) {
        return _fkin(n)->type();
    }

    integer DLL_EXPORT kin_start_(integer* n, integer* p) {
        return _fkin(n)->start(*p);
    }

    integer DLL_EXPORT kin_speciesindex_(integer* n, const char* nm, const char* ph, 
        ftnlen lennm, ftnlen lenph) {
        return _fkin(n)->kineticsSpeciesIndex(f2string(nm, lennm), f2string(ph, lenph));
    }

    //---------------------------------------

    integer DLL_EXPORT kin_ntotalspecies_(integer* n) {
        return _fkin(n)->nTotalSpecies();
    }

    integer DLL_EXPORT kin_nreactions_(integer* n) {
        return _fkin(n)->nReactions();
    }

    doublereal DLL_EXPORT kin_reactantstoichcoeff_(integer* n, integer* k, integer* i) {
        return _fkin(n)->reactantStoichCoeff(*k,*i);
    }

    doublereal DLL_EXPORT kin_productstoichcoeff_(integer* n, integer* k, integer* i) {
        return _fkin(n)->productStoichCoeff(*k,*i);
    }

    integer DLL_EXPORT kin_reactiontype_(integer* n, integer* i) {
        return _fkin(n)->reactionType(*i);
    }

    integer DLL_EXPORT kin_getfwdratesofprogress_(integer* n, doublereal* fwdROP) {
        Kinetics* k = _fkin(n);
        try {
            k->getFwdRatesOfProgress(fwdROP);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_getrevratesofprogress_(integer* n, doublereal* revROP) {
        Kinetics* k = _fkin(n);
        try {
            k->getRevRatesOfProgress(revROP);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_isreversible_(integer* n, integer* i) {
        return (int)_fkin(n)->isReversible(*i);
    }

    integer DLL_EXPORT kin_getnetratesofprogress_(integer* n, doublereal* netROP) {
        try {
            Kinetics* k = _fkin(n);
            k->getNetRatesOfProgress(netROP);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_getcreationrates_(integer* n, doublereal* cdot) {
        try {
            Kinetics* k = _fkin(n);
            k->getCreationRates(cdot);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_getdestructionrates_(integer* n, doublereal* ddot) {
        try {
            Kinetics* k = _fkin(n);
            k->getDestructionRates(ddot);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_getnetproductionrates_(integer* n, doublereal* wdot) {
        try {
            Kinetics* k = _fkin(n);
            k->getNetProductionRates(wdot);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    doublereal DLL_EXPORT kin_multiplier_(integer* n, integer* i) {
        return _fkin(n)->multiplier(*i);
    }

    //integer DLL_EXPORT kin_phase_(integer* n, integer* i) {
    //    return thermo_index(_fkin(n)->thermo(*i).id());
    //}

    integer DLL_EXPORT kin_getequilibriumconstants_(integer* n, doublereal* kc) {
        try {
            Kinetics* k = _fkin(n);
            k->getEquilibriumConstants(kc);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_getreactionstring_(integer* n, integer* i, char* buf, ftnlen lenbuf) {
        try {
            Kinetics* k = _fkin(n);
            string r = k->reactionString(*i);
            int lout = min(lenbuf,r.size());
            copy(r.c_str(), r.c_str() + lout, buf);
            for (int nn = lout; nn < lenbuf; nn++) buf[nn] = '\0';
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_setmultiplier_(integer* n, integer* i, doublereal* v) {
        try {
            _fkin(n)->setMultiplier(*i,*v);
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    integer DLL_EXPORT kin_advancecoverages_(integer* n, doublereal* tstep) {
        try {
            Kinetics* k = _fkin(n);
            if (k->type() == cInterfaceKinetics) {
                ((InterfaceKinetics*)k)->advanceCoverages(*tstep);
            }
            else {
                throw CanteraError("kin_advanceCoverages",
                    "wrong kinetics manager type");
            }
            return 0;
        }
        catch (CanteraError) {return -1;}
    }

    //------------------- Transport ---------------------------

    integer DLL_EXPORT newtransport(char* model,  
        integer* ith, integer* loglevel, ftnlen lenmodel) {
        string mstr = f2string(model, lenmodel);
        thermo_t* t = _fth(ith);
        try {
            Transport* tr = newTransportMgr(mstr, t, *loglevel);
            return Cabinet<Transport>::cabinet()->add(tr);
        }
        catch (CanteraError) { return -1; }
    }
    
    doublereal DLL_EXPORT trans_viscosity_(integer* n) {
        try {return _ftrans(n)->viscosity();}
        catch (CanteraError) { return -1.0; }
    }

    doublereal DLL_EXPORT trans_thermalConductivity_(integer* n) {
        try {return _ftrans(n)->thermalConductivity();}
        catch (CanteraError) { return -1.0; }
    }

    integer DLL_EXPORT trans_getThermalDiffCoeffs_(integer* n, doublereal* dt) {
        try { _ftrans(n)->getThermalDiffCoeffs(dt); return 0; }
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT trans_getMixDiffCoeffs_(integer* n, doublereal* d) {
        try { _ftrans(n)->getMixDiffCoeffs(d); return 0;}
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT trans_getBinDiffCoeffs_(integer* n, integer* ld, doublereal* d) {
        try { _ftrans(n)->getBinaryDiffCoeffs(*ld,d); return 0;}
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT trans_getMultiDiffCoeffs_(integer* n, integer* ld, doublereal* d) {
        try { _ftrans(n)->getMultiDiffCoeffs(*ld,d); return 0;}
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT trans_setParameters_(integer* n, integer* type, integer* k, doublereal* d) {
        try { _ftrans(n)->setParameters(*type, *k, d); return 0;}
        catch (CanteraError) { return -1; }
    }

    //-------------------- Functions ---------------------------

//     integer DLL_EXPORT import_phase_(integer* nth, integer* nxml, char* id, ftnlen lenid) {
//         thermo_t* thrm = th(nth);
//         XML_Node* node = _xml(nxml);
//         string idstr = f2string(id, lenid);
//         try {
//             importPhase(*node, thrm);
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//     }

//     integer DLL_EXPORT import_kinetics_(integer* nxml, char* id, 
//         integer* nphases, integer* ith, integer* nkin, ftnlen lenid) {
//         vector<thermo_t*> phases;
//         for (int i = 0; i < nphases; i++) {
//             phases.push_back(th(ith[i]));
//         }
//         XML_Node* node = _xml(nxml);
//         Kinetics* k = kin(nkin);
//         string idstr = f2string(id, lenid);
//         try {
//             importKinetics(*node, phases, k);
//             return 0;
//         }
//         catch (CanteraError) { return -1; }
//     }


    integer DLL_EXPORT phase_report_(integer* nth, 
        char* buf, integer* show_thermo, ftnlen buflen) {
        try {
            bool stherm = (*show_thermo != 0);
            string s = report(*_fth(nth), stherm);
            if (int(s.size()) > buflen - 1) {
                return -(s.size() + 1);
            }
            copy(s.begin(), s.end(), buf);
            for (int nn = s.size(); nn < buflen; nn++) buf[nn] = '\0';
            return 0;
            
        }
        catch (CanteraError) { return -1; }
    }

    integer DLL_EXPORT getCanteraError_(char* buf, ftnlen buflen) {
        string e; // = "<no error>";
        //if (nErrors() > 0)
        e = lastErrorMessage();
        int n = min(e.size(), buflen-1);
        copy(e.begin(), e.begin() + n, buf);
        for (int nn = n; nn < buflen; nn++) buf[nn] = '\0';
        return 0;
    }

    integer DLL_EXPORT addCanteraDirectory_(integer* buflen, char* buf) {
        addDirectory(string(buf));
        return 0;
    }

//     integer DLL_EXPORT readlog_(integer* n, char* buf) {
//         string s;
//         writelog("function readlog is deprecated!");
//         //getlog(s);
//         int nlog = s.size();
//         if (n < 0) return nlog; 
//         int nn = min(n-1, nlog);
//         copy(s.begin(), s.begin() + nn,
//             buf);
//         buf[min(nlog, n-1)] = '\0';
//         //clearlog();
//         return 0;

//     } 

    integer DLL_EXPORT delThermo_(integer* n) {
        try {
            Cabinet<ThermoPhase>::cabinet()->del(*n);
            return 0;
        }
        catch (CanteraError) {
            return -1;
        }
    }

    integer DLL_EXPORT delKinetics_(integer* n) {
        Cabinet<Kinetics>::cabinet()->del(*n);
        return 0;
    }

    integer DLL_EXPORT delTransport_(integer* n) {
        Cabinet<Transport>::cabinet()->del(*n);
        return 0;
    }

    integer DLL_EXPORT buildSolutionFromXML(char* src, integer* ixml, char* id, 
        integer* ith, integer* ikin, ftnlen lensrc, ftnlen lenid) {

        XML_Node* root = 0;
        if (*ixml > 0) root = _xml(ixml);

        thermo_t* t = _fth(ith);
        kinetics_t* k = _fkin(ikin);

        Kinetics& kin = *k;
        XML_Node *x, *r=0;
        if (root) r = &root->root();
        x = find_XML(f2string(src, lensrc), r, f2string(id,lenid), "", "phase");
        if (!x) return false;
        importPhase(*x, t);
        kin.addPhase(*t);
        kin.init();
        installReactionArrays(*x, kin, x->id());
        t->setState_TP(300.0, OneAtm);
        if (r) {
            if (&x->root() != &r->root()) delete &x->root();
        }
        else delete &x->root();
        return 0;
    }


 //    integer DLL_EXPORT ck_to_ctml(char* in_file, char* db_file,
//         char* tr_file, char* out_file, char* id_tag) {
//         return convert_ck(in_file, db_file, tr_file, out_file, id_tag);
//     }

}
