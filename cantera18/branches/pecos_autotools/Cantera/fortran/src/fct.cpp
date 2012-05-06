/**
 *   Cantera Fortran interface library. This library of functions is designed
 *   to encapsulate Cantera functionality and make it available for
 *   use in languages and applications other than C++. A set of
 *   library functions is provided that are declared "extern C". All
 *   Cantera objects are stored and referenced by integers - no
 *   pointers are passed to or from the calling application.
 */
/*
 * $Id: fct.cpp,v 1.14 2009/07/23 16:56:48 hkmoffa Exp $
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// Cantera includes
#include "equil.h"
#include "KineticsFactory.h"
#include "TransportFactory.h"
#include "ThermoFactory.h"
#include "ctml.h"
#include "importKinetics.h"
#include "../../clib/src/Storage.h"
#include "../../clib/src/Cabinet.h"
#include "InterfaceKinetics.h"
#include "PureFluidPhase.h"

#include "flib_defs.h"

inline XML_Node* _xml(const integer* n) {
    return Cabinet<XML_Node>::cabinet()->item(*n);
}

inline ThermoPhase* _fph(const integer* n) {
    return th(*n);
}

//inline Kinetics* _fkin(const integer* n) {
//    return kin(*n);
//}

static Kinetics* _fkin(const integer* n) {
    if (*n >= 0) 
        return kin(*n);
    else {
        error("_fkin: negative kinetics index");
        return kin(0);
    }
}

inline ThermoPhase* _fth(const integer* n) {
    return th(*n);
}

inline Transport* _ftrans(const integer* n) {
    return trans(*n);
}

std::string f2string(const char* s, ftnlen n) {
    int k;
    std::string ss = "";
    for (k = 0; k < n; k++) {
        if (s[k] == '\0') break;
        ss += s[k];
    }
    return ss;
}

static void handleError() {
    error(lastErrorMessage());
}
 
/**
 * Exported functions.
 */
extern "C" {

    status_t cantera_error_(const char* proc, const char* msg, 
        ftnlen proclen, ftnlen msglen) {
        try {
            std::string sproc = f2string(proc, proclen);
            std::string smsg = f2string(msg, msglen);
            throw CanteraError(sproc, smsg);
        }
        catch (CanteraError) {
            handleError(); return -1;
        }
        return -1;
    }

    //--------------- Phase ---------------------//

    status_t DLL_EXPORT phase_getname_(const integer* n, char* nm, 
        ftnlen lennm) {
        try {
            std::string pnm = _fph(n)->name();
            int lout = min(lennm,pnm.size());
            copy(pnm.c_str(), pnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) nm[nn] = ' ';            
            return 0;
        }
        catch (CanteraError) { handleError(); return -1; }
    }        

    integer DLL_EXPORT phase_nelements_(const integer* n) {
        return _fph(n)->nElements();
    }

    integer DLL_EXPORT phase_nspecies_(const integer* n) {
        return _fph(n)->nSpecies();
    }

    doublereal DLL_EXPORT phase_temperature_(const integer* n) {
        return _fph(n)->temperature();
    }

    status_t DLL_EXPORT phase_settemperature_(const integer* n, doublereal* t) {
        _fph(n)->setTemperature(*t);
        return 0;
    }

    doublereal DLL_EXPORT phase_density_(const integer* n) {
        return _fph(n)->density();
    }

    status_t DLL_EXPORT phase_setdensity_(const integer* n, doublereal* rho) {
        _fph(n)->setDensity(*rho);
        return 0;
    }

    doublereal DLL_EXPORT phase_molardensity_(const integer* n) {
        return _fph(n)->molarDensity();
    }

    doublereal DLL_EXPORT phase_meanmolecularweight_(const integer* n) {
        return _fph(n)->meanMolecularWeight();
    }

    integer DLL_EXPORT phase_elementindex_(const integer* n, char* nm, ftnlen lennm) {
        std::string elnm = f2string(nm, lennm);
        return _fph(n)->elementIndex(elnm) + 1;
    }

    integer DLL_EXPORT phase_speciesindex_(const integer* n, char* nm, ftnlen lennm) {
        std::string spnm = f2string(nm, lennm);
        return _fph(n)->speciesIndex(spnm) + 1;
    }

    status_t DLL_EXPORT phase_getmolefractions_(const integer* n, doublereal* x) {
        _fph(n)->getMoleFractions(x);
        return 0;
    }

    doublereal DLL_EXPORT phase_molefraction_(const integer* n, integer* k) {
        return _fph(n)->moleFraction(*k-1);
    }

    status_t DLL_EXPORT phase_getmassfractions_(const integer* n, doublereal* y) {
        ThermoPhase* p = _fph(n);
        p->getMassFractions(y);
        return 0;
    } 

    doublereal DLL_EXPORT phase_massfraction_(const integer* n, integer* k) {
        return _fph(n)->massFraction(*k-1);
    }

    status_t DLL_EXPORT phase_setmolefractions_(const integer* n, double* x, const integer* norm) {
        ThermoPhase* p = _fph(n);
        if (*norm) p->setMoleFractions(x);
        else p->setMoleFractions_NoNorm(x);
        return 0;
    }

    status_t DLL_EXPORT phase_setmolefractionsbyname_(const integer* n, char* x, ftnlen lx) {
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
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT phase_setmassfractions_(const integer* n, doublereal* y, const integer* norm) {
        ThermoPhase* p = _fph(n);
        if (*norm) p->setMassFractions(y);
        else p->setMassFractions_NoNorm(y);
        return 0;
    }

    status_t DLL_EXPORT phase_setmassfractionsbyname_(const integer* n, char* y, ftnlen leny) {
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
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT phase_getatomicweights_(const integer* n, doublereal* atw) {
        ThermoPhase* p = _fph(n);
        const vector_fp& wt = p->atomicWeights();
        copy(wt.begin(), wt.end(), atw);
        return 0;
    }

    status_t DLL_EXPORT phase_getmolecularweights_(const integer* n, doublereal* mw) {
        ThermoPhase* p = _fph(n);
        const vector_fp& wt = p->molecularWeights();
        copy(wt.begin(), wt.end(), mw);
        return 0;
    }


    status_t DLL_EXPORT phase_getspeciesname_(const integer* n, integer* k, char* nm, ftnlen lennm) {
        try {
            std::string spnm = _fph(n)->speciesName(*k-1);
            int lout = min(lennm,spnm.size());
            copy(spnm.c_str(), spnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) nm[nn] = ' ';            
            return 0;
        }
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT phase_getelementname_(const integer* n, integer* m, char* nm, ftnlen lennm) {
        try {
            std::string elnm = _fph(n)->elementName(*m-1);
            int lout = min(lennm,elnm.size());
            copy(elnm.c_str(), elnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) nm[nn] = ' ';
            return 0;
        }
        catch (CanteraError) { handleError(); return -1; }
    }


    doublereal DLL_EXPORT phase_natoms_(const integer* n, integer* k, integer* m) {
        try {
            return _fph(n)->nAtoms(*k-1,*m-1);
        }
        catch (CanteraError) { handleError(); return -1; }
    }

        


    //-------------- Thermo --------------------//


    //    status_t DLL_EXPORT th_thermoIndex(char* id) {
    //    return thermo_index(id);
    //}

    integer DLL_EXPORT newthermofromxml_(integer* mxml) {
        try {
            XML_Node* x = _xml(mxml);
            thermo_t* th = newPhase(*x);
            return Storage::storage()->addThermo(th);
        }
        catch (CanteraError) { handleError(); return -1; }
    }

    integer DLL_EXPORT th_nspecies_(const integer* n) {
        return _fth(n)->nSpecies();
    }

    integer DLL_EXPORT th_eostype_(const integer* n) {
        return _fth(n)->eosType();
    }

    doublereal DLL_EXPORT th_enthalpy_mole_(const integer* n) {
        try {return _fth(n)->enthalpy_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_intenergy_mole_(const integer* n) {
        try {return _fth(n)->intEnergy_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_entropy_mole_(const integer* n) {
        try {return _fth(n)->entropy_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_gibbs_mole_(const integer* n) {
        try {return _fth(n)->gibbs_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_cp_mole_(const integer* n) {
        try {return _fth(n)->cp_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_cv_mole_(const integer* n) {
        try {return _fth(n)->cv_mole();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_pressure_(const integer* n) {
        try {return _fth(n)->pressure();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_enthalpy_mass_(const integer* n) {
        try {return _fth(n)->enthalpy_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_intenergy_mass_(const integer* n) {
        try {return _fth(n)->intEnergy_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_entropy_mass_(const integer* n) {
        try {return _fth(n)->entropy_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_gibbs_mass_(const integer* n) {
        try {return _fth(n)->gibbs_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_cp_mass_(const integer* n) {
        try {return _fth(n)->cp_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    doublereal DLL_EXPORT th_cv_mass_(const integer* n) {
        try {return _fth(n)->cv_mass();}
        catch (CanteraError) {handleError(); return DERR;}
    }

    status_t DLL_EXPORT th_chempotentials_(const integer* n, doublereal* murt) {
        thermo_t* thrm = _fth(n);
        thrm->getChemPotentials(murt);
        return 0;
    }

    status_t DLL_EXPORT th_setpressure_(const integer* n, doublereal* p) {
        try {
            _fth(n)->setPressure(*p);
            return 0; 
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT th_set_hp_(const integer* n, doublereal* v1, doublereal* v2) {
        try { 
            _fth(n)->setState_HP(*v1, *v2);
            return 0; 
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT th_set_uv_(const integer* n, doublereal* v1, doublereal* v2) {
        try { 
            _fth(n)->setState_UV(*v1, *v2);
            return 0; 
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT th_set_sv_(const integer* n, doublereal* v1, doublereal* v2) {
        try { _fth(n)->setState_SV(*v1, *v2);
        return 0; }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT th_set_sp_(const integer* n, doublereal* v1, doublereal* v2) {
        try { 
            _fth(n)->setState_SP(*v1, *v2);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT th_equil_(const integer* n, char* XY, ftnlen lenxy) {
        try { 
            equilibrate(*_fth(n), f2string(XY,lenxy).c_str()); return 0; 
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    doublereal DLL_EXPORT th_refpressure_(const integer* n) {
        return _fth(n)->refPressure();
    }

    doublereal DLL_EXPORT th_mintemp_(const integer* n, integer* k) {
        return _fth(n)->minTemp(*k-1);
    }

    doublereal DLL_EXPORT th_maxtemp_(const integer* n, integer* k) {
        return _fth(n)->maxTemp(*k-1);
    }


    status_t DLL_EXPORT th_getenthalpies_rt_(const integer* n, doublereal* h_rt) {
        thermo_t* thrm = _fth(n);
        thrm->getEnthalpy_RT(h_rt);
        return 0;
    }

    status_t DLL_EXPORT th_getentropies_r_(const integer* n, doublereal* s_r) {
        thermo_t* thrm = _fth(n);
        thrm->getEntropy_R(s_r);
        return 0;
    }

    status_t DLL_EXPORT th_getcp_r_(const integer* n, integer* lenm, doublereal* cp_r) {
        thermo_t* thrm = _fth(n);
        thrm->getCp_R(cp_r);
        return 0;
    }


    //-------------- Kinetics ------------------//

    integer DLL_EXPORT newkineticsfromxml_(integer* mxml, integer* iphase, 
        const integer* neighbor1, const integer* neighbor2, const integer* neighbor3, 
        const integer* neighbor4) {
        try {
            XML_Node* x = _xml(mxml);
            vector<thermo_t*> phases;
            phases.push_back(_fth(iphase));
            if (*neighbor1 >= 0) {
                phases.push_back(_fth(neighbor1));
                if (*neighbor2 >= 0) {
                    phases.push_back(_fth(neighbor2));
                    if (*neighbor3 >= 0) {
                        phases.push_back(_fth(neighbor3));
                        if (*neighbor4 >= 0) {
                            phases.push_back(_fth(neighbor4));
                        }
                    }
                }
            }
            Kinetics* kin = newKineticsMgr(*x, phases);
            if (kin) {
                int k = Storage::storage()->addKinetics(kin);
                return k; //Storage::storage()->addKinetics(kin);
            }
            else {
                return 0;
            }
        }
        catch (CanteraError) { handleError(); return 999; }
    }

//     status_t DLL_EXPORT installRxnArrays_(integer* pxml, integer* ikin, 
//         char* default_phase) {
//         try {
//             XML_Node* p = _xml(pxml);
//             kinetics_t* k = kin(ikin);
//             string defphase = string(default_phase);
//             installReactionArrays(*p, *k, defphase);
//             return 0;
//         }
//         catch (CanteraError) { handleError(); return -1; }
//     }

    //-------------------------------------
    integer DLL_EXPORT kin_type_(const integer* n) {
        return _fkin(n)->type();
    }

    integer DLL_EXPORT kin_start_(const integer* n, integer* p) {
        return _fkin(n)->start(*p)+1;
    }

    integer DLL_EXPORT kin_speciesindex_(const integer* n, const char* nm, const char* ph, 
        ftnlen lennm, ftnlen lenph) {
        return _fkin(n)->kineticsSpeciesIndex(f2string(nm, lennm), f2string(ph, lenph))+1;
    }

    //---------------------------------------

    integer DLL_EXPORT kin_ntotalspecies_(const integer* n) {
        return _fkin(n)->nTotalSpecies();
    }

    integer DLL_EXPORT kin_nreactions_(const integer* n) {
        return _fkin(n)->nReactions();
    }

    integer DLL_EXPORT kin_nphases_(const integer* n) {
        return _fkin(n)->nPhases();
    }

    integer DLL_EXPORT kin_phaseindex_(const integer* n, const char* ph,
        ftnlen lenph) {
        return _fkin(n)->phaseIndex(f2string(ph, lenph));
    }

    doublereal DLL_EXPORT kin_reactantstoichcoeff_(const integer* n, integer* k, integer* i) {
        return _fkin(n)->reactantStoichCoeff(*k-1,*i-1);
    }

    doublereal DLL_EXPORT kin_productstoichcoeff_(const integer* n, integer* k, integer* i) {
        return _fkin(n)->productStoichCoeff(*k-1,*i-1);
    }

    integer DLL_EXPORT kin_reactiontype_(const integer* n, integer* i) {
        return _fkin(n)->reactionType(*i-1);
    }

    status_t DLL_EXPORT kin_getfwdratesofprogress_(const integer* n, doublereal* fwdROP) {
        Kinetics* k = _fkin(n);
        try {
            k->getFwdRatesOfProgress(fwdROP);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_getrevratesofprogress_(const integer* n, doublereal* revROP) {
        Kinetics* k = _fkin(n);
        try {
            k->getRevRatesOfProgress(revROP);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    integer DLL_EXPORT kin_isreversible_(const integer* n, integer* i) {
        return (int)_fkin(n)->isReversible(*i);
    }

    status_t DLL_EXPORT kin_getnetratesofprogress_(const integer* n, doublereal* netROP) {
        try {
            Kinetics* k = _fkin(n);
            k->getNetRatesOfProgress(netROP);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_getcreationrates_(const integer* n, doublereal* cdot) {
        try {
            Kinetics* k = _fkin(n);
            k->getCreationRates(cdot);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_getdestructionrates_(const integer* n, doublereal* ddot) {
        try {
            Kinetics* k = _fkin(n);
            k->getDestructionRates(ddot);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_getnetproductionrates_(const integer* n, doublereal* wdot) {
        try {
            Kinetics* k = _fkin(n);
            k->getNetProductionRates(wdot);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    doublereal DLL_EXPORT kin_multiplier_(const integer* n, integer* i) {
        return _fkin(n)->multiplier(*i);
    }

    //status_t DLL_EXPORT kin_phase_(const integer* n, integer* i) {
    //    return thermo_index(_fkin(n)->thermo(*i).id());
    //}

    status_t DLL_EXPORT kin_getequilibriumconstants_(const integer* n, doublereal* kc) {
        try {
            Kinetics* k = _fkin(n);
            k->getEquilibriumConstants(kc);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_getreactionstring_(const integer* n, integer* i, char* buf, ftnlen lenbuf) {
        try {
            Kinetics* k = _fkin(n);
            std::string r = k->reactionString(*i-1);
            int lout = min(lenbuf,r.size());
            copy(r.c_str(), r.c_str() + lout, buf);
            for (int nn = lout; nn < lenbuf; nn++) buf[nn] = ' ';
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_setmultiplier_(const integer* n, integer* i, doublereal* v) {
        try {
            _fkin(n)->setMultiplier(*i-1,*v);
            return 0;
        }
        catch (CanteraError) {handleError(); return -1;}
    }

    status_t DLL_EXPORT kin_advancecoverages_(const integer* n, doublereal* tstep) {
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
        catch (CanteraError) {handleError(); return -1;}
    }

    //------------------- Transport ---------------------------

    integer DLL_EXPORT newtransport_(char* model,  
        integer* ith, integer* loglevel, ftnlen lenmodel) {
        std::string mstr = f2string(model, lenmodel);
        thermo_t* t = _fth(ith);
        try {
            Transport* tr = newTransportMgr(mstr, t, *loglevel);
            return Storage::storage()->addTransport(tr);
        }
        catch (CanteraError) { handleError(); return -1; }
    }
    
    doublereal DLL_EXPORT trans_viscosity_(const integer* n) {
        try {return _ftrans(n)->viscosity();}
        catch (CanteraError) { handleError(); return DERR; }
    }

    doublereal DLL_EXPORT trans_thermalconductivity_(const integer* n) {
        try {return _ftrans(n)->thermalConductivity();}
        catch (CanteraError) { handleError(); return DERR; }
    }

    status_t DLL_EXPORT trans_getthermaldiffcoeffs_(const integer* n, doublereal* dt) {
        try { _ftrans(n)->getThermalDiffCoeffs(dt); return 0; }
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT trans_getmixdiffcoeffs_(const integer* n, doublereal* d) {
        try { _ftrans(n)->getMixDiffCoeffs(d); return 0;}
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT trans_getbindiffcoeffs_(const integer* n, integer* ld, doublereal* d) {
        try { _ftrans(n)->getBinaryDiffCoeffs(*ld,d); return 0;}
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT trans_getmultidiffcoeffs_(const integer* n, integer* ld, doublereal* d) {
        try { _ftrans(n)->getMultiDiffCoeffs(*ld,d); return 0;}
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT trans_setparameters_(const integer* n, integer* type, integer* k, doublereal* d) {
        try { _ftrans(n)->setParameters(*type, *k, d); return 0;}
        catch (CanteraError) { handleError(); return -1; }
    }

    //-------------------- Functions ---------------------------

//     status_t DLL_EXPORT import_phase_(const integer* nth, const integer* nxml, char* id, ftnlen lenid) {
//         thermo_t* thrm = th(nth);
//         XML_Node* node = _xml(nxml);
//         string idstr = f2string(id, lenid);
//         try {
//             importPhase(*node, thrm);
//             return 0;
//         }
//         catch (CanteraError) { handleError(); return -1; }
//     }

//     status_t DLL_EXPORT import_kinetics_(const integer* nxml, char* id, 
//         const integer* nphases, integer* ith, const integer* nkin, ftnlen lenid) {
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
//         catch (CanteraError) { handleError(); return -1; }
//     }


    status_t DLL_EXPORT ctphase_report_(const integer* nth, 
        char* buf, integer* show_thermo, ftnlen buflen) {
        try {
            bool stherm = (*show_thermo != 0);
            std::string s = report(*_fth(nth), stherm);
            if (int(s.size()) > buflen - 1) {
                return -(s.size() + 1);
            }
            copy(s.begin(), s.end(), buf);
            for (int nn = s.size(); nn < buflen; nn++) buf[nn] = ' ';
            return 0;
            
        }
        catch (CanteraError) { handleError(); return -1; }
    }

    status_t DLL_EXPORT ctgetcanteraerror_(char* buf, ftnlen buflen) {
        std::string e; // = "<no error>";
        //if (nErrors() > 0)
        e = lastErrorMessage();
        int n = min(e.size(), buflen-1);
        copy(e.begin(), e.begin() + n, buf);
        for (int nn = n; nn < buflen; nn++) buf[nn] = ' ';
        return 0;
    }

    status_t DLL_EXPORT ctaddcanteradirectory_(integer* buflen, char* buf) {
        addDirectory(std::string(buf));
        return 0;
    }


    status_t DLL_EXPORT ctbuildsolutionfromxml(char* src, integer* ixml, char* id, 
        integer* ith, integer* ikin, ftnlen lensrc, ftnlen lenid) {

        XML_Node* root = 0;
        if (*ixml > 0) root = _xml(ixml);

        thermo_t* t = _fth(ith);
        kinetics_t* k = _fkin(ikin);

        Kinetics& kin = *k;
        XML_Node *x, *r=0;
        if (root) r = &root->root();
	std::string srcS = f2string(src, lensrc);
	std::string idS  = f2string(id, lenid);
	if (srcS != "") {
           x = get_XML_Node(srcS, r);
	} else {
           x = get_XML_Node(idS, r);
	}
        // x = find_XML(f2string(src, lensrc), r, f2string(id,lenid), "", "phase");
        if (!x) return 0;
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
}
