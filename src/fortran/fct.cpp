/**
 *   Cantera Fortran interface library. This library of functions is designed
 *   to encapsulate Cantera functionality and make it available for
 *   use in languages and applications other than C++. A set of
 *   library functions is provided that are declared "extern C". All
 *   Cantera objects are stored and referenced by integers - no
 *   pointers are passed to or from the calling application.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// Cantera includes
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/ctml.h"
#include "cantera/kinetics/importKinetics.h"
#include "clib/Cabinet.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include "cantera/clib/clib_defs.h"
#include "cantera/clib/ct.h"

using namespace Cantera;

typedef Cabinet<XML_Node, false> XmlCabinet;
typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Transport> TransportCabinet;

typedef integer status_t;

namespace {

XML_Node* _xml(const integer* n)
{
    return &XmlCabinet::item(*n);
}

ThermoPhase* _fph(const integer* n)
{
    return &ThermoCabinet::item(*n);
}

static Kinetics* _fkin(const integer* n)
{
    return &KineticsCabinet::item(*n);
}

ThermoPhase* _fth(const integer* n)
{
    return &ThermoCabinet::item(*n);
}

Transport* _ftrans(const integer* n)
{
    return &TransportCabinet::item(*n);
}

} // unnamed namespace

std::string f2string(const char* s, ftnlen n)
{
    int k;
    std::string ss = "";
    for (k = 0; k < n; k++) {
        if (s[k] == '\0') {
            break;
        }
        ss += s[k];
    }
    return ss;
}

/**
 * Exported functions.
 */
extern "C" {

    status_t cantera_error_(const char* proc, const char* msg,
                            ftnlen proclen, ftnlen msglen)
    {
        std::string sproc = f2string(proc, proclen);
        std::string smsg = f2string(msg, msglen);
        throw CanteraError(sproc, smsg);
        return -1;
    }

    //--------------- Phase ---------------------//

    status_t phase_getname_(const integer* n, char* nm,
                            ftnlen lennm)
    {
        try {
            std::string pnm = _fph(n)->name();
            int lout = std::min(lennm, (int) pnm.size());
            std::copy(pnm.c_str(), pnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    integer phase_nelements_(const integer* n)
    {
        try {
            return _fph(n)->nElements();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer phase_nspecies_(const integer* n)
    {
        try {
            return _fph(n)->nSpecies();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal phase_temperature_(const integer* n)
    {
        try {
            return _fph(n)->temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_settemperature_(const integer* n, doublereal* t)
    {
        try {
            _fph(n)->setTemperature(*t);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal phase_density_(const integer* n)
    {
        try {
            return _fph(n)->density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_setdensity_(const integer* n, doublereal* rho)
    {
        try {
            _fph(n)->setDensity(*rho);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal phase_molardensity_(const integer* n)
    {
        try {
            return _fph(n)->molarDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal phase_meanmolecularweight_(const integer* n)
    {
        try {
            return _fph(n)->meanMolecularWeight();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    integer phase_elementindex_(const integer* n, char* nm, ftnlen lennm)
    {
        try {
            std::string elnm = f2string(nm, lennm);
            return _fph(n)->elementIndex(elnm) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer phase_speciesindex_(const integer* n, char* nm, ftnlen lennm)
    {
        try {
            std::string spnm = f2string(nm, lennm);
            return _fph(n)->speciesIndex(spnm) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t phase_getmolefractions_(const integer* n, doublereal* x)
    {
        try {
            _fph(n)->getMoleFractions(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal phase_molefraction_(const integer* n, integer* k)
    {
        try {
            return _fph(n)->moleFraction(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }

    }

    status_t phase_getmassfractions_(const integer* n, doublereal* y)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->getMassFractions(y);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal phase_massfraction_(const integer* n, integer* k)
    {
        try {
            return _fph(n)->massFraction(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_setmolefractions_(const integer* n, double* x, const integer* norm)
    {
        try {
            ThermoPhase* p = _fph(n);
            if (*norm) {
                p->setMoleFractions(x);
            } else {
                p->setMoleFractions_NoNorm(x);
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmolefractionsbyname_(const integer* n, char* x, ftnlen lx)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->setMoleFractionsByName(f2string(x, lx));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmassfractions_(const integer* n, doublereal* y, const integer* norm)
    {
        try {
            ThermoPhase* p = _fph(n);
            if (*norm) {
                p->setMassFractions(y);
            } else {
                p->setMassFractions_NoNorm(y);
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmassfractionsbyname_(const integer* n, char* y, ftnlen leny)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->setMassFractionsByName(f2string(y, leny));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getatomicweights_(const integer* n, doublereal* atw)
    {
        try {
            ThermoPhase* p = _fph(n);
            const vector_fp& wt = p->atomicWeights();
            copy(wt.begin(), wt.end(), atw);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getmolecularweights_(const integer* n, doublereal* mw)
    {
        try {
            ThermoPhase* p = _fph(n);
            const vector_fp& wt = p->molecularWeights();
            copy(wt.begin(), wt.end(), mw);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getspeciesname_(const integer* n, integer* k, char* nm, ftnlen lennm)
    {
        try {
            std::string spnm = _fph(n)->speciesName(*k-1);
            int lout = std::min(lennm, (int) spnm.size());
            std::copy(spnm.c_str(), spnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getelementname_(const integer* n, integer* m, char* nm, ftnlen lennm)
    {
        try {
            std::string elnm = _fph(n)->elementName(*m-1);
            int lout = std::min(lennm, (int) elnm.size());
            std::copy(elnm.c_str(), elnm.c_str() + lout, nm);
            for (int nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal phase_natoms_(const integer* n, integer* k, integer* m)
    {
        try {
            return _fph(n)->nAtoms(*k-1,*m-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    //-------------- Thermo --------------------//

    integer th_newfromfile_(char* filename, char* phasename, ftnlen lenf, ftnlen lenp)
    {
        try {
            thermo_t* th = newPhase(f2string(filename, lenf),
                                    f2string(phasename, lenp));
            return ThermoCabinet::add(th);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer newthermofromxml_(integer* mxml)
    {
        try {
            XML_Node* x = _xml(mxml);
            thermo_t* th = newPhase(*x);
            return ThermoCabinet::add(th);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer th_geteostype_(const integer* n, char* buf, ftnlen lenbuf)
    {
        try {
            return copyString(_fth(n)->type(), buf, lenbuf);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal th_enthalpy_mole_(const integer* n)
    {
        try {
            return _fth(n)->enthalpy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_intenergy_mole_(const integer* n)
    {
        try {
            return _fth(n)->intEnergy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_entropy_mole_(const integer* n)
    {
        try {
            return _fth(n)->entropy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_gibbs_mole_(const integer* n)
    {
        try {
            return _fth(n)->gibbs_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_cp_mole_(const integer* n)
    {
        try {
            return _fth(n)->cp_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_cv_mole_(const integer* n)
    {
        try {
            return _fth(n)->cv_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_pressure_(const integer* n)
    {
        try {
            return _fth(n)->pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_enthalpy_mass_(const integer* n)
    {
        try {
            return _fth(n)->enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_intenergy_mass_(const integer* n)
    {
        try {
            return _fth(n)->intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_entropy_mass_(const integer* n)
    {
        try {
            return _fth(n)->entropy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_gibbs_mass_(const integer* n)
    {
        try {
            return _fth(n)->gibbs_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_cp_mass_(const integer* n)
    {
        try {
            return _fth(n)->cp_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_cv_mass_(const integer* n)
    {
        try {
            return _fth(n)->cv_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t th_chempotentials_(const integer* n, doublereal* murt)
    {
        try {
            thermo_t* thrm = _fth(n);
            thrm->getChemPotentials(murt);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_setpressure_(const integer* n, doublereal* p)
    {
        try {
            _fth(n)->setPressure(*p);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_hp_(const integer* n, doublereal* v1, doublereal* v2)
    {
        try {
            _fth(n)->setState_HP(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_uv_(const integer* n, doublereal* v1, doublereal* v2)
    {
        try {
            _fth(n)->setState_UV(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_sv_(const integer* n, doublereal* v1, doublereal* v2)
    {
        try {
            _fth(n)->setState_SV(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_sp_(const integer* n, doublereal* v1, doublereal* v2)
    {
        try {
            _fth(n)->setState_SP(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_equil_(const integer* n, char* XY, char* solver, doublereal* rtol, int* max_steps, int* max_iter, int* estimate_equil, int* log_level, ftnlen lenxy, ftnlen lensolver)
    {
        try {
            _fth(n)->equilibrate(f2string(XY,lenxy), f2string(solver,lensolver), *rtol, *max_steps, *max_iter, *estimate_equil, *log_level);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal th_refpressure_(const integer* n)
    {
        try {
            return _fth(n)->refPressure();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal th_mintemp_(const integer* n, integer* k)
    {
        try {
            return _fth(n)->minTemp(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal th_maxtemp_(const integer* n, integer* k)
    {
        try {
            return _fth(n)->maxTemp(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t th_getenthalpies_rt_(const integer* n, doublereal* h_rt)
    {
        try {
            thermo_t* thrm = _fth(n);
            thrm->getEnthalpy_RT(h_rt);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getentropies_r_(const integer* n, doublereal* s_r)
    {
        try {
            thermo_t* thrm = _fth(n);
            thrm->getEntropy_R(s_r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getcp_r_(const integer* n, integer* lenm, doublereal* cp_r)
    {
        try {
            thermo_t* thrm = _fth(n);
            thrm->getCp_R(cp_r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getpartialmolarintenergies_r_(const integer* n, doublereal* ie)
    {
        try {
            thermo_t* thrm = _fth(n);
            thrm->getPartialMolarIntEnergies(ie);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    //-------------- Kinetics ------------------//

    integer kin_newfromfile_(const char* filename, const char* phasename,
                             integer* reactingPhase, const integer* neighbor1,
                             const integer* neighbor2, const integer* neighbor3,
                             const integer* neighbor4, ftnlen nlen, ftnlen plen)
    {
        try {
            std::vector<thermo_t*> phases;
            phases.push_back(_fth(reactingPhase));
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
            auto kin = newKinetics(phases, f2string(filename, nlen),
                                   f2string(phasename, plen));
            return KineticsCabinet::add(kin.release());
        } catch (...) {
            return handleAllExceptions(999, ERR);
        }
    }

    integer newkineticsfromxml_(integer* mxml, integer* iphase,
                                const integer* neighbor1, const integer* neighbor2, const integer* neighbor3,
                                const integer* neighbor4)
    {
        try {
            XML_Node* x = _xml(mxml);
            std::vector<thermo_t*> phases;
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
                return KineticsCabinet::add(kin);
            } else {
                return 0;
            }
        } catch (...) {
            return handleAllExceptions(999, ERR);
        }
    }

    //-------------------------------------
    integer kin_gettype_(const integer* n, char* buf, ftnlen buflen)
    {
        try {
            return copyString(_fkin(n)->kineticsType(), buf, buflen);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_start_(const integer* n, integer* p)
    {
        try {
            return _fkin(n)->kineticsSpeciesIndex(0, *p)+1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_speciesindex_(const integer* n, const char* nm, const char* ph,
                              ftnlen lennm, ftnlen lenph)
    {
        try {
            return _fkin(n)->kineticsSpeciesIndex(f2string(nm, lennm),
                                                  f2string(ph, lenph)) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //---------------------------------------

    integer kin_ntotalspecies_(const integer* n)
    {
        try {
            return _fkin(n)->nTotalSpecies();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_nreactions_(const integer* n)
    {
        try {
            return _fkin(n)->nReactions();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_nphases_(const integer* n)
    {
        try {
            return _fkin(n)->nPhases();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_phaseindex_(const integer* n, const char* ph, ftnlen lenph)
    {
        try {
            return _fkin(n)->phaseIndex(f2string(ph, lenph));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal kin_reactantstoichcoeff_(const integer* n, integer* k, integer* i)
    {
        try {
            return _fkin(n)->reactantStoichCoeff(*k-1,*i-1);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal kin_productstoichcoeff_(const integer* n, integer* k, integer* i)
    {
        try {
            return _fkin(n)->productStoichCoeff(*k-1,*i-1);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_reactiontype_(const integer* n, integer* i)
    {
        try {
            return _fkin(n)->reactionType(*i-1);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t kin_getfwdratesofprogress_(const integer* n, doublereal* fwdROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getFwdRatesOfProgress(fwdROP);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getrevratesofprogress_(const integer* n, doublereal* revROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getRevRatesOfProgress(revROP);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    integer kin_isreversible_(const integer* n, integer* i)
    {
        try {
            return (int)_fkin(n)->isReversible(*i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t kin_getnetratesofprogress_(const integer* n, doublereal* netROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getNetRatesOfProgress(netROP);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getcreationrates_(const integer* n, doublereal* cdot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getCreationRates(cdot);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getdestructionrates_(const integer* n, doublereal* ddot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getDestructionRates(ddot);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getnetproductionrates_(const integer* n, doublereal* wdot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getNetProductionRates(wdot);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal kin_multiplier_(const integer* n, integer* i)
    {
        try {
            return _fkin(n)->multiplier(*i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t kin_getequilibriumconstants_(const integer* n, doublereal* kc)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getEquilibriumConstants(kc);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getreactionstring_(const integer* n, integer* i, char* buf, ftnlen lenbuf)
    {
        try {
            Kinetics* k = _fkin(n);
            std::string r = k->reactionString(*i-1);
            int lout = std::min(lenbuf, (int) r.size());
            std::copy(r.c_str(), r.c_str() + lout, buf);
            for (int nn = lout; nn < lenbuf; nn++) {
                buf[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_setmultiplier_(const integer* n, integer* i, doublereal* v)
    {
        try {
            _fkin(n)->setMultiplier(*i-1,*v);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_advancecoverages_(const integer* n, doublereal* tstep)
    {
        try {
            Kinetics* k = _fkin(n);
            if (k->kineticsType() == "Surf" || k->kineticsType() == "Edge") {
                ((InterfaceKinetics*)k)->advanceCoverages(*tstep);
            } else {
                throw CanteraError("kin_advanceCoverages",
                                   "wrong kinetics manager type");
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    //------------------- Transport ---------------------------

    integer newtransport_(char* model, integer* ith,
                          integer* loglevel, ftnlen lenmodel)
    {
        try {
            std::string mstr = f2string(model, lenmodel);
            thermo_t* t = _fth(ith);
            Transport* tr = newTransportMgr(mstr, t, *loglevel);
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer trans_newdefault_(integer* ith, integer* loglevel, ftnlen lenmodel)
    {
        try {
            thermo_t* t = _fth(ith);
            Transport* tr = newDefaultTransportMgr(t, *loglevel);
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal trans_viscosity_(const integer* n)
    {
        try {
            return _ftrans(n)->viscosity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal trans_electricalconductivity_(const integer* n)
    {
        try {
            return _ftrans(n)->electricalConductivity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal trans_thermalconductivity_(const integer* n)
    {
        try {
            return _ftrans(n)->thermalConductivity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t trans_getthermaldiffcoeffs_(const integer* n, doublereal* dt)
    {
        try {
            _ftrans(n)->getThermalDiffCoeffs(dt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffs_(const integer* n, doublereal* d)
    {
        try {
            _ftrans(n)->getMixDiffCoeffs(d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffsmass_(const integer* n, doublereal* d)
    {
        try {
            _ftrans(n)->getMixDiffCoeffsMass(d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffsmole_(const integer* n, doublereal* d)
    {
        try {
            _ftrans(n)->getMixDiffCoeffsMole(d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getbindiffcoeffs_(const integer* n, integer* ld, doublereal* d)
    {
        try {
            _ftrans(n)->getBinaryDiffCoeffs(*ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmultidiffcoeffs_(const integer* n, integer* ld, doublereal* d)
    {
        try {
            _ftrans(n)->getMultiDiffCoeffs(*ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_setparameters_(const integer* n, integer* type, integer* k, doublereal* d)
    {
        try {
            _ftrans(n)->setParameters(*type, *k, d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------------- Functions ---------------------------

    status_t ctphase_report_(const integer* nth,
                             char* buf, integer* show_thermo, ftnlen buflen)
    {
        try {
            bool stherm = (*show_thermo != 0);
            std::string s = _fth(nth)->report(stherm);
            if (int(s.size()) > buflen - 1) {
                return -(s.size() + 1);
            }
            copy(s.begin(), s.end(), buf);
            for (int nn = s.size(); nn < buflen; nn++) {
                buf[nn] = ' ';
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t ctgetcanteraerror_(char* buf, ftnlen buflen)
    {
        try {
            std::string e;
            e = Application::Instance()->lastErrorMessage();
            int n = std::min((int) e.size(), buflen-1);
            copy(e.begin(), e.begin() + n, buf);
            for (int nn = n; nn < buflen; nn++) {
                buf[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t ctaddcanteradirectory_(integer* buflen, char* buf)
    {
        try {
            addDirectory(std::string(buf));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t ctbuildsolutionfromxml(char* src, integer* ixml, char* id,
                                    integer* ith, integer* ikin, ftnlen lensrc, ftnlen lenid)
    {
        try {
            XML_Node* root = 0;
            if (*ixml > 0) {
                root = _xml(ixml);
            }

            thermo_t* t = _fth(ith);
            Kinetics* k = _fkin(ikin);

            XML_Node* x, *r=0;
            if (root) {
                r = &root->root();
            }
            std::string srcS = f2string(src, lensrc);
            std::string idS = f2string(id, lenid);
            if (srcS != "") {
                x = get_XML_Node(srcS, r);
            } else {
                x = get_XML_Node(idS, r);
            }
            if (!x) {
                return 0;
            }
            importPhase(*x, t);
            k->addPhase(*t);
            k->init();
            installReactionArrays(*x, *k, x->id());
            t->setState_TP(300.0, OneAtm);
            if (r) {
                if (&x->root() != &r->root()) {
                    delete &x->root();
                }
            } else {
                delete &x->root();
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }
}
