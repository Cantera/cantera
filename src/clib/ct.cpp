/**
 *  @file ct.cpp
 *   Cantera interface library. This library of functions is designed
 *   to encapsulate Cantera functionality and make it available for
 *   use in languages and applications other than C++. A set of
 *   library functions is provided that are declared "extern C". All
 *   Cantera objects are stored and referenced by integers - no
 *   pointers are passed to or from the calling application.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ct.h"

// Cantera includes
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Interface.h"
#include "cantera/thermo/ThermoFactory.h"
#include "clib_utils.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/base/ExternalLogger.h"

using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Transport> TransportCabinet;
typedef Cabinet<Solution> SolutionCabinet;
typedef Cabinet<Reaction> ReactionCabinet;

template<> ThermoCabinet* ThermoCabinet::s_storage = 0;
template<> KineticsCabinet* KineticsCabinet::s_storage = 0;
template<> TransportCabinet* TransportCabinet::s_storage = 0;
template<> SolutionCabinet* SolutionCabinet::s_storage = 0;
template<> ReactionCabinet* ReactionCabinet::s_storage = 0;

/**
 * Exported functions.
 */
extern "C" {

    int ct_appdelete()
    {
        try {
            appdelete();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //--------------- Solution ------------------//

    int soln_newSolution(const char* infile,
                         const char* name,
                         const char* transport)
    {
        try {
            auto soln = newSolution(infile, name, transport);
            int id = SolutionCabinet::add(soln);
            // add associated objects
            ThermoCabinet::add(soln->thermo(), id);
            if (soln->kinetics()) {
                KineticsCabinet::add(soln->kinetics(), id);
            }
            if (soln->transport()) {
                TransportCabinet::add(soln->transport(), id);
            }
            return id;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int soln_newInterface(const char* infile,
                          const char* name,
                          int na,
                          const int* adjacent)
    {
        try {
            shared_ptr<Solution> soln;
            if (na) {
                vector<shared_ptr<Solution>> adj;
                for (int i = 0; i < na; i++) {
                    adj.push_back(SolutionCabinet::at(adjacent[i]));
                }
                soln = newInterface(infile, name, adj);
            } else {
                soln = newInterface(infile, name);
                // adjacent phases can be retrieved via soln_adjacent
            }
            // add associated objects
            int id = SolutionCabinet::add(soln);
            ThermoCabinet::add(soln->thermo(), id);
            if (soln->kinetics()) {
                KineticsCabinet::add(soln->kinetics(), id);
            }
            if (soln->transport()) {
                TransportCabinet::add(soln->transport(), id);
            }
            return id;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int soln_del(int n)
    {
        try {
            if (n >= 0 && n < SolutionCabinet::size()) {
                // remove all associated objects
                auto soln = SolutionCabinet::at(n);
                int index = ThermoCabinet::index(*(soln->thermo()), n);
                if (index >= 0) {
                    ThermoCabinet::del(index);
                }
                if (soln->kinetics()) {
                    index = KineticsCabinet::index(*(soln->kinetics()), n);
                    if (index >= 0) {
                        KineticsCabinet::del(index);
                    }
                }
                if (soln->transport()) {
                    index = TransportCabinet::index(*(soln->transport()), n);
                    if (index >= 0) {
                        TransportCabinet::del(index);
                    }
                }
            }
            SolutionCabinet::del(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int soln_name(int n, int buflen, char* buf)
    {
        try {
            string name = SolutionCabinet::at(n)->name();
            copyString(name, buf, buflen);
            return int(name.size());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int soln_thermo(int n)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            return ThermoCabinet::index(*soln->thermo(), n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int soln_kinetics(int n)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            if (!soln->kinetics()) {
                return -1;
            }
            return KineticsCabinet::index(*(soln->kinetics()), n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int soln_transport(int n)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            if (!soln->transport()) {
                return -1;
            }
            return TransportCabinet::index(*(soln->transport()), n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int soln_setTransportModel(int n, const char* model)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            TransportCabinet::del(
                TransportCabinet::index(*(soln->transport()), n));
            soln->setTransportModel(model);
            return TransportCabinet::add(soln->transport(), n);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t soln_nAdjacent(int n)
    {
        try {
            return SolutionCabinet::at(n)->nAdjacent();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int soln_adjacent(int n, int a)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            if (a < 0 || a >= static_cast<int>(soln->nAdjacent())) {
                return -1;
            }
            auto adj = soln->adjacent(a);
            // add associated objects
            int id = SolutionCabinet::add(adj);
            ThermoCabinet::add(adj->thermo(), id);
            if (adj->kinetics()) {
                KineticsCabinet::add(adj->kinetics(), id);
            }
            if (adj->transport()) {
                TransportCabinet::add(adj->transport(), id);
            }
            return id;
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int soln_adjacentName(int n, int a, int lennm, char* nm)
    {
        try {
            auto soln = SolutionCabinet::at(n);
            return static_cast<int>(copyString(soln->adjacentName(a), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //--------------- Phase ---------------------//

    size_t thermo_nElements(int n)
    {
        try {
            return ThermoCabinet::at(n)->nElements();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t thermo_nSpecies(int n)
    {
        try {
            return ThermoCabinet::at(n)->nSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double thermo_temperature(int n)
    {
        try {
            return ThermoCabinet::at(n)->temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setTemperature(int n, double t)
    {
        try {
            ThermoCabinet::at(n)->setTemperature(t);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    double thermo_density(int n)
    {
        try {
            return ThermoCabinet::at(n)->density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setDensity(int n, double rho)
    {
        try {
            ThermoCabinet::at(n)->setDensity(rho);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    double thermo_molarDensity(int n)
    {
        try {
            return ThermoCabinet::at(n)->molarDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_meanMolecularWeight(int n)
    {
        try {
            return ThermoCabinet::at(n)->meanMolecularWeight();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t thermo_elementIndex(int n, const char* nm)
    {
        try {
             size_t k = ThermoCabinet::at(n)->elementIndex(nm);
             if (k == npos) {
                throw CanteraError("thermo_elementIndex",
                                   "No such element {}.", nm);
             }
             return k;
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t thermo_speciesIndex(int n, const char* nm)
    {
        try {
             size_t k = ThermoCabinet::at(n)->speciesIndex(nm);
             if (k == npos) {
                throw CanteraError("thermo_speciesIndex",
                                   "No such species {}.", nm);
             }
             return k;
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int thermo_getMoleFractions(int n, size_t lenx, double* x)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(lenx);
            p->getMoleFractions(x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_moleFraction(int n, size_t k)
    {
        try {
            return ThermoCabinet::at(n)->moleFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_getMassFractions(int n, size_t leny, double* y)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(leny);
            p->getMassFractions(y);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_massFraction(int n, size_t k)
    {
        try {
            return ThermoCabinet::at(n)->massFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setMoleFractions(int n, size_t lenx, double* x, int norm)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(lenx);
            if (norm) {
                p->setMoleFractions(x);
            } else {
                p->setMoleFractions_NoNorm(x);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMoleFractionsByName(int n, const char* x)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->setMoleFractionsByName(x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMassFractions(int n, size_t leny,
                               double* y, int norm)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(leny);
            if (norm) {
                p->setMassFractions(y);
            } else {
                p->setMassFractions_NoNorm(y);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMassFractionsByName(int n, const char* y)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->setMassFractionsByName(y);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getAtomicWeights(int n, size_t lenm, double* atw)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkElementArraySize(lenm);
            const vector<double>& wt = p->atomicWeights();
            copy(wt.begin(), wt.end(), atw);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getMolecularWeights(int n, size_t lenm, double* mw)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(lenm);
            const vector<double>& wt = p->molecularWeights();
            copy(wt.begin(), wt.end(), mw);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getCharges(int n, size_t lenm, double* sc)
    {
        try {
            auto& p = ThermoCabinet::at(n);
            p->checkSpeciesArraySize(lenm);
            p->getCharges(sc);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getName(int n, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::at(n)->name(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setName(int n, const char* nm)
    {
        try {
            ThermoCabinet::at(n)->setName(nm);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getSpeciesName(int n, size_t k, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::at(n)->speciesName(k), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getElementName(int n, size_t m, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::at(n)->elementName(m), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_nAtoms(int n, size_t k, size_t m)
    {
        try {
            return ThermoCabinet::at(n)->nAtoms(k,m);
        } catch (...) {
            return handleAllExceptions(ERR, ERR);
        }
    }

    int thermo_addElement(int n, const char* name, double weight)
    {
        try {
            ThermoCabinet::at(n)->addElement(name, weight);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------- Thermo --------------------//

    int thermo_getEosType(int n, size_t leneos, char* eos)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::at(n)->type(), eos, leneos));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_enthalpy_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->enthalpy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_intEnergy_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->intEnergy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_entropy_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->entropy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_gibbs_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->gibbs_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cp_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->cp_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cv_mole(int n)
    {
        try {
            return ThermoCabinet::at(n)->cv_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_pressure(int n)
    {
        try {
            return ThermoCabinet::at(n)->pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_enthalpy_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_intEnergy_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_entropy_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->entropy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_gibbs_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->gibbs_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cp_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->cp_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cv_mass(int n)
    {
        try {
            return ThermoCabinet::at(n)->cv_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_electricPotential(int n)
    {
        try {
            return ThermoCabinet::at(n)->electricPotential();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_chemPotentials(int n, size_t lenm, double* murt)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getChemPotentials(murt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_electrochemPotentials(int n, size_t lenm, double* emu)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getElectrochemPotentials(emu);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setPressure(int n, double p)
    {
        try {
            if (p < 0.0) throw CanteraError("thermo_setPressure",
                                                "pressure cannot be negative");
            ThermoCabinet::at(n)->setPressure(p);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TP(int n, double* vals)
    {
        try{
            ThermoCabinet::at(n)->setState_TP(vals[0], vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TD(int n, double* vals)
    {
        try{
            ThermoCabinet::at(n)->setState_TD(vals[0], vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_DP(int n, double* vals)
    {
        try{
            ThermoCabinet::at(n)->setState_DP(vals[0], vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_HP(int n, double* vals)
    {
        try {
            if (vals[1] < 0.0) {
                throw CanteraError("thermo_set_HP",
                                   "pressure cannot be negative");
            }
            ThermoCabinet::at(n)->setState_HP(vals[0],vals[1]);
            if (ThermoCabinet::at(n)->temperature() < 0.0) {
                throw CanteraError("thermo_set_HP",
                                   "temperature cannot be negative");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_UV(int n, double* vals)
    {
        try {
            if (vals[1] < 0.0) {
                throw CanteraError("thermo_set_UV",
                                   "specific volume cannot be negative");
            }
            ThermoCabinet::at(n)->setState_UV(vals[0],vals[1]);
            if (ThermoCabinet::at(n)->temperature() < 0.0) {
                throw CanteraError("thermo_set_UV",
                                   "temperature cannot be negative");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_SV(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_SV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_SP(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_SP(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_ST(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_ST(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TV(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_TV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_PV(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_PV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_UP(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_UP(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_VH(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_VH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TH(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_TH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_SH(int n, double* vals)
    {
        try {
            ThermoCabinet::at(n)->setState_SH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_equilibrate(int n, const char* XY, const char* solver,
                           double rtol, int maxsteps, int maxiter, int loglevel)
    {
        try {
            ThermoCabinet::at(n)->equilibrate(XY, solver, rtol, maxsteps,
                                              maxiter, 0, loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_refPressure(int n)
    {
        try {
            return ThermoCabinet::at(n)->refPressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_minTemp(int n, int k)
    {
        try {
            auto& ph = ThermoCabinet::at(n);
            if (k != -1) {
                ph->checkSpeciesIndex(k);
                return ph->minTemp(k);
            } else {
                return ph->minTemp();
            }
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_maxTemp(int n, int k)
    {
        try {
            auto& ph = ThermoCabinet::at(n);
            if (k != -1) {
                ph->checkSpeciesIndex(k);
                return ph->maxTemp(k);
            } else {
                return ph->maxTemp();
            }
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }


    int thermo_getEnthalpies_RT(int n, size_t lenm, double* h_rt)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getEnthalpy_RT_ref(h_rt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getEntropies_R(int n, size_t lenm, double* s_r)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getEntropy_R_ref(s_r);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getCp_R(int n, size_t lenm, double* cp_r)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getCp_R_ref(cp_r);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setElectricPotential(int n, double v)
    {
        try {
            ThermoCabinet::at(n)->setElectricPotential(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getPartialMolarEnthalpies(int n, size_t lenm, double* pmh)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getPartialMolarEnthalpies(pmh);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getPartialMolarEntropies(int n, size_t lenm, double* pms)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getPartialMolarEntropies(pms);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getPartialMolarIntEnergies(int n, size_t lenm, double* pmu)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getPartialMolarIntEnergies(pmu);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getPartialMolarCp(int n, size_t lenm, double* pmcp)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getPartialMolarCp(pmcp);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getPartialMolarVolumes(int n, size_t lenm, double* pmv)
    {
        try {
            auto& thrm = ThermoCabinet::at(n);
            thrm->checkSpeciesArraySize(lenm);
            thrm->getPartialMolarVolumes(pmv);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_thermalExpansionCoeff(int n)
    {
        try {
            return ThermoCabinet::at(n)->thermalExpansionCoeff();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_isothermalCompressibility(int n)
    {
        try {
            return ThermoCabinet::at(n)->isothermalCompressibility();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    //-------------- pure fluids ---------------//

    double thermo_critTemperature(int n)
    {
        try {
            return ThermoCabinet::at(n)->critTemperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_critPressure(int n)
    {
        try {
            return ThermoCabinet::at(n)->critPressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_critDensity(int n)
    {
        try {
            return ThermoCabinet::at(n)->critDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_vaporFraction(int n)
    {
        try {
            return ThermoCabinet::as<PureFluidPhase>(n)->vaporFraction();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_satTemperature(int n, double p)
    {
        try {
            return ThermoCabinet::at(n)->satTemperature(p);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_satPressure(int n, double t)
    {
        try {
            return ThermoCabinet::at(n)->satPressure(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setState_Psat(int n, double p, double x)
    {
        try {
            ThermoCabinet::as<PureFluidPhase>(n)->setState_Psat(p, x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setState_Tsat(int n, double t, double x)
    {
        try {
            ThermoCabinet::as<PureFluidPhase>(n)->setState_Tsat(t, x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------- Kinetics ------------------//

    int kin_getType(int n, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(KineticsCabinet::at(n)->kineticsType(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t kin_start(int n, int p)
    {
        try {
            return KineticsCabinet::at(n)->kineticsSpeciesIndex(0,p);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_speciesIndex(int n, const char* nm)
    {
        try {
            return KineticsCabinet::at(n)->kineticsSpeciesIndex(nm);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    //---------------------------------------

    size_t kin_nSpecies(int n)
    {
        try {
            return KineticsCabinet::at(n)->nTotalSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_nReactions(int n)
    {
        try {
            return KineticsCabinet::at(n)->nReactions();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_nPhases(int n)
    {
        try {
            return KineticsCabinet::at(n)->nPhases();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_phaseIndex(int n, const char* ph)
    {
        try {
            size_t k = KineticsCabinet::at(n)->phaseIndex(ph);
            if (k == npos) {
                throw CanteraError("kin_phaseIndex",
                                   "No such phase {}.", ph);
            }
            return k;
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double kin_reactantStoichCoeff(int n, int k, int i)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkSpeciesIndex(k);
            kin->checkReactionIndex(i);
            return kin->reactantStoichCoeff(k,i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double kin_productStoichCoeff(int n, int k, int i)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkSpeciesIndex(k);
            kin->checkReactionIndex(i);
            return kin->productStoichCoeff(k,i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int kin_getReactionType(int n, int i, size_t len, char* buf)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkReactionIndex(i);
            return static_cast<int>(copyString(kin->reaction(i)->type(), buf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getFwdRatesOfProgress(int n, size_t len, double* fwdROP)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getFwdRatesOfProgress(fwdROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getRevRatesOfProgress(int n, size_t len, double* revROP)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getRevRatesOfProgress(revROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_isReversible(int n, int i)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkReactionIndex(i);
            return (int) kin->isReversible(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getNetRatesOfProgress(int n, size_t len, double* netROP)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getNetRatesOfProgress(netROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getFwdRateConstants(int n, size_t len, double* kfwd)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getFwdRateConstants(kfwd);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getRevRateConstants(int n, int doIrreversible, size_t len, double* krev)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getRevRateConstants(krev, doIrreversible != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getDelta(int n, int job, size_t len, double* delta)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            switch (job) {
            case 0:
                k->getDeltaEnthalpy(delta);
                break;
            case 1:
                k->getDeltaGibbs(delta);
                break;
            case 2:
                k->getDeltaEntropy(delta);
                break;
            case 3:
                k->getDeltaSSEnthalpy(delta);
                break;
            case 4:
                k->getDeltaSSGibbs(delta);
                break;
            case 5:
                k->getDeltaSSEntropy(delta);
                break;
            default:
                return ERR;
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getCreationRates(int n, size_t len, double* cdot)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkSpeciesArraySize(len);
            k->getCreationRates(cdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getDestructionRates(int n, size_t len, double* ddot)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkSpeciesArraySize(len);
            k->getDestructionRates(ddot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getNetProductionRates(int n, size_t len, double* wdot)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkSpeciesArraySize(len);
            k->getNetProductionRates(wdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getSourceTerms(int n, size_t len, double* ydot)
    {
        try {
            // @todo This function only works for single phase kinetics
            auto& k = KineticsCabinet::at(n);
            ThermoPhase& p = k->thermo();
            size_t nsp = p.nSpecies();
            k->checkSpeciesArraySize(len);
            k->checkSpeciesArraySize(nsp);
            k->getNetProductionRates(ydot);
            double rho_inv = 1.0 / p.density();
            for (size_t k = 0; k < nsp; k++) {
                ydot[k] *= p.molecularWeight(k) * rho_inv;
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double kin_multiplier(int n, int i)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkReactionIndex(i);
            return kin->multiplier(i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t kin_phase(int n, size_t i)
    {
        try {
            auto& kin = KineticsCabinet::at(n);
            kin->checkPhaseIndex(i);
            return ThermoCabinet::index(kin->thermo(i));
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int kin_getEquilibriumConstants(int n, size_t len, double* kc)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionArraySize(len);
            k->getEquilibriumConstants(kc);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getReactionString(int n, int i, int len, char* buf)
    {
        try {
            auto& k = KineticsCabinet::at(n);
            k->checkReactionIndex(i);
            return static_cast<int>(copyString(k->reaction(i)->equation(), buf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_setMultiplier(int n, int i, double v)
    {
        try {
            if (v >= 0.0) {
                auto& kin = KineticsCabinet::at(n);
                kin->checkReactionIndex(i);
                kin->setMultiplier(i,v);
                return 0;
            } else {
                return ERR;
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_advanceCoverages(int n, double tstep)
    {
        try {
            KineticsCabinet::as<InterfaceKinetics>(n)->advanceCoverages(tstep);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //------------------- Transport ---------------------------

    int trans_transportModel(int i, int lennm, char* nm)
    {
        try {
            return static_cast<int>(
                copyString(TransportCabinet::at(i)->transportModel(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_viscosity(int n)
    {
        try {
            return TransportCabinet::at(n)->viscosity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_electricalConductivity(int n)
    {
        try {
            return TransportCabinet::at(n)->electricalConductivity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_thermalConductivity(int n)
    {
        try {
            return TransportCabinet::at(n)->thermalConductivity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getThermalDiffCoeffs(int n, int ldt, double* dt)
    {
        try {
            auto& tr = TransportCabinet::at(n);
            tr->checkSpeciesArraySize(ldt);
            tr->getThermalDiffCoeffs(dt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMixDiffCoeffs(int n, int ld, double* d)
    {
        try {
            auto& tr = TransportCabinet::at(n);
            tr->checkSpeciesArraySize(ld);
            tr->getMixDiffCoeffs(d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getBinDiffCoeffs(int n, int ld, double* d)
    {
        try {
            // @todo length of d should be passed for bounds checking
            auto& tr = TransportCabinet::at(n);
            tr->checkSpeciesArraySize(ld);
            tr->getBinaryDiffCoeffs(ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMultiDiffCoeffs(int n, int ld, double* d)
    {
        try {
            // @todo length of d should be passed for bounds checking
            auto& tr = TransportCabinet::at(n);
            tr->checkSpeciesArraySize(ld);
            tr->getMultiDiffCoeffs(ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMolarFluxes(int n, const double* state1,
                             const double* state2, double delta, double* fluxes)
    {
        try {
            TransportCabinet::at(n)->getMolarFluxes(state1, state2, delta, fluxes);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMassFluxes(int n, const double* state1,
                            const double* state2, double delta, double* fluxes)
    {
        try {
            TransportCabinet::at(n)->getMassFluxes(state1, state2, delta, fluxes);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------------- Functions ---------------------------

    int thermo_report(int nth, int show_thermo, double threshold, int ibuf, char* buf)
    {
        try {
            bool stherm = (show_thermo != 0);
            return static_cast<int>(
                copyString(ThermoCabinet::at(nth)->report(stherm, threshold),
                buf, ibuf));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_print(int nth, int show_thermo, double threshold)
    {
        try {
            bool stherm = (show_thermo != 0);
            writelog(ThermoCabinet::at(nth)->report(stherm, threshold));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_getCanteraError(int buflen, char* buf)
    {
        try {
            string e = Application::Instance()->lastErrorMessage();
            copyString(e, buf, buflen);
            return int(e.size());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_addCanteraDirectory(size_t buflen, const char* buf)
    {
        try {
            addDirectory(buf);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_getDataDirectories(const char* sep, int buflen, char* buf)
    {
        try {
            return static_cast<int>(copyString(getDataDirectories(sep), buf, buflen));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_getCanteraVersion(int buflen, char* buf)
    {
        try {
            return static_cast<int>(copyString(CANTERA_VERSION, buf, buflen));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_getGitCommit(int buflen, char* buf)
    {
        try {
            return static_cast<int>(copyString(gitCommit(), buf, buflen));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_suppress_thermo_warnings(int suppress)
    {
        try {
            suppress_thermo_warnings(suppress != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_use_legacy_rate_constants(int legacy)
    {
        try {
            use_legacy_rate_constants(legacy != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_setLogWriter(void* logger)
    {
        try {
            Logger* logwriter = (Logger*)logger;
            setLogger(logwriter);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_setLogCallback(LogCallback writer)
    {
        static unique_ptr<Logger> logwriter;
        try {
            logwriter = make_unique<ExternalLogger>(writer);
            setLogger(logwriter.get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_resetStorage()
    {
        try {
            SolutionCabinet::reset();
            ThermoCabinet::reset();
            KineticsCabinet::reset();
            TransportCabinet::reset();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearStorage()
    {
        try {
            SolutionCabinet::clear();
            ThermoCabinet::clear();
            KineticsCabinet::clear();
            TransportCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_parent(int n)
    {
        try {
            return ThermoCabinet::parent(n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int thermo_size()
    {
        try {
            return ThermoCabinet::size();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_parent(int n)
    {
        try {
            return KineticsCabinet::parent(n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int trans_parent(int n)
    {
        try {
            return TransportCabinet::parent(n);
        } catch (...) {
            return handleAllExceptions(-2, ERR);
        }
    }

    int thermo_del(int n)
    {
        // function is needed for sourcegen
        // no-op: object is managed by Solution
        return 0;
    }

    int kin_del(int n)
    {
        // function is needed for sourcegen
        // no-op: object is managed by Solution
        return 0;
    }

    int trans_del(int n)
    {
        // function is needed for sourcegen
        // no-op: object is managed by Solution
        return 0;
    }
}
