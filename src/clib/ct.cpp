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

#define CANTERA_USE_INTERNAL
#include "cantera/clib/ct.h"

// Cantera includes
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/ctml.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/ThermoFactory.h"
#include "Cabinet.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/PureFluidPhase.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Transport> TransportCabinet;
typedef Cabinet<XML_Node, false> XmlCabinet;

template<> ThermoCabinet* ThermoCabinet::s_storage = 0;
template<> KineticsCabinet* KineticsCabinet::s_storage = 0;
template<> TransportCabinet* TransportCabinet::s_storage = 0;

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

    //--------------- Phase ---------------------//

    size_t thermo_nElements(int n)
    {
        try {
            return ThermoCabinet::item(n).nElements();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t thermo_nSpecies(int n)
    {
        try {
            return ThermoCabinet::item(n).nSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    doublereal thermo_temperature(int n)
    {
        try {
            return ThermoCabinet::item(n).temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setTemperature(int n, double t)
    {
        try {
            ThermoCabinet::item(n).setTemperature(t);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal thermo_density(int n)
    {
        try {
            return ThermoCabinet::item(n).density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setDensity(int n, double rho)
    {
        if (rho < 0.0) {
            return -1;
        }
        try {
            ThermoCabinet::item(n).setDensity(rho);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal thermo_molarDensity(int n)
    {
        try {
            return ThermoCabinet::item(n).molarDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setMolarDensity(int n, double ndens)
    {
        if (ndens < 0.0) {
            return -1;
        }
        try {
            ThermoCabinet::item(n).setMolarDensity(ndens);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    doublereal thermo_meanMolecularWeight(int n)
    {
        try {
            return ThermoCabinet::item(n).meanMolecularWeight();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t thermo_elementIndex(int n, const char* nm)
    {
        try {
            return ThermoCabinet::item(n).elementIndex(nm);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t thermo_speciesIndex(int n, const char* nm)
    {
        try {
            return ThermoCabinet::item(n).speciesIndex(nm);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int thermo_getMoleFractions(int n, size_t lenx, double* x)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(lenx);
            p.getMoleFractions(x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal thermo_moleFraction(int n, size_t k)
    {
        try {
            return ThermoCabinet::item(n).moleFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_getMassFractions(int n, size_t leny, double* y)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(leny);
            p.getMassFractions(y);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal thermo_massFraction(int n, size_t k)
    {
        try {
            return ThermoCabinet::item(n).massFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setMoleFractions(int n, size_t lenx, double* x, int norm)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(lenx);
            if (norm) {
                p.setMoleFractions(x);
            } else {
                p.setMoleFractions_NoNorm(x);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMoleFractionsByName(int n, const char* x)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.setMoleFractionsByName(x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMassFractions(int n, size_t leny,
                               double* y, int norm)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(leny);
            if (norm) {
                p.setMassFractions(y);
            } else {
                p.setMassFractions_NoNorm(y);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setMassFractionsByName(int n, const char* y)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.setMassFractionsByName(y);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getAtomicWeights(int n, size_t lenm, double* atw)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkElementArraySize(lenm);
            const vector_fp& wt = p.atomicWeights();
            copy(wt.begin(), wt.end(), atw);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getMolecularWeights(int n, size_t lenm, double* mw)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(lenm);
            const vector_fp& wt = p.molecularWeights();
            copy(wt.begin(), wt.end(), mw);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getCharges(int n, size_t lenm, double* sc)
    {
        try {
            ThermoPhase& p = ThermoCabinet::item(n);
            p.checkSpeciesArraySize(lenm);
            p.getCharges(sc);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getName(int n, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::item(n).name(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setName(int n, const char* nm)
    {
        try {
            ThermoCabinet::item(n).setName(nm);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getSpeciesName(int n, size_t k, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::item(n).speciesName(k), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getElementName(int n, size_t m, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::item(n).elementName(m), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    doublereal thermo_nAtoms(int n, size_t k, size_t m)
    {
        try {
            return ThermoCabinet::item(n).nAtoms(k,m);
        } catch (...) {
            return handleAllExceptions(ERR, ERR);
        }
    }

    int thermo_addElement(int n, const char* name, doublereal weight)
    {
        try {
            ThermoCabinet::item(n).addElement(name, weight);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------- Thermo --------------------//

    int thermo_newFromFile(const char* filename, const char* phasename) {
        try {
            return ThermoCabinet::add(newPhase(filename, phasename));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_newFromXML(int mxml)
    {
        try {
            XML_Node& x = XmlCabinet::item(mxml);
            thermo_t* th = newPhase(x);
            return ThermoCabinet::add(th);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getEosType(int n, size_t leneos, char* eos)
    {
        try {
            return static_cast<int>(copyString(ThermoCabinet::item(n).type(), eos, leneos));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double thermo_enthalpy_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).enthalpy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_intEnergy_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).intEnergy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_entropy_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).entropy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_gibbs_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).gibbs_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cp_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).cp_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cv_mole(int n)
    {
        try {
            return ThermoCabinet::item(n).cv_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_pressure(int n)
    {
        try {
            return ThermoCabinet::item(n).pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_enthalpy_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_intEnergy_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_entropy_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).entropy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_gibbs_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).gibbs_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cp_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).cp_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_cv_mass(int n)
    {
        try {
            return ThermoCabinet::item(n).cv_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_electricPotential(int n)
    {
        try {
            return ThermoCabinet::item(n).electricPotential();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_chemPotentials(int n, size_t lenm, double* murt)
    {
        try {
            ThermoPhase& thrm = ThermoCabinet::item(n);
            thrm.checkSpeciesArraySize(lenm);
            thrm.getChemPotentials(murt);
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
            ThermoCabinet::item(n).setPressure(p);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_RP(int n, double* vals)
    {
        try{
            ThermoCabinet::item(n).setState_RP(vals[0], vals[1]);
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
            ThermoCabinet::item(n).setState_HP(vals[0],vals[1]);
            if (ThermoCabinet::item(n).temperature() < 0.0) {
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
            ThermoCabinet::item(n).setState_UV(vals[0],vals[1]);
            if (ThermoCabinet::item(n).temperature() < 0.0) {
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
            ThermoCabinet::item(n).setState_SV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_SP(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_SP(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_ST(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_ST(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TV(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_TV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_PV(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_PV(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_UP(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_UP(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_VH(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_VH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_TH(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_TH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_set_SH(int n, double* vals)
    {
        try {
            ThermoCabinet::item(n).setState_SH(vals[0],vals[1]);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_equilibrate(int n, const char* XY, int solver,
                           double rtol, int maxsteps, int maxiter, int loglevel)
    {
        try {
            string ssolver;
            if (solver < 0) {
                ssolver = "auto";
            } else if (solver == 0) {
                ssolver = "element_potential";
            } else if (solver == 1) {
                ssolver = "gibbs";
            } else if (solver == 2) {
                ssolver = "vcs";
            } else {
                throw CanteraError("thermo_equilibrate",
                    "Invalid equilibrium solver specified.");
            }
            ThermoCabinet::item(n).equilibrate(XY, ssolver, rtol, maxsteps,
                                               maxiter, 0, loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal thermo_refPressure(int n)
    {
        try {
            return ThermoCabinet::item(n).refPressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal thermo_minTemp(int n, int k)
    {
        try {
            ThermoPhase& ph = ThermoCabinet::item(n);
            if (k != -1) {
                ph.checkSpeciesIndex(k);
                return ph.minTemp(k);
            } else {
                return ph.minTemp();
            }
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal thermo_maxTemp(int n, int k)
    {
        try {
            ThermoPhase& ph = ThermoCabinet::item(n);
            if (k != -1) {
                ph.checkSpeciesIndex(k);
                return ph.maxTemp(k);
            } else {
                return ph.maxTemp();
            }
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }


    int thermo_getEnthalpies_RT(int n, size_t lenm, double* h_rt)
    {
        try {
            ThermoPhase& thrm = ThermoCabinet::item(n);
            thrm.checkSpeciesArraySize(lenm);
            thrm.getEnthalpy_RT_ref(h_rt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getEntropies_R(int n, size_t lenm, double* s_r)
    {
        try {
            ThermoPhase& thrm = ThermoCabinet::item(n);
            thrm.checkSpeciesArraySize(lenm);
            thrm.getEntropy_R_ref(s_r);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_getCp_R(int n, size_t lenm, double* cp_r)
    {
        try {
            ThermoPhase& thrm = ThermoCabinet::item(n);
            thrm.checkSpeciesArraySize(lenm);
            thrm.getCp_R_ref(cp_r);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setElectricPotential(int n, double v)
    {
        try {
            ThermoCabinet::item(n).setElectricPotential(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal thermo_thermalExpansionCoeff(int n)
    {
        try {
            return ThermoCabinet::item(n).thermalExpansionCoeff();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal thermo_isothermalCompressibility(int n)
    {
        try {
            return ThermoCabinet::item(n).isothermalCompressibility();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    //-------------- pure fluids ---------------//

    double thermo_critTemperature(int n)
    {
        try {
            return ThermoCabinet::item(n).critTemperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_critPressure(int n)
    {
        try {
            return ThermoCabinet::item(n).critPressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_critDensity(int n)
    {
        try {
            return ThermoCabinet::item(n).critDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_vaporFraction(int n)
    {
        try {
            return ThermoCabinet::get<PureFluidPhase>(n).vaporFraction();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_satTemperature(int n, double p)
    {
        try {
            return ThermoCabinet::item(n).satTemperature(p);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double thermo_satPressure(int n, double t)
    {
        try {
            return ThermoCabinet::item(n).satPressure(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int thermo_setState_Psat(int n, double p, double x)
    {
        try {
            ThermoCabinet::get<PureFluidPhase>(n).setState_Psat(p, x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_setState_Tsat(int n, double t, double x)
    {
        try {
            ThermoCabinet::get<PureFluidPhase>(n).setState_Tsat(t, x);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------- Kinetics ------------------//

    int kin_newFromFile(const char* filename, const char* phasename,
                        int reactingPhase, int neighbor1, int neighbor2,
                        int neighbor3, int neighbor4)
    {
        try {
            vector<thermo_t*> phases;
            phases.push_back(&ThermoCabinet::item(reactingPhase));
            if (neighbor1 >= 0) {
                phases.push_back(&ThermoCabinet::item(neighbor1));
                if (neighbor2 >= 0) {
                    phases.push_back(&ThermoCabinet::item(neighbor2));
                    if (neighbor3 >= 0) {
                        phases.push_back(&ThermoCabinet::item(neighbor3));
                        if (neighbor4 >= 0) {
                            phases.push_back(&ThermoCabinet::item(neighbor4));
                        }
                    }
                }
            }
            unique_ptr<Kinetics> kin = newKinetics(phases, filename, phasename);
            return KineticsCabinet::add(kin.release());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_newFromXML(int mxml, int iphase,
                       int neighbor1, int neighbor2, int neighbor3,
                       int neighbor4)
    {
        try {
            XML_Node& x = XmlCabinet::item(mxml);
            vector<thermo_t*> phases;
            phases.push_back(&ThermoCabinet::item(iphase));
            if (neighbor1 >= 0) {
                phases.push_back(&ThermoCabinet::item(neighbor1));
                if (neighbor2 >= 0) {
                    phases.push_back(&ThermoCabinet::item(neighbor2));
                    if (neighbor3 >= 0) {
                        phases.push_back(&ThermoCabinet::item(neighbor3));
                        if (neighbor4 >= 0) {
                            phases.push_back(&ThermoCabinet::item(neighbor4));
                        }
                    }
                }
            }
            Kinetics* kin = newKineticsMgr(x, phases);
            if (kin) {
                return KineticsCabinet::add(kin);
            } else {
                return 0;
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------------------------------
    int kin_getType(int n, size_t lennm, char* nm)
    {
        try {
            return static_cast<int>(copyString(KineticsCabinet::item(n).kineticsType(), nm, lennm));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t kin_start(int n, int p)
    {
        try {
            return KineticsCabinet::item(n).kineticsSpeciesIndex(0,p);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_speciesIndex(int n, const char* nm, const char* ph)
    {
        try {
            return KineticsCabinet::item(n).kineticsSpeciesIndex(nm, ph);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    //---------------------------------------

    size_t kin_nSpecies(int n)
    {
        try {
            return KineticsCabinet::item(n).nTotalSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_nReactions(int n)
    {
        try {
            return KineticsCabinet::item(n).nReactions();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_nPhases(int n)
    {
        try {
            return KineticsCabinet::item(n).nPhases();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_phaseIndex(int n, const char* ph)
    {
        try {
            return KineticsCabinet::item(n).phaseIndex(ph);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t kin_reactionPhaseIndex(int n)
    {
        try {
            return KineticsCabinet::item(n).reactionPhaseIndex();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double kin_reactantStoichCoeff(int n, int k, int i)
    {
        try {
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkSpeciesIndex(k);
            kin.checkReactionIndex(i);
            return kin.reactantStoichCoeff(k,i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double kin_productStoichCoeff(int n, int k, int i)
    {
        try {
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkSpeciesIndex(k);
            kin.checkReactionIndex(i);
            return kin.productStoichCoeff(k,i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int kin_reactionType(int n, int i)
    {
        try {
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkReactionIndex(i);
            return kin.reactionType(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getFwdRatesOfProgress(int n, size_t len, double* fwdROP)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getFwdRatesOfProgress(fwdROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getRevRatesOfProgress(int n, size_t len, double* revROP)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getRevRatesOfProgress(revROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_isReversible(int n, int i)
    {
        try {
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkReactionIndex(i);
            return (int) kin.isReversible(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getNetRatesOfProgress(int n, size_t len, double* netROP)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getNetRatesOfProgress(netROP);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getFwdRateConstants(int n, size_t len, double* kfwd)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getFwdRateConstants(kfwd);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getRevRateConstants(int n, int doIrreversible, size_t len, double* krev)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getRevRateConstants(krev, doIrreversible != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getDelta(int n, int job, size_t len, double* delta)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            switch (job) {
            case 0:
                k.getDeltaEnthalpy(delta);
                break;
            case 1:
                k.getDeltaGibbs(delta);
                break;
            case 2:
                k.getDeltaEntropy(delta);
                break;
            case 3:
                k.getDeltaSSEnthalpy(delta);
                break;
            case 4:
                k.getDeltaSSGibbs(delta);
                break;
            case 5:
                k.getDeltaSSEntropy(delta);
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
            Kinetics& k = KineticsCabinet::item(n);
            k.checkSpeciesArraySize(len);
            k.getCreationRates(cdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getDestructionRates(int n, size_t len, double* ddot)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkSpeciesArraySize(len);
            k.getDestructionRates(ddot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getNetProductionRates(int n, size_t len, double* wdot)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkSpeciesArraySize(len);
            k.getNetProductionRates(wdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getSourceTerms(int n, size_t len, double* ydot)
    {
        try {
            // @todo This function only works for single phase kinetics
            Kinetics& k = KineticsCabinet::item(n);
            ThermoPhase& p = k.thermo();
            size_t nsp = p.nSpecies();
            k.checkSpeciesArraySize(len);
            k.checkSpeciesArraySize(nsp);
            k.getNetProductionRates(ydot);
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
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkReactionIndex(i);
            return kin.multiplier(i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t kin_phase(int n, size_t i)
    {
        try {
            Kinetics& kin = KineticsCabinet::item(n);
            kin.checkPhaseIndex(i);
            return ThermoCabinet::index(kin.thermo(i));
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int kin_getEquilibriumConstants(int n, size_t len, double* kc)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionArraySize(len);
            k.getEquilibriumConstants(kc);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_getReactionString(int n, int i, int len, char* buf)
    {
        try {
            Kinetics& k = KineticsCabinet::item(n);
            k.checkReactionIndex(i);
            return static_cast<int>(copyString(k.reactionString(i), buf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_setMultiplier(int n, int i, double v)
    {
        try {
            if (v >= 0.0) {
                Kinetics& kin = KineticsCabinet::item(n);
                kin.checkReactionIndex(i);
                kin.setMultiplier(i,v);
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
            KineticsCabinet::get<InterfaceKinetics>(n).advanceCoverages(tstep);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //------------------- Transport ---------------------------

    int trans_newDefault(int ith, int loglevel)
    {
        try {
            Transport* tr = newDefaultTransportMgr(&ThermoCabinet::item(ith),
                                                   loglevel);
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_new(const char* model, int ith, int loglevel)
    {
        try {
            Transport* tr = newTransportMgr(model, &ThermoCabinet::item(ith),
                                            loglevel);
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_viscosity(int n)
    {
        try {
            return TransportCabinet::item(n).viscosity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_electricalConductivity(int n)
    {
        try {
            return TransportCabinet::item(n).electricalConductivity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double trans_thermalConductivity(int n)
    {
        try {
            return TransportCabinet::item(n).thermalConductivity();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getThermalDiffCoeffs(int n, int ldt, double* dt)
    {
        try {
            Transport& tr = TransportCabinet::item(n);
            tr.checkSpeciesArraySize(ldt);
            tr.getThermalDiffCoeffs(dt);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMixDiffCoeffs(int n, int ld, double* d)
    {
        try {
            Transport& tr = TransportCabinet::item(n);
            tr.checkSpeciesArraySize(ld);
            tr.getMixDiffCoeffs(d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getBinDiffCoeffs(int n, int ld, double* d)
    {
        try {
            // @todo length of d should be passed for bounds checking
            Transport& tr = TransportCabinet::item(n);
            tr.checkSpeciesArraySize(ld);
            tr.getBinaryDiffCoeffs(ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMultiDiffCoeffs(int n, int ld, double* d)
    {
        try {
            // @todo length of d should be passed for bounds checking
            Transport& tr = TransportCabinet::item(n);
            tr.checkSpeciesArraySize(ld);
            tr.getMultiDiffCoeffs(ld,d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_setParameters(int n, int type, int k, double* d)
    {
        try {
            TransportCabinet::item(n).setParameters(type, k, d);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMolarFluxes(int n, const double* state1,
                             const double* state2, double delta, double* fluxes)
    {
        try {
            TransportCabinet::item(n).getMolarFluxes(state1, state2, delta, fluxes);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_getMassFluxes(int n, const double* state1,
                            const double* state2, double delta, double* fluxes)
    {
        try {
            TransportCabinet::item(n).getMassFluxes(state1, state2, delta, fluxes);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------------- Functions ---------------------------

    int thermo_report(int nth, int ibuf, char* buf, int show_thermo)
    {
        try {
            bool stherm = (show_thermo != 0);
            string s = ThermoCabinet::item(nth).report(stherm);
            if (int(s.size()) > ibuf - 1) {
                return -(static_cast<int>(s.size()) + 1);
            }
            copyString(s, buf, ibuf);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_print(int nth, int show_thermo, double threshold)
    {
        try {
            bool stherm = (show_thermo != 0);
            writelog(ThermoCabinet::item(nth).report(stherm, threshold)+"\n");
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

    int ct_getDataDirectories(int buflen, char* buf, const char* sep)
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

    int ct_clearStorage()
    {
        try {
            ThermoCabinet::clear();
            KineticsCabinet::clear();
            TransportCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int thermo_del(int n)
    {
        try {
            ThermoCabinet::del(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int kin_del(int n)
    {
        try {
            KineticsCabinet::del(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int trans_del(int n)
    {
        try {
            TransportCabinet::del(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_ck2cti(const char* in_file, const char* db_file, const char* tr_file,
                  const char* id_tag, int debug, int validate)
    {
        try {
            ck2cti(in_file, db_file, tr_file, id_tag);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
