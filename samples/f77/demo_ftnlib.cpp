/*!
      A simple Fortran 77 interface


 This file is an example of how to write an interface to use Cantera
 in Fortran 77 programs. The basic idea is to store pointers to
 Cantera objects in global storage, and then create Fortran-callable
 functions that access the objects through the pointers.

 This particular example defines functions that return thermodynamic
 properties, transport properties, and kinetic rates for reacting
 ideal gas mixtures. Only a single pointer to an IdealGasPhase object is
 stored, so only one reaction mechanism may be used at any one time in
 the application.  Of course, it is a simple modification to store
 multiple objects if it is desired to use multiple reaction
 mechanisms.

 The functions defined here are ones commonly needed in application
 programs that simulate gas-phase combustion or similar
 processes. Similar libraries to access other capabilities of Cantera
 (surface chemistry, etc.) could be written in the same way.

 This library is designed for Fortran compilers that expect external
 procedure names to be lowercase with a trailing underscore. If this
 is not the case, the procedure names must be edited before use.

 */

// add any other Cantera header files you need here
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"

#include <iostream>

using namespace Cantera;
using std::string;

// store a pointer to a Solution object
// provides access to the pointers for functions in other libraries
static shared_ptr<Solution> _sol = NULL;

// store a pointer to the thermophase object
static shared_ptr<ThermoPhase> _gas = NULL;
shared_ptr<ThermoPhase> _gasptr()
{
    return _gas;
}

// store a pointer to the kinetics object
static shared_ptr<Kinetics> _kin = NULL;
shared_ptr<Kinetics> _kinptr()
{
    return _kin;
}

// store a pointer to a transport manager
static shared_ptr<Transport> _trans = NULL;
shared_ptr<Transport> _transptr()
{
    return _trans;
}

// error handler
void handleError(CanteraError& err)
{
    std::cout << err.what() << std::endl;
    exit(-1);
}

// extern "C" turns off C++ name-mangling, so that the procedure names
// in the object file are exactly as shown here.

extern "C" {

    /// This is the Fortran main program. This works for g77; it may
    /// need to be modified for other Fortran compilers
#ifdef NEED_ALT_MAIN
    extern int MAIN__();
#endif

    /**
     * Read in a reaction mechanism file and create a Solution
     * object. The file may be in Cantera input format or in YAML. (If
     * you have a file in Chemkin-compatible format, use utility
     * program ck2yaml first to convert it into Cantera format.)
     */
    void newidealgasmix_(char* file, char* id, char* transport,
                         ftnlen lenfile, ftnlen lenid, ftnlen lentr)
    {
        string trmodel = "";
        try {
            string fin = string(file, lenfile);
            string fth = string(id, lenid);
            trmodel = string(transport, lentr);
            _sol = newSolution(fin, fth, trmodel);
            _gas = _sol->thermo();
            _kin = _sol->kinetics();
            _trans = _sol->transport();
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    /// integer function nElements()
    integer nelements_()
    {
        return _gas->nElements();
    }

    /// integer function nSpecies()
    integer nspecies_()
    {
        return _gas->nSpecies();
    }

    /// integer function nReactions()
    integer nreactions_()
    {
        return _kin->nReactions();
    }

    void getspeciesname_(integer* k, char* name, ftnlen n)
    {
        int ik = *k - 1;
        std::fill(name, name + n, ' ');
        string spnm = _gas->speciesName(ik);
        int ns = spnm.size();
        unsigned int nmx = (ns > n ? n : ns);
        copy(spnm.begin(), spnm.begin()+nmx, name);
    }

    //-------------- setting the state ----------------------------

    /// subroutine setState_TPX(T, P, X)
    void setstate_tpx_(double* T, double* P, double* X)
    {
        try {
            _gas->setState_TPX(*T, *P, X);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    /// subroutine setState_TPX_String(T, P, X)
    void setstate_tpx_string_(double* T, double* P,
                              char* X, ftnlen lenx)
    {
        try {
            _gas->setState_TPX(*T, *P, string(X, lenx));
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    void setstate_try_(double* T, double* rho, double* Y)
    {
        try {
            _gas->setState_TRY(*T, *rho, Y);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    void setstate_tpy_(double* T, double* p, double* Y)
    {
        try {
            _gas->setState_TPY(*T, *p, Y);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    void setstate_sp_(double* s, double* p)
    {
        try {
            _gas->setState_SP(*s, *p);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    //-------------- thermodynamic properties ----------------------

    /// Temperature (K)
    double temperature_()
    {
        return _gas->temperature();
    }

    /// Pressure (Pa)
    double pressure_()
    {
        return _gas->pressure();
    }

    /// Density (kg/m^3)
    double density_()
    {
        return _gas->density();
    }

    /// Mean molar mass (kg/kmol).
    double meanmolarmass_()
    {
        return _gas->meanMolecularWeight();
    }

    /// Molar enthalpy (J/kmol)
    double enthalpy_mole_()
    {
        return _gas->enthalpy_mole();
    }

    /// Molar internal energy (J/kmol)
    double intenergy_mole_()
    {
        return _gas->intEnergy_mole();
    }

    /// Molar entropy (J/kmol-K)
    double entropy_mole_()
    {
        return _gas->entropy_mole();
    }

    /// Molar heat capacity at constant P (J/kmol-K)
    double cp_mole_()
    {
        return _gas->cp_mole();
    }

    /// Molar Gibbs function (J/kmol)
    double gibbs_mole_()
    {
        return _gas->gibbs_mole();
    }

    double enthalpy_mass_()
    {
        return _gas->enthalpy_mass();
    }

    double intenergy_mass_()
    {
        return _gas->intEnergy_mass();
    }

    double entropy_mass_()
    {
        return _gas->entropy_mass();
    }

    double cp_mass_()
    {
        return _gas->cp_mass();
    }

    double cv_mass_()
    {
        return _gas->cv_mass();
    }

    double gibbs_mass_()
    {
        return _gas->gibbs_mass();
    }

    void gotmolefractions_(double* x)
    {
        _gas->getMoleFractions(x);
    }

    void gotmassfractions_(double* y)
    {
        _gas->getMassFractions(y);
    }

    void equilibrate_(char* opt, ftnlen lenopt)
    {
        try {
            if (lenopt != 2) {
                throw CanteraError("equilibrate",
                                   "two-character string required.");
            }
            string optstr = string(opt, 2);
            _gas->equilibrate(optstr);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    //---------------- kinetics -------------------------

    void getreactioneqn_(integer* i, char* eqn, ftnlen n)
    {
        int irxn = *i - 1;
        std::fill(eqn, eqn + n, ' ');
        string e = _kin->reactionString(irxn);
        int ns = e.size();
        unsigned int nmx = (ns > n ? n : ns);
        copy(e.begin(), e.begin()+nmx, eqn);
    }

    void getnetproductionrates_(double* wdot)
    {
        _kin->getNetProductionRates(wdot);
    }

    void getcreationrates_(double* cdot)
    {
        _kin->getCreationRates(cdot);
    }

    void getdestructionrates_(double* ddot)
    {
        _kin->getDestructionRates(ddot);
    }

    void getnetratesofprogress_(double* q)
    {
        _kin->getNetRatesOfProgress(q);
    }

    void getfwdratesofprogress_(double* q)
    {
        _kin->getFwdRatesOfProgress(q);
    }

    void getrevratesofprogress_(double* q)
    {
        _kin->getRevRatesOfProgress(q);
    }

    //-------------------- transport properties --------------------

    double viscosity_()
    {
        try {
            return _trans->viscosity();
        } catch (CanteraError& err) {
            handleError(err);
            return 0.0;
        }
    }

    double thermalconductivity_()
    {
        try {
            return _trans->thermalConductivity();
        } catch (CanteraError& err) {
            handleError(err);
            return 0.0;
        }
    }

    void getmixdiffcoeffs_(double* diff)
    {
        try {
            _trans->getMixDiffCoeffs(diff);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

    void getthermaldiffcoeffs_(double* dt)
    {
        try {
            _trans->getThermalDiffCoeffs(dt);
        } catch (CanteraError& err) {
            handleError(err);
        }
    }

}

/*
 *  HKM 7/22/09:
 *    I'm skeptical that you need this for any system.
 *    Definitely creates an error (dupl main()) for the solaris
 *    system
 */
#ifdef NEED_ALT_MAIN
/**
 * This C++ main program simply calls the Fortran main program.
 */
int main()
{
    try {
        return MAIN__();
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        exit(-1);
    } catch (...) {
        cout << "An exception was trapped. Program terminating." << endl;
        exit(-1);
    }
}
#endif
