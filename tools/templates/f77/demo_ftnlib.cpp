/*!
      A simple Fortran 77 interface

 This file is an example of how to write an interface to use Cantera
 in Fortran 77 programs. The basic idea is to store pointers to
 Cantera objects in global storage, and then create Fortran-callable
 functions that access the objects through the pointers.

 This particular example defines functions that return thermodynamic
 properties and kinetic rates for reacting ideal gas mixtures. Only a
 single pointer to an IdealGasMix object is stored, so only one
 reaction mechanism may be used at any one time in the application.
 Of course, it is a simple modification to store multiple objects if 
 it is desired to use multiple reaction mechanisms.

 The functions defined here are ones commonly needed in application
 programs that simulate gas-phase combustion or similar processes.

 */

// add any other Cantera header files you need here
#include "IdealGasMix.h"
#include "equilibrium.h"

// store a pointer to an IdealGasMix object. The object itself will
// be created by the call to init_.
static IdealGasMix* _gas = 0;


// provides access to the pointer for functions in other libraries
IdealGasMix* _gasptr() { return _gas; }

// error handler 
void handleError() {
    showErrors(cout);
    exit(-1);
}

// extern "C" turns off C++ name-mangling, so that the procedure names
// in the object file are exactly as shown here.

extern "C" {

    /// This is the Fortran main program
    extern int MAIN__();

    /**
     * Read in a reaction mechanism file and create an IdealGasMix
     * object. The file may be in Cantera input format or in CTML. (If
     * you have a file in Chemkin-compatible format, use utility
     * program ck2cti first to convert it into Cantera format.)
     */
    void newidealgasmix_(char* file, char* id, 
        ftnlen lenfile, ftnlen lenid) {
        try {
            string fin = string(file, lenfile);
            string fth = string(id, lenid);
            if (_gas) delete _gas;
            _gas = new IdealGasMix(fin, fth);
        }
        catch (CanteraError) {
            handleError();
        }
    }

    ///   integer function nElements() 
    integer nelements_() { return _gas->nElements(); }

    /// integer function nSpecies()
    integer nspecies_() { return _gas->nSpecies(); }

    /// integer function nReactions()
    integer nreactions_() { return _gas->nReactions(); }


    //-------------- setting the state ----------------------------

    // subroutine setState_TPX(T, P, X)
    void setstate_tpx_(doublereal* T, doublereal* P, doublereal* X) {
        try {
            _gas->setState_TPX(*T, *P, X);
        }
        catch (CanteraError) { handleError(); }
    } 

    /// subroutine setState_TPX_String(T, P, X)
    void setstate_tpx_string_(doublereal* T, doublereal* P, 
        char* X, ftnlen lenx) {
        try {
            _gas->setState_TPX(*T, *P, string(X, lenx));
        }
        catch (CanteraError) { handleError(); }
    } 

    void setstate_try_(doublereal* T, doublereal* rho, doublereal* Y) {
        try {
            _gas->setState_TRY(*T, *rho, Y);
        }
        catch (CanteraError) { handleError(); }
    } 

    //-------------- thermodynamic properties ----------------------

    /// Temperature (K)
    doublereal temperature_() { 
        return _gas->temperature();
    }

    /// Pressure (Pa)
    doublereal pressure_() { 
        return _gas->pressure();
    }
    
    /// Density (kg/m^3)
    doublereal density_() { 
        return _gas->density();
    }

    /// Mean molar mass (kg/kmol).
    doublereal meanmolarmass_() { 
        return _gas->meanMolecularWeight();
    }

    /// Molar enthalpy (J/kmol) 
    doublereal enthalpy_mole_() { 
        return _gas->enthalpy_mole();
    }

    /// Molar internal energy (J/kmol)
    doublereal intenergy_mole_() { 
        return _gas->intEnergy_mole();
    }
    
    /// Molar entropy (J/kmol-K)
    doublereal entropy_mole_() { 
        return _gas->entropy_mole();
    }

    /// Molar heat capacity at constant P (J/kmol-K)
    doublereal cp_mole_() { 
        return _gas->cp_mole();
    }

    /// Molar Gibbs function (J/kmol)
    doublereal gibbs_mole_() { 
        return _gas->gibbs_mole();
    }

    doublereal enthalpy_mass_() { 
        return _gas->enthalpy_mass();
    }

    doublereal intenergy_mass_() { 
        return _gas->intEnergy_mass();
    }
    
    doublereal entropy_mass_() { 
        return _gas->entropy_mass();
    }

    doublereal cp_mass_() { 
        return _gas->cp_mass();
    }

    doublereal gibbs_mass_() { 
        return _gas->gibbs_mass();
    }
    
    void gotmolefractions_(doublereal* x) {
        _gas->getMoleFractions(x);
    }

    void gotmassfractions_(doublereal* y) {
        _gas->getMassFractions(y);
    }

    void equilibrate_(char* opt, ftnlen lenopt) {
        try {
            if (lenopt != 2) {
                throw CanteraError("equilibrate",
                    "two-character string required.");
            }
            string optstr = string(opt, 2);
            equilibrate(*_gas, optstr.c_str());
        }
        catch (CanteraError) { handleError(); }
    }


    //---------------- kinetics -------------------------

    void getreactioneqn_(integer* i, char* eqn, ftnlen n) {
        int irxn = *i - 1;
        fill(eqn, eqn + n, ' ');
        string e = _gas->reactionString(irxn);
        int ns = e.size();
        unsigned int nmx = (ns > n ? n : ns);
        copy(e.begin(), e.begin()+nmx, eqn);
    }
        
    void getnetproductionrates_(doublereal* wdot) {
        _gas->getNetProductionRates(wdot);
    }

    void getcreationrates_(doublereal* cdot) {
        _gas->getCreationRates(cdot);
    }

    void getdestructionrates_(doublereal* ddot) {
        _gas->getDestructionRates(ddot);
    }

    void getnetratesofprogress_(doublereal* q) {
        _gas->getNetRatesOfProgress(q);
    }

    void getfwdratesofprogress_(doublereal* q) {
        _gas->getFwdRatesOfProgress(q);
    }

    void getrevratesofprogress_(doublereal* q) {
        _gas->getRevRatesOfProgress(q);
    }

}


/**
 * This C++ main program simply calls the Fortran main program. 
 */
int main() {
    try {
        return MAIN__();
    }
    catch (CanteraError) {
        showErrors(cerr);
        exit(-1);
    }
    catch (...) {
        cout << "An exception was trapped. Program terminating." << endl;
        exit(-1);
    }
}

