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

 The functions defined here are ones commonly needed in application
 programs that simulate gas-phase combustion or similar processes.
 You may find that you can use this library without modification.  If
 you need to access additional capabilities of Cantera from Fortran,
 you can use this file as a guide. Also take a look at the Fortran 77
 examples in the demos/f77 subdirectory within the directory where
 Cantera is installed.

 */

// add any other Cantera header files you need here
#include "IdealGasMix.h"
#include "equilibrium.h"


// store a pointer to an IdealGasMix object. The object itself will
// be created by the call to init_.
static IdealGasMix* _gas = 0;



extern "C" {

    /// This is the Fortran main program
    extern int MAIN__();


    /**
     * Read in a reaction mechanism file and create an IdealGasMix
     * object. The file may be in Chemkin-compatible format or in
     * CTML. The name of a thermodynamic database may be supplied as a
     * second argument. If none is required, enter an empty string as
     * the second argument.
     */
    void readmechanism_(char* file, char* thermo, 
        ftnlen lenfile, ftnlen lenthermo) {
        string fin = string(file, lenfile);
        string fth = string(thermo, lenthermo);
        if (_gas) delete _gas;
        _gas = new IdealGasMix(fin, fth);
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
        _gas->setState_TPX(*T, *P, X);
    } 

    /// subroutine setState_TPX_AsString(T, P, X)
    void setstate_tpx_asstring_(doublereal* T, doublereal* P, 
        char* X, ftnlen lenx) {
        _gas->setState_TPX(*T, *P, string(X, lenx));
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
        return _gas->enthalpy_mole();
    }

    doublereal intenergy_mass_() { 
        return _gas->intEnergy_mole();
    }
    
    doublereal entropy_mass_() { 
        return _gas->entropy_mole();
    }

    doublereal cp_mass_() { 
        return _gas->cp_mole();
    }

    doublereal gibbs_mass_() { 
        return _gas->gibbs_mole();
    }

    
    void equilibrate_(integer* opt) {
        int option = *opt;
        equilibrate(*_gas, option);
    }


    //---------------- kinetics -------------------------

    void getreactioneqn_(integer* i, char* eqn, ftnlen n) {
        int irxn = *i - 1;
        fill(eqn, eqn + n, ' ');
        string e = _gas->reactionString(irxn);
        unsigned int nmx = (e.size() > n ? n : e.size());
        copy(e.begin(), e.begin()+nmx, eqn);
    }
        
    void getnetproductionrates_(doublereal* wdot) {
        _gas->getNetProductionRates(wdot);
    }

    void getnetratesofprogress_(doublereal* q) {
        _gas->getNetRatesOfProgress(q);
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
        return -1;
    }
    catch (...) {
        cout << "An exception was trapped. Program terminating." << endl;
        return -1;
    }
}

