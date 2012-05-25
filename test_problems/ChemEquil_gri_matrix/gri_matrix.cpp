/*
 *  Copyright 2002 California Institute of Technology
 */

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    int numSucc = 0;
    int numFail = 0;
    try {
        int retnSub;
        double T = 500.;
        int solver = -1;

        IdealGasMix g("gri30.xml", "gri30_mix");
        double pres = OneAtm;

        size_t kk = g.nSpecies();
        vector_fp Xmol(kk, 0.0);
        size_t iCH4 = g.speciesIndex("CH4");
        size_t iO2 = g.speciesIndex("O2");
        size_t iN2 = g.speciesIndex("N2");

        /*
         * Do an initial calculation that can be debugged
         * easily
         */
        Xmol[iCH4] = 0.6;
        Xmol[iO2] = 0.0;
        Xmol[iN2] = 0.4;
        g.setState_TPX(T, pres, DATA_PTR(Xmol));
        try {
            retnSub = equilibrate(g, "TP", solver);
            cout << g;
            if (retnSub != 1) {
                cerr << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                cout << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                exit(-1);
            }
        } catch (CanteraError& err) {
            cout << g;
            std::cerr << err.what() << std::endl;
            cerr << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            cout << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            exit(-1);
        }

        /*
         * Do an initial calculation that can be debugged
         * easily
         */
        Xmol[iCH4] = 0.0;
        Xmol[iO2] = 0.6;
        Xmol[iN2] = 0.4;
        g.setState_TPX(T, pres, DATA_PTR(Xmol));
        try {
            retnSub = equilibrate(g, "TP", solver);
            cout << g;
            if (retnSub != 1) {
                cerr << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                cout << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                exit(-1);
            }
        } catch (CanteraError& err) {
            cout << g;
            std::cerr << err.what() << std::endl;
            cerr << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            cout << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            exit(-1);
        }

        /*
          * Do an initial calculation that can be debugged
          * easily
          */
        Xmol[iCH4] = 0.3;
        Xmol[iO2] = 0.3;
        Xmol[iN2] = 0.4;
        T = 2000.;
        pres = OneAtm;
        g.setState_TPX(T, pres, DATA_PTR(Xmol));
        try {
            retnSub = equilibrate(g, "TP", solver);
            cout << g;
            if (retnSub != 1) {
                cerr << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                cout << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                exit(-1);
            }
        } catch (CanteraError& err) {
            cout << g;
            std::cerr << err.what() << std::endl;
            cerr << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            cout << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            exit(-1);
        }

        /*
           * Do an initial calculation that can be debugged
           * easily
           */
        Xmol[iCH4] = 0.3;
        Xmol[iO2] = 0.3;
        Xmol[iN2] = 0.0;
        T = 2000.;
        pres = 1.0;
        g.setState_TPX(T, pres, DATA_PTR(Xmol));
        try {
            retnSub = equilibrate(g, "TP", solver);
            cout << g;
            if (retnSub != 1) {
                cerr << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                cout << "ERROR: ChemEquil equilibration step failed at "
                     << " T    = " << T
                     << " Pres = " << pres
                     << endl;
                exit(-1);
            }
        } catch (CanteraError& err) {
            cout << g;
            std::cerr << err.what() << std::endl;
            cerr << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            cout << "ERROR: equilibration step failed at "
                 << " T    = " << T
                 << " Pres = " << pres
                 << endl;
            exit(-1);
        }

        /*
         * OK do the matrix.
         */
        bool showResults = false;
        for (int iS = 0; iS < 6; iS++) {
            Xmol[iCH4] = iS * 0.6 / 5.0;
            Xmol[iO2] = 1.0 - Xmol[iCH4] - Xmol[iN2];
            for (int iP = 0; iP < 10; iP++) {
                pres = pow(10.0, iP) *1.0E-2;
                for (int iT = 0; iT < 25; iT++) {
                    double T = 100. * iT + 300.;

                    /*
                     * Reset mole fractions to a base state
                     */
                    g.setState_TPX(T, pres, DATA_PTR(Xmol));

                    retnSub = 0;
                    try {
                        retnSub = equilibrate(g, "TP", solver);
                        if (retnSub != 1) {
                            cerr << "ERROR: ChemEquil equilibration step failed at "
                                 << " T    = " << T
                                 << " Pres = " << pres
                                 << endl;
                            cout << "ERROR: ChemEquil equilibration step failed at "
                                 << " T    = " << T
                                 << " Pres = " << pres
                                 << endl;
                        }
                    } catch (CanteraError& err) {
                        std::cerr << err.what() << std::endl;
                        cerr << "ERROR: equilibration step failed at "
                             << " T    = " << T
                             << " Pres = " << pres
                             << endl;
                        cout << "ERROR: equilibration step failed at "
                             << " T    = " << T
                             << " Pres = " << pres
                             << endl;
                    }
                    if (showResults) {
                        cout << g;
                    }
                    if (retnSub ==1) {
                        numSucc++;
                    } else {
                        numFail++;
                    }

                }
            }
        }

        cout << "NUMBER OF SUCCESSES =  " << numSucc << endl;
        cout << "NUMBER OF FAILURES  =  " << numFail << endl;

        return numFail;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "ERROR: program terminating due to unforeseen circumstances." << endl;
        return -1;
    }
}
