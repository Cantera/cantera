/**
 *  @file surfaceSolver2.cpp
 *
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
//  Example
//
//  Read a surface growth mechanism and calculate the solution
//  using Placid.
//

#define MSSIZE 200

#ifdef DEBUG_HKM
int iDebug_HKM = 0;
#endif

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{

}

#include "cantera/Interface.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/kinetics/solveSP.h"
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Cantera;

void printGas(ostream& oooo, ThermoPhase* gasTP, InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    double C[MSSIZE];
    char buf[256];
    oooo.precision(3);
    string gasPhaseName          = "gas";
    gasTP->getMoleFractions(x);
    gasTP->getConcentrations(C);
    double Temp = gasTP->temperature();
    double p = gasTP->pressure();
    oooo << "Gas Temperature = " << Temp << endl;
    oooo << "Gas Pressure    = " << p << endl;
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, 0);
    oooo << "Gas Phase:  " << gasPhaseName << "   "
         << "(" << kstart << ")" << endl;
    oooo << "                       Name      "
         << "     Conc              MoleF       SrcRate " << endl;
    oooo << "                                 "
         << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
    double sum = 0.0;
    size_t nspGas = gasTP->nSpecies();
    for (size_t k = 0; k < nspGas; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(k, 0);
        sprintf(buf, "%4d %24s   %14g %14g  %14e\n",
                (int) k, gasTP->speciesName(k).c_str(),
                C[k], x[k], src[kstart]);
        oooo << buf;
        sum += x[k];
    }
    oooo << "Sum of gas mole fractions= " << sum << endl;
    oooo << endl;


}

void printBulk(ostream& oooo,ThermoPhase* bulkPhaseTP, InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    double C[MSSIZE];
    char buf[256];
    oooo.precision(3);
    string bulkParticlePhaseName = bulkPhaseTP->id();
    bulkPhaseTP->getMoleFractions(x);
    bulkPhaseTP->getConcentrations(C);
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, 1);
    double dens = bulkPhaseTP->density();
    oooo << "Bulk Phase:  " << bulkParticlePhaseName << "   "
         << "(" << kstart << ")" << endl;
    double Temp = bulkPhaseTP->temperature();
    double p = bulkPhaseTP->pressure();
    oooo << "Bulk Temperature = " << Temp << endl;
    oooo << "Bulk Pressure    = " << p << endl;
    oooo << "                       Name      "
         << "     Conc              MoleF       SrcRate " << endl;
    oooo << "                                 "
         << "   (kmol/m^3)                   (kmol/m^2/s) " << endl;
    double sum = 0.0;
    double Wsum = 0.0;
    const vector_fp& molecW = bulkPhaseTP->molecularWeights();
    size_t nspBulk = bulkPhaseTP->nSpecies();
    for (size_t k = 0; k < nspBulk; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(k, 1);
        sprintf(buf, "%4d %24s   %14g %14g  %14e\n",
                (int) k, bulkPhaseTP->speciesName(k).c_str(),
                C[k], x[k], src[kstart]);
        oooo << buf;
        sum += x[k];
        Wsum += src[kstart] * molecW[k];
    }
    oooo << "Bulk Weight Growth Rate = " << Wsum << " kg/m^2/s" << endl;
    double gr = Wsum / dens;
    oooo << "Bulk Growth Rate = " << gr << " m/s" << endl;
    oooo << "Bulk Growth Rate = " << gr * 1.0E6 * 3600.
         << " microns / hour" << endl;
    oooo << "Density of bulk phase = " << dens << " kg / m^3 "<< endl;
    oooo << "                      = " << dens / 1.0E3
         <<" gm / cm^3 " << endl;
    oooo << "Sum of bulk mole fractions= " << sum << endl;
    oooo << endl;
}

void printSurf(ostream& oooo, ThermoPhase* surfPhaseTP,
               InterfaceKinetics* iKin_ptr, double* src)
{
    double x[MSSIZE];
    char buf[256];
    oooo.precision(3);
    string surfParticlePhaseName = surfPhaseTP->id();
    surfPhaseTP->getMoleFractions(x);
    size_t kstart = iKin_ptr->kineticsSpeciesIndex(0, 2);
    oooo << "Surface Phase:  " << surfParticlePhaseName
         << " (" << kstart << ")" << endl;
    double Temp = surfPhaseTP->temperature();
    double p = surfPhaseTP->pressure();
    oooo << "Surface Temperature = " << Temp << endl;
    oooo << "Surface Pressure    = " << p << endl;
    oooo << "                       Name      "
         << "   Coverage         SrcRate " << endl;
    double sum = 0.0;
    size_t nspSurf = surfPhaseTP->nSpecies();
    for (size_t k = 0; k < nspSurf; k++) {
        kstart = iKin_ptr->kineticsSpeciesIndex(0, 2);
        double srcK = src[kstart];
        if (fabs(srcK) < 1.0E-8) {
            srcK = 0.0;
        }
        sprintf(buf, "%4d %24s   %14g   %14e\n",
                (int) k, surfPhaseTP->speciesName(k).c_str(),
                x[k], srcK);
        oooo << buf;
        sum += x[k];
    }
    oooo << "Sum of coverages = " << sum << endl;
}

int main(int argc, char** argv)
{
    string infile;
    int i, k;
    int ioflag = 1;
    // look for command-line options
    if (argc > 1) {
        string tok;
        for (int j = 1; j < argc; j++) {
            tok = string(argv[j]);
            if (tok[0] == '-') {
                int nopt = static_cast<int>(tok.size());
                for (int n = 1; n < nopt; n++) {
                    if (tok[n] == 'h') {
                        printUsage();
                        exit(0);
                    } else if (tok[n] == 'd') {
                        int lvl = 0;
                        if (j < (argc - 1)) {
                            string tokla = string(argv[j+1]);
                            if (strlen(tokla.c_str()) > 0) {
                                lvl = atoi(tokla.c_str());
                                n = nopt - 1;
                                j += 1;
                                ioflag = lvl;
                            }
                        }
                    } else {
                        printUsage();
                        exit(1);
                    }
                }
            } else if (infile == "") {
                infile = tok;
            } else {
                printUsage();
                exit(1);
            }
        }
    }
    if (infile == "") {
        infile = "diamond.cti";
    }

    try {
        /*************************************************************/

        /*
         *  FILL IN THESE NAMES FOR EACH PROBLEM
         */
        /*
         * ProblemNumber = 0 : diamond.cti
         *               = 1 : haca.cti
         */
        int ProblemNumber = 1;
        string gasPhaseName          = "gas";
        string bulkParticlePhaseName = "diamond";
        string surfParticlePhaseName = "diamond_100";
        if (ProblemNumber == 1) {
            gasPhaseName          = "gas";
            bulkParticlePhaseName = "soot";
            surfParticlePhaseName = "soot_interface";
        }

        /************************************************************/
        XML_Node* xc = new XML_Node();
        string path = findInputFile(infile);
        ctml::get_CTML_Tree(xc, path);

        XML_Node* const xg = (XML_Node*) findXMLPhase(xc, gasPhaseName);
        if (!xg) {
            printf("ERROR: Could not find gas phase named, %s, in file\n",
                   gasPhaseName.c_str());
            exit(-1);
        }
        ThermoPhase* gasTP = newPhase(*xg);
        size_t nspGas = gasTP->nSpecies();
        cout << "Number of species = " << nspGas << endl;

        XML_Node* const xd =
            (XML_Node*) findXMLPhase(xc, bulkParticlePhaseName);
        if (!xd) {
            printf("ERROR: Could not find bulk phase named, %s, in file\n",
                   bulkParticlePhaseName.c_str());
            exit(-1);
        }
        ThermoPhase* bulkPhaseTP = newPhase(*xd);
        size_t nspBulk = bulkPhaseTP->nSpecies();
        cout << "Number of species in bulk phase named " <<
             bulkParticlePhaseName << " = " << nspBulk << endl;


        XML_Node* const xs =
            (XML_Node*) findXMLPhase(xc, surfParticlePhaseName);
        if (!xs) {
            printf("ERROR: Could not find surf Particle phase named,"
                   "%s, in file\n",
                   surfParticlePhaseName.c_str());
            exit(-1);
        }
        ThermoPhase* surfPhaseTP = newPhase(*xs);
        size_t nsp_d100 = surfPhaseTP->nSpecies();
        cout << "Number of species in surface phase, " << surfParticlePhaseName
             << " = " << nsp_d100 << endl;

        vector<ThermoPhase*> phaseList;
        phaseList.push_back(gasTP);
        phaseList.push_back(bulkPhaseTP);
        phaseList.push_back(surfPhaseTP);



        InterfaceKinetics* iKin_ptr = new InterfaceKinetics();
        importKinetics(*xs, phaseList, iKin_ptr);
        size_t nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        ofstream ofile("results2.txt");


        // create a second copy of the same surface phase
        // (this is a made up problem btw to check the software capability)
        ThermoPhase* surfPhaseTP2 = newPhase(*xs);
        size_t nsp2 = surfPhaseTP2->nSpecies();
        string pname = surfPhaseTP2->id();
        cout << "Number of species in 2nd surface phase, " << pname
             << " = " << nsp2 << endl;

        vector<ThermoPhase*> phaseList2;
        phaseList2.push_back(gasTP);
        phaseList2.push_back(bulkPhaseTP);
        phaseList2.push_back(surfPhaseTP2);

        // create the second  InterfaceKinetics object based on the
        // second surface phase.
        InterfaceKinetics* iKin2_ptr = new InterfaceKinetics();
        importKinetics(*xs, phaseList2, iKin2_ptr);
        nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        double x[MSSIZE], p = OneAtm;

        /*
         * Set the Gas State:
         * -> note that the states are set in the xml files too
         */
        for (i = 0; i < MSSIZE; i++) {
            x[i] = 0.0;
        }
        if (ProblemNumber == 0) {
            x[0] = 0.0010;
            x[1] = 0.9888;
            x[2] = 0.0002;
            x[3] = 0.0100;
            p = 20.0*OneAtm/760.0;
            gasTP->setState_TPX(1200., p, x);
        }

        /*
         * Set the surface initial state
         *  other problem numbers take their initial state from the xml files.
         */
        for (i = 0; i < MSSIZE; i++) {
            x[i] = 0.0;
        }
        if (ProblemNumber == 0) {
            size_t i0 = surfPhaseTP->speciesIndex("c6H*");
            if (i0 != npos) {
                x[i0] = 0.1;
            }
            size_t i1 = surfPhaseTP->speciesIndex("c6HH");
            if (i1 != npos) {
                x[i1] = 0.9;
            }
            surfPhaseTP->setState_TX(1200., x);
        }

        /*
         * Set the bulk Phase State
         */
        for (i = 0; i < MSSIZE; i++) {
            x[i] = 0.0;
        }
        if (ProblemNumber == 0) {
            x[0] = 1.0;
            bulkPhaseTP->setState_TPX(1200., p, x);
        }

        /*
         *  Set-up the Surface Problem
         *    This problem will consist of 2 identical InterfaceKinetics objects
         */
        vector<InterfaceKinetics*> vecKinPtrs;
        vecKinPtrs.push_back(iKin_ptr);
        vecKinPtrs.push_back(iKin2_ptr);

        // Create the ImplicitSurfChem problem
        // Initialize it and call the pseudo steadystate capability.
        ImplicitSurfChem* surfaceProb = new ImplicitSurfChem(vecKinPtrs);
        surfaceProb->initialize();
        surfaceProb->setIOFlag(ioflag);
        surfaceProb->solvePseudoSteadyStateProblem();

        /*
         * Download the source terms for the rate equations
         */
        double src[MSSIZE];
        double src2[MSSIZE];
        iKin_ptr->getNetProductionRates(src);
        iKin2_ptr->getNetProductionRates(src2);

        double sum = 0.0;
        if (ProblemNumber == 0) {
            double naH;
            for (k = 0; k < 13; k++) {
                if (k < 4) {
                    naH = gasTP->nAtoms(k, 0);
                } else if (k == 4) {
                    naH = 0;
                } else if (k > 4) {
                    int itp = k - 5;
                    naH = surfPhaseTP->nAtoms(itp, 0);
                }
                cout << k << "  " << naH << "  " ;
                if (fabs(src[k]) < 2.0E-17) {
                    cout << " nil" << endl;
                } else {
                    cout << src[k] << endl;
                }
                sum += naH * src[k];
            }
            cout << "sum = " << sum << endl;
        }


        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;
        printSurf(cout, surfPhaseTP2, iKin2_ptr, src2) ;


        printGas(ofile, gasTP, iKin_ptr, src);
        printBulk(ofile, bulkPhaseTP, iKin_ptr, src);
        printSurf(ofile, surfPhaseTP, iKin_ptr, src) ;
        printSurf(ofile, surfPhaseTP2, iKin2_ptr, src2) ;

        /*****************************************************************************/
        /*  Now Tweak the inputs and do a quick calculation */
        /****************************************************************************/

        /*
         * Set the Gas State:
         * -> note that the states are set in the xml files too
         */
        double pres = gasTP->pressure();
        gasTP->getMoleFractions(x);
        double tmp = 0.3 * std::min(x[0], x[1]);
        x[0] += tmp;
        x[1] -= tmp;
        gasTP->setState_PX(pres, x);

        surfaceProb->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);
        iKin2_ptr->getNetProductionRates(src2);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;
        printSurf(cout, surfPhaseTP2, iKin2_ptr, src2) ;

        /*****************************************************************************/
        /*  Now Tweak the inputs and do a quick calculation */
        /****************************************************************************/

        /*
         * Set the Gas State:
         * -> note that the states are set in the xml files too
         */

        /*
         * Set the Gas State:
         * -> note that the states are set in the xml files too
         */
        pres = gasTP->pressure();
        double temp = gasTP->temperature();
        temp += 95;
        gasTP->setState_TP(temp, pres);

        surfaceProb->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);
        iKin2_ptr->getNetProductionRates(src2);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;
        printSurf(cout, surfPhaseTP2, iKin2_ptr, src2) ;

        /*****************************************************************************/
        /*  Now Don't Tweak the inputs at all */
        /****************************************************************************/
        gasTP->setState_TP(temp, pres);

        surfaceProb->solvePseudoSteadyStateProblem();
        iKin_ptr->getNetProductionRates(src);
        iKin2_ptr->getNetProductionRates(src2);

        printGas(cout, gasTP, iKin_ptr, src);
        printBulk(cout, bulkPhaseTP, iKin_ptr, src);
        printSurf(cout, surfPhaseTP, iKin_ptr, src) ;
        printSurf(cout, surfPhaseTP2, iKin2_ptr, src2) ;

        delete surfaceProb;
        surfaceProb = 0;
        delete iKin_ptr;
        iKin_ptr = 0;
        delete iKin2_ptr;
        iKin2_ptr = 0;
        delete gasTP;
        gasTP = 0;
        delete bulkPhaseTP;
        bulkPhaseTP = 0;
        delete surfPhaseTP;
        surfPhaseTP = 0;
        delete surfPhaseTP2;
        surfPhaseTP2 = 0;
        delete xc;
        appdelete();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }

    return 0;
}
/***********************************************************/
