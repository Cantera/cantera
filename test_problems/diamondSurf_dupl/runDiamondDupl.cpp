/**
 *  @file runDiamond.cpp
 *
 */

//  Example
//
// Note that this example needs updating. It works fine, but is
// written in a way that is less than transparent or
// user-friendly. This could be rewritten using class Interface to
// make things simpler.

#include "cantera/kinetics.h"

using namespace std;
using namespace Cantera;

void printDbl(double val)
{
    if (fabs(val) < 5.0E-17) {
        cout << " nil";
    } else {
        cout << val;
    }
}

int main(int argc, char** argv)
{
    int i, k;
    string infile = "diamond.xml";

    try {
        XML_Node* xc = new XML_Node();
        string path = findInputFile(infile);
        ctml::get_CTML_Tree(xc, path);
        cout.precision(3);

        XML_Node* const xg = xc->findNameID("phase", "gas");
        ThermoPhase* gasTP = newPhase(*xg);
        int nsp = gasTP->nSpecies();
        cout << "Number of species = " << nsp << endl;

        XML_Node* const xd = xc->findNameID("phase", "diamond");
        ThermoPhase* diamondTP = newPhase(*xd);
        int nsp_diamond = diamondTP->nSpecies();
        cout << "Number of species in diamond = " << nsp_diamond << endl;


        XML_Node* const xs = xc->findNameID("phase", "diamond_100");
        ThermoPhase* diamond100TP = newPhase(*xs);
        //SurfPhase *diamond100TP = new SurfPhase(*xs);
        int nsp_d100 = diamond100TP->nSpecies();
        cout << "Number of species in diamond_100 = " << nsp_d100 << endl;

        vector<ThermoPhase*> phaseList;
        phaseList.push_back(gasTP);
        phaseList.push_back(diamondTP);
        phaseList.push_back(diamond100TP);
        InterfaceKinetics* iKin_ptr = new InterfaceKinetics();
        importKinetics(*xs, phaseList, iKin_ptr);
        int nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        double x[20];
        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 0.0010;
        x[1] = 0.9888;
        x[2] = 0.0002;
        x[3] = 0.0100;
        double p = 20.0*OneAtm/760.0;

        gasTP->setState_TPX(1200., p, x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        int i0 = diamond100TP->speciesIndex("c6H*");
        x[i0] = 0.1;
        int i1 = diamond100TP->speciesIndex("c6HH");
        x[i1] = 0.9;
        diamond100TP->setState_TX(1200., x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 1.0;
        diamondTP->setState_TPX(1200., p, x);

        iKin_ptr->advanceCoverages(100.);

        // Throw some asserts in here to test that they compile
        AssertTrace(p == p);
        AssertThrow(p == p, "main");
        AssertThrowMsg(i == 20, "main", "are you kidding");

        double src[20];
        for (i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        iKin_ptr->getNetProductionRates(src);
        double sum = 0.0;
        double naH = 0.0;
        for (k = 0; k < 13; k++) {
            if (k < 4) {
                naH = gasTP->nAtoms(k, 0);
            } else if (k == 4) {
                naH = 0;
            } else if (k > 4) {
                int itp = k - 5;
                naH = diamond100TP->nAtoms(itp, 0);
            }
            cout << k << "  " << naH << "  " ;
            printDbl(src[k]);
            cout << endl;
            sum += naH * src[k];

        }

        cout << "sum = ";
        printDbl(sum);
        cout << endl;
        double mwd = diamondTP->molecularWeight(0);
        double dens = diamondTP->density();
        double gr = src[4] * mwd / dens;
        gr *= 1.0E6 * 3600.;
        cout << "growth rate = " << gr << " microns per hour" << endl;


        diamond100TP->getMoleFractions(x);
        cout << "Coverages:" << endl;
        for (k = 0; k < 8; k++) {
            cout << k << "   " << diamond100TP->speciesName(k)
                 << "   "
                 << x[k] << endl;
        }


        /*********************************************************************************/
        /*
         *  OK NOW DUPLICATE EVERYTHING AND RECALCULATE
         */
        ThermoPhase* gasTP_dupl         = gasTP->duplMyselfAsThermoPhase();
        ThermoPhase* diamondTP_dupl     = diamondTP->duplMyselfAsThermoPhase();
        ThermoPhase* diamond100TP_dupl  = diamond100TP->duplMyselfAsThermoPhase();


        vector<ThermoPhase*> phaseList_dupl;
        phaseList_dupl.push_back(gasTP_dupl);
        phaseList_dupl.push_back(diamondTP_dupl);
        phaseList_dupl.push_back(diamond100TP_dupl);
        InterfaceKinetics* iKin_ptr_dupl = new InterfaceKinetics();
        importKinetics(*xs, phaseList_dupl, iKin_ptr_dupl);
        int nr_dupl = iKin_ptr_dupl->nReactions();
        cout << "Number of reactions = " << nr_dupl << endl;


        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 0.0010;
        x[1] = 0.9888;
        x[2] = 0.0002;
        x[3] = 0.0100;
        p = 20.0*OneAtm/760.0;

        gasTP_dupl->setState_TPX(1200., p, x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        i0 = diamond100TP_dupl->speciesIndex("c6H*");
        x[i0] = 0.1;
        i1 = diamond100TP_dupl->speciesIndex("c6HH");
        x[i1] = 0.9;
        diamond100TP_dupl->setState_TX(1200., x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 1.0;
        diamondTP_dupl->setState_TPX(1200., p, x);

        iKin_ptr_dupl->advanceCoverages(100.);

        // Throw some asserts in here to test that they compile
        AssertTrace(p == p);
        AssertThrow(p == p, "main");
        AssertThrowMsg(i == 20, "main", "are you kidding");


        for (i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        iKin_ptr_dupl->getNetProductionRates(src);
        sum = 0.0;
        naH = 0.0;
        for (k = 0; k < 13; k++) {
            if (k < 4) {
                naH = gasTP_dupl->nAtoms(k, 0);
            } else if (k == 4) {
                naH = 0;
            } else if (k > 4) {
                int itp = k - 5;
                naH = diamond100TP_dupl->nAtoms(itp, 0);
            }
            cout << k << "  " << naH << "  " ;
            printDbl(src[k]);
            cout << endl;
            sum += naH * src[k];

        }

        cout << "sum = ";
        printDbl(sum);
        cout << endl;
        mwd = diamondTP_dupl->molecularWeight(0);
        dens = diamondTP_dupl->density();
        gr = src[4] * mwd / dens;
        gr *= 1.0E6 * 3600.;
        cout << "growth rate = " << gr << " microns per hour" << endl;


        diamond100TP_dupl->getMoleFractions(x);
        cout << "Coverages:" << endl;
        for (k = 0; k < 8; k++) {
            cout << k << "   " << diamond100TP_dupl->speciesName(k)
                 << "   "
                 << x[k] << endl;
        }



    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }

    return 0;
}
/***********************************************************/
