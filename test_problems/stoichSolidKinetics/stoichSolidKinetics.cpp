#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include <iomanip>
#include <sstream>

using namespace Cantera;
using namespace std;

std::stringstream ss;

void printValue(const std::string& label, double value)
{
    ss.str("");
    ss.clear();
    ss << std::setw(13) << value;
    std::cout << label << ss.str() << std::endl;
}

void printRates(InterfaceKinetics& iKin)
{
    vector_fp work(iKin.nReactions(), 0.0);
    iKin.getNetRatesOfProgress(&work[0]);
    printValue("ROP_net       = ", work[0]);

    iKin.getFwdRatesOfProgress(&work[0]);
    printValue("ROP_forward   = ", work[0]);

    iKin.getRevRatesOfProgress(&work[0]);
    printValue("ROP_reverse   = ", work[0]);

    iKin.getFwdRateConstants(&work[0]);
    printValue("    kfwd      = ", work[0]);

    iKin.getRevRateConstants(&work[0]);
    printValue("    krev      = ", work[0]);
}

void testProblem()
{
    ss << std::scientific << std::setprecision(3) << std::uppercase;
    XML_Node* xc = get_XML_File("ReactionSurf.xml");
    XML_Node* xg = xc->findNameID("phase", "reaction_surface");
    if (!xg) {
        throw CanteraError("testProblem", "couldn't find file");
    }
    unique_ptr<ThermoPhase> surfTP(newPhase(*xg));
    unique_ptr<ThermoPhase> gasTP(newPhase("gas.xml"));

    unique_ptr<ThermoPhase> cao_s(newPhase("solidPhases.xml", "CaO(S)"));
    unique_ptr<ThermoPhase> caco3_s(newPhase("solidPhases.xml", "CaCO3(S)"));
    unique_ptr<ThermoPhase> c_s(newPhase("solidPhases.xml", "C(S)"));
    unique_ptr<ThermoPhase> fe3o4_s(newPhase("solidPhases.xml", "Fe3O4(S)"));
    unique_ptr<ThermoPhase> feo_s(newPhase("solidPhases.xml", "FeO(S)"));
    unique_ptr<ThermoPhase> fe_s(newPhase("solidPhases.xml", "Fe(S)"));

    vector<ThermoPhase*> phaseList {
        gasTP.get(), cao_s.get(), caco3_s.get(), c_s.get(), fe3o4_s.get(),
        feo_s.get(), fe_s.get(), surfTP.get()
    };
    InterfaceKinetics iKin;
    importKinetics(*xg, phaseList, &iKin);

    vector_fp mll(gasTP->nSpecies(), 0.0);
    size_t igco2 = gasTP->speciesIndex("CO2");
    size_t igo2  = gasTP->speciesIndex("O2");
    size_t ign2  = gasTP->speciesIndex("N2");

    mll[igco2] = 0.2;
    mll[igo2] = 0.1;
    mll[ign2] = 0.7;

    // Set the bath gas of 1000 K and 1 atm
    double Temp = 1000.;
    gasTP->setState_TPX(Temp, OneAtm, &mll[0]);
    cao_s->setState_TP(Temp, OneAtm);
    caco3_s->setState_TP(Temp, OneAtm);
    c_s->setState_TP(Temp, OneAtm);
    fe3o4_s->setState_TP(Temp, OneAtm);
    surfTP->setState_TP(Temp, OneAtm);

    vector_fp work(gasTP->nSpecies(), 0.0);

    cout << "*** StoichSolidKinetics Test ***" << endl;

    cout << "Tests for the proper behavior of heterogeneous reactions\n"
         << "when phases may or may not exist:\n"
         << "        CaCO3(s) =   CO2(g) +  CaO(s)  \n" << endl;

    for (int ktrials = 0; ktrials < 2; ktrials++) {
        iKin.getDeltaSSGibbs(&work[0]);
        printValue("   deltaGSS      = ", work[0]);

        iKin.getDeltaGibbs(&work[0]);
        printValue("   deltaG        = ", work[0]);

        gasTP->getChemPotentials(&work[0]);
        double mu_CO2 = work[igco2];
        printValue("   mu_CO2(g)     = ", mu_CO2);

        cao_s->getGibbs_RT(&work[0]);
        double mu_cao = work[0] * GasConstant * Temp;
        printValue("   mu_cao(s)     = ", mu_cao);

        caco3_s->getChemPotentials(&work[0]);
        double mu_caco3  = work[0];
        printValue("   mu_caco3      = ", mu_caco3);

        double deltaG_calc =  mu_CO2 +  mu_cao - mu_caco3;
        printValue("   deltaG_calc   = ", deltaG_calc);

        gasTP->getActivities(&work[0]);
        double act_CO2 = work[igco2];
        printValue("   act_CO2       = ", act_CO2);

        cao_s->getActivities(&work[0]);
        printValue("   act_cao(s)    = ", work[0]);

        caco3_s->getActivities(&work[0]);
        printValue("   act_caco3(s)  = ", work[0]);

        cout << "*** Base problem assuming that all phases exist:" << endl;
        printRates(iKin);

        cout << "*** Setting CaO(S) phase to nonexistent:" << endl;
        size_t ip_cao = iKin.phaseIndex("CaO(S)");
        iKin.setPhaseExistence(ip_cao, false);
        iKin.setPhaseStability(ip_cao, true);
        printRates(iKin);

        cout << "*** Setting CaCO3(S) phase to nonexistent:" << endl;
        size_t ip_caco3 = iKin.phaseIndex("CaCO3(S)");
        iKin.setPhaseExistence(ip_caco3, false);
        iKin.setPhaseStability(ip_caco3, true);
        printRates(iKin);

        cout << "*** OK Setting CaO(S) phase to existent, CaCO3 nonexistent:" << endl;
        iKin.setPhaseExistence(ip_cao, true);
        printRates(iKin);

        cout << "*** Setting Gas phase to nonexistent, CaCO3 nonexistent:" << endl;
        size_t ip_gas = iKin.phaseIndex("air");
        iKin.setPhaseExistence(ip_gas, false);
        iKin.setPhaseStability(ip_gas, true);
        printRates(iKin);

        cout << "*** Setting to all phases existing:" << endl;
        iKin.setPhaseExistence(ip_gas, true);
        iKin.setPhaseExistence(ip_cao, true);
        iKin.setPhaseExistence(ip_caco3, true);
        printRates(iKin);

        if (ktrials == 0) {
            cout << "*** Setting so that forward rate if faster:" << endl;
            mll[igco2] = 0.002;
            mll[igo2] =  0.1;
            mll[ign2] =  0.898;
            gasTP->setState_TPX(Temp, OneAtm, &mll[0]);
        }
    }
}

int main(int argc, char** argv)
{
    try {
        testProblem();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 0;
    }
}
