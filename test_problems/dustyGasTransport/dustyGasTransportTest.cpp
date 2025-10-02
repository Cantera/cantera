#include "cantera/thermo/ThermoFactory.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/transport/DustyGasTransport.h"
#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    try {
        auto g = newThermo("h2o2.yaml");
        auto tran = newTransport(g, "DustyGas");
        auto tranDusty = std::dynamic_pointer_cast<DustyGasTransport>(tran);

        size_t nsp = g->nSpecies();
        vector<double> multiD(nsp*nsp);

        double T = 500;
        g->setState_TPX(T, OneAtm,
                        "OH:1, H:2, O2:3, O:1.0E-8, H2:1.0E-8, H2O:1.0E-8, H2O2:1.0E-8, HO2:1.0E-8, AR:1.0E-8");

        tranDusty->setPorosity(0.2);
        tranDusty->setTortuosity(4.0);
        tranDusty->setMeanPoreRadius(1.5E-7);
        tranDusty->setMeanParticleDiameter(1.5E-6);

        tranDusty->getMultiDiffCoeffs(nsp, multiD.data());
        printf("MultiDiffusion coefficients: \n");
        for (size_t i = 0; i < nsp; i++) {
            for (size_t j = 0; j < nsp; j++) {
                printf(" %15.8E,", multiD[nsp*j + i]);
            }
            printf("\n");
        }

        vector<double> state1;
        g->saveState(state1);
        g->setState_TP(T, 1.2 * OneAtm);
        vector<double> state2;
        g->saveState(state2);
        double delta = 0.001;
        vector<double> fluxes;
        fluxes.resize(nsp);

        tranDusty->getMolarFluxes(&state1[0], &state1[0], delta, &fluxes[0]);
        for (size_t i = 0; i < nsp; i++) {
            printf(" flux[%d] = %13.8E\n", int(i), fluxes[i]);
        }
        tranDusty->getMolarFluxes(&state1[0], &state2[0], delta, &fluxes[0]);
        for (size_t i = 0; i < nsp; i++) {
            printf(" flux[%d] = %13.8E\n", int(i), fluxes[i]);
        }

        Cantera::appdelete();
        return 0;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
