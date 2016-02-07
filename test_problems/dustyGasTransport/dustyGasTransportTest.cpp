#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/transport/DustyGasTransport.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    try {
        int log_level = 0;

        unique_ptr<ThermoPhase> g(newPhase("h2o2.xml"));
        unique_ptr<Transport> tran(newTransportMgr("DustyGas", g.get(), log_level));
        DustyGasTransport* tranDusty = dynamic_cast<DustyGasTransport*>(tran.get());

        size_t nsp = g->nSpecies();
        vector_fp multiD(nsp*nsp);

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

        vector_fp state1;
        g->saveState(state1);
        g->setState_TP(T, 1.2 * OneAtm);
        vector_fp state2;
        g->saveState(state2);
        double delta = 0.001;
        vector_fp fluxes;
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
