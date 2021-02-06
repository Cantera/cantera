#include "gtest/gtest.h"

#include "cantera/transport/TransportData.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

#include "../thermo/thermo_data.h"

using namespace Cantera;

class TransportFromScratch : public testing::Test
{
public:
    TransportFromScratch()
        : sH2(new Species("H2", parseCompString("H:2")))
        , sO2(new Species("O2", parseCompString("O:2")))
        , sH2O(new Species("H2O", parseCompString("H:2 O:1")))
        , tH2(new GasTransportData())
        , tO2(new GasTransportData())
        , tH2O(new GasTransportData())
    {
        sH2->thermo.reset(new NasaPoly2(200, 3500, 101325, h2_nasa_coeffs));
        sO2->thermo.reset(new NasaPoly2(200, 3500, 101325, o2_nasa_coeffs));
        sH2O->thermo.reset(new NasaPoly2(200, 3500, 101325, h2o_nasa_coeffs));

        tH2->setCustomaryUnits("linear", 2.92, 38.0, 0.0, 0.79, 280.0);
        tO2->setCustomaryUnits("linear", 3.458, 107.40, 0.0, 1.60, 3.80);
        tH2O->setCustomaryUnits("nonlinear", 2.605, 572.4, 1.844, 0.0, 4.00);

        sH2->transport = tH2;
        sO2->transport = tO2;
        sH2O->transport = tH2O;

        std::string phase_def = "ideal_gas(name='test', elements='O H',"
            "species='gri30: H2 O2 H2O')";

        XML_Node* fxml = get_XML_from_string(phase_def);
        ref.reset(newPhase(*fxml->findByName("phase")));
        test.reset(new IdealGasPhase());

        test->addElement("O");
        test->addElement("H");
        test->addSpecies(sH2);
        test->addSpecies(sO2);
        test->addSpecies(sH2O);
        test->initThermo();

        ref->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        test->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    }

    shared_ptr<Species> sH2, sO2, sH2O;
    shared_ptr<GasTransportData> tH2, tO2, tH2O;
    shared_ptr<ThermoPhase> ref;
    shared_ptr<ThermoPhase> test;
};

TEST_F(TransportFromScratch, binaryDiffCoeffs)
{
    Transport* trRef = newTransportMgr("Mix", ref.get());
    MixTransport trTest;
    trTest.init(test.get());

    size_t K = ref->nSpecies();
    Array2D bdiffRef(3,3);
    Array2D bdiffTest(3,3);
    ref->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    test->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    trRef->getBinaryDiffCoeffs(K, &bdiffRef(0,0));
    trTest.getBinaryDiffCoeffs(K, &bdiffTest(0,0));

    for (size_t i=0; i < K; i++) {
        for (size_t j=0; j < K; j++) {
            EXPECT_DOUBLE_EQ(bdiffRef(i,j), bdiffTest(i,j)) << "i = " << i << ", j = " << j;
        }
    }
}

TEST_F(TransportFromScratch, mixDiffCoeffs)
{
    Transport* trRef = newTransportMgr("Mix", ref.get());
    MixTransport trTest;
    trTest.init(test.get());

    size_t K = ref->nSpecies();
    vector_fp Dref(3);
    vector_fp Dtest(3);
    ref->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    test->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    trRef->getMixDiffCoeffs(&Dref[0]);
    trTest.getMixDiffCoeffs(&Dtest[0]);

    for (size_t k=0; k < K; k++) {
        EXPECT_DOUBLE_EQ(Dref[k], Dtest[k]) << "k = " << k;
    }
}

TEST_F(TransportFromScratch, viscosity)
{
    Transport* trRef = newTransportMgr("Mix", ref.get());
    MixTransport trTest;
    trTest.init(test.get());

    for (int i = 0; i < 10; i++) {
        double T = 300 + 111*i;
        ref->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        test->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        EXPECT_DOUBLE_EQ(trRef->viscosity(), trTest.viscosity()) << "T = " << T;
    }
}

TEST_F(TransportFromScratch, thermalConductivityMix)
{
    Transport* trRef = newTransportMgr("Mix", ref.get());
    MixTransport trTest;
    trTest.init(test.get());

    for (int i = 0; i < 10; i++) {
        double T = 300 + 111*i;
        ref->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        test->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        EXPECT_DOUBLE_EQ(trRef->thermalConductivity(),
                         trTest.thermalConductivity()) << "T = " << T;
    }
}

TEST_F(TransportFromScratch, multiDiffCoeffs)
{
    Transport* trRef = newTransportMgr("Multi", ref.get());
    MultiTransport trTest;
    trTest.init(test.get());

    size_t K = ref->nSpecies();
    Array2D Dref(3,3);
    Array2D Dtest(3,3);
    ref->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    test->setState_TPX(400, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
    trRef->getMultiDiffCoeffs(K, &Dref(0,0));
    trTest.getMultiDiffCoeffs(K, &Dtest(0,0));

    for (size_t i=0; i < K; i++) {
        for (size_t j=0; j < K; j++) {
            EXPECT_DOUBLE_EQ(Dref(i,j), Dtest(i,j)) << "i = " << i << ", j = " << j;
        }
    }
}

TEST_F(TransportFromScratch, thermalConductivityMulti)
{
    Transport* trRef = newTransportMgr("Multi", ref.get());
    MultiTransport trTest;
    trTest.init(test.get());

    for (int i = 0; i < 10; i++) {
        double T = 300 + 111*i;
        ref->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        test->setState_TPX(T, 5e5, "H2:0.5, O2:0.3, H2O:0.2");
        EXPECT_DOUBLE_EQ(trRef->thermalConductivity(),
                         trTest.thermalConductivity()) << "T = " << T;
    }
}

int main(int argc, char** argv)
{
    printf("Running main() from transportFromScratch.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    appdelete();
    return result;
}
