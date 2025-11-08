classdef ctTestWaterTransport < ctTestCase

    properties
        phase
    end

    methods (TestMethodSetup)

        function createPhase(self)
            self.phase = ct.Water();
        end

    end

    properties (TestParameter)
        TSet1 = {400, 400, 620, 620};
        PSet1 = {1e6, 8e6, 1.6e7, 2.8e7};
        muSet1 = {2.1880e-4, 2.2061e-4, 6.7489e-5, 7.5684e-5};
        rtolSet1 = {1e-3, 1e-3, 2e-3, 2e-3};

        TSet2 = {600, 620, 620};
        PSet2 = {1e6, 5e6, 1.5e7};
        muSet2 = {2.1329e-5, 2.1983e-5, 2.2858e-5};
        rtolSet2 = {1e-3, 1e-3, 2e-3};

        TSet3 = {660, 660, 660};
        PSet3 = {2.2e7, 2.54e7, 2.8e7};
        muSet3 = {2.7129e-5, 3.8212e-5, 5.3159e-5};
        rtolSet3 = {2e-3, 1e-2, 2e-2};

        TSet4 = {647.43, 647.43, 648.23, 648.23};
        DSet4 = {280.34, 318.89, 301.34, 330.59};
        muSet4 = {3.7254e-5, 4.2286e-5, 3.9136e-5, 4.2102e-5};
        rtolSet4 = {6e-3, 6e-3, 6e-3, 6e-3};

        TSet5 = {400, 400, 620, 620};
        PSet5 = {1e6, 8e6, 1.6e7, 2.8e7};
        kSet5 = {0.68410, 0.68836, 0.45458, 0.49705};
        rtolSet5 = {1e-3, 1e-3, 2e-3, 2e-3};

        TSet6 = {600, 620, 620};
        PSet6 = {1e6, 5e6, 1.5e7};
        kSet6 = {0.047636, 0.055781, 0.10524};
        rtolSet6 = {1e-3, 1e-3, 2e-3};

        TSet7 = {660, 660, 660};
        PSet7 = {2.2e7, 2.54e7, 2.8e7};
        kSet7 = {0.14872, 0.35484, 0.38479};
        rtolSet7 = {1e-2, 2e-2, 1e-2};
    end

    methods (Test, ParameterCombination = 'sequential')

        function testViscosityLiquid(self, TSet1, PSet1, muSet1, rtolSet1)
            self.phase.TP = {TSet1, PSet1};
            self.verifyEqual(self.phase.viscosity, muSet1, 'RelTol', rtolSet1);
        end

        function testViscosityVapor(self, TSet2, PSet2, muSet2, rtolSet2)
            self.phase.TP = {TSet2, PSet2};
            self.verifyEqual(self.phase.viscosity, muSet2, 'RelTol', rtolSet2);
        end

        function testViscositySuperCritical(self, TSet3, PSet3, muSet3, rtolSet3)
            self.phase.TP = {TSet3, PSet3};
            self.verifyEqual(self.phase.viscosity, muSet3, 'RelTol', rtolSet3);
        end

        function testViscosityNearCritical(self, TSet4, DSet4, muSet4, rtolSet4)
            self.phase.basis = 'mass';
            self.phase.TD = {TSet4, DSet4};
            self.verifyEqual(self.phase.viscosity, muSet4, 'RelTol', rtolSet4);
        end

        function testThermalConductivityLiquid(self, TSet5, PSet5, kSet5, rtolSet5)
            self.phase.TP = {TSet5, PSet5};
            self.verifyEqual(self.phase.thermalConductivity, kSet5, 'RelTol', rtolSet5);
        end

        function testThermalConductivityVapor(self, TSet6, PSet6, kSet6, rtolSet6)
            self.phase.TP = {TSet6, PSet6};
            self.verifyEqual(self.phase.thermalConductivity, kSet6, 'RelTol', rtolSet6);
        end

        function testThermalConductivitySuperCrit(self, TSet7, PSet7, kSet7, rtolSet7)
            self.phase.TP = {TSet7, PSet7};
            self.verifyEqual(self.phase.thermalConductivity, kSet7, 'RelTol', rtolSet7);
        end

    end

end
