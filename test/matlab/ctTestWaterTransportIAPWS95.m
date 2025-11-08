classdef ctTestWaterTransportIAPWS95 < ctTestCase
    properties
        phase
    end

    methods (TestMethodSetup)

        function createPhase(self)
            self.phase = ct.Solution('../data/thermo-models.yaml', 'liquid-water');
        end

    end

    properties (TestParameter)
        TSet1 = {400, 400, 620, 620};
        PSet1 = {1e6, 8e6, 1.6e7, 2.8e7};
        muSet1 = {2.1880e-4, 2.2061e-4, 6.7489e-5, 7.5684e-5};
        rtolSet1 = {2e-4, 2e-4, 1e-4, 1e-4};

        TSet2 = {660, 660, 660};
        PSet2 = {2.2e7, 2.54e7, 2.8e7};
        muSet2 = {2.7129e-5, 3.8212e-5, 5.3159e-5};
        rtolSet2 = {1e-4, 1e-4, 1e-4};

        TSet3 = {400, 400, 620, 620};
        PSet3 = {1e6, 8e6, 1.6e7, 2.8e7};
        kSet3 = {0.68410, 0.68836, 0.45458, 0.49705};
        rtolSet3 = {1e-4, 1e-4, 1e-4, 1e-4};

        TSet4 = {660, 660, 660};
        PSet4 = {2.2e7, 2.54e7, 2.8e7};
        kSet4 = {0.14872, 0.35484, 0.38479};
        rtolSet4 = {1e-4, 1e-4, 1e-4};
    end

    methods (Test, ParameterCombination = 'sequential')

        function testViscosityLiquid(self, TSet1, PSet1, muSet1, rtolSet1)
            self.phase.TP = {TSet1, PSet1};
            self.verifyEqual(self.phase.viscosity, muSet1, 'RelTol', rtolSet1);
        end

        function testViscositySuperCritical(self, TSet2, PSet2, muSet2, rtolSet2)
            self.phase.TP = {TSet2, PSet2};
            self.verifyEqual(self.phase.viscosity, muSet2, 'RelTol', rtolSet2);
        end


        function testThermalConductivityLiquid(self, TSet3, PSet3, kSet3, rtolSet3)
            self.phase.TP = {TSet3, PSet3};
            self.verifyEqual(self.phase.thermalConductivity, kSet3, 'RelTol', rtolSet3);
        end

        function testThermalConductivitySuperCrit(self, TSet4, PSet4, kSet4, rtolSet4)
            self.phase.TP = {TSet4, PSet4};
            self.verifyEqual(self.phase.thermalConductivity, kSet4, 'RelTol', rtolSet4);
        end

    end

end
