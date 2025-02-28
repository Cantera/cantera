classdef ctTestPureFluid < matlab.unittest.TestCase

    properties
        phase
    end

    properties (SetAccess = immutable)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodSetup)

        function createPhase(self)
            self.phase = Water;
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods (Test)

        function testTemperatureLimits(self)
            CO2 = CarbonDioxide;
            self.verifyEqual(CO2.minTemp, 216.54, 'RelTol', self.rtol);
            self.verifyEqual(CO2.maxTemp, 1500.0, 'RelTol', self.rtol);

            clear CO2
        end

        function testCriticalProperties(self)
            self.verifyEqual(self.phase.critPressure, 22.089e6, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.critTemperature, 647.286, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.critPressure, 317.0, 'RelTol', self.rtol);
        end

        function testSetState(self)
            self.phase.PQ = {101325, 0.5};
            self.verifyEqual(self.phase.P, 101325, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.Q, 0.5, 'RelTol', self.rtol);

            self.phase.TQ = {500, 0.8};
            self.verifyEqual(self.phase.T, 500, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.Q, 0.8, 'RelTol', self.rtol);
        end

        function testSubstanceSet(self)
            self.phase.TV = {400, 1.45};
            self.verifyEqual(self.phase.T, 400, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.V, 1.45, 'RelTol', self.rtol);
            try
                self.phase.TV = {300, -1};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Negative specific volume');
            end

            self.phase.PV = {101325, 1.45};
            self.verifyEqual(self.phase.P, 101325, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.V, 1.45, 'RelTol', self.rtol);

            self.phase.UP = {-1.45e7, 101325};
            self.verifyEqual(self.phase.U, -1.45e7, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.P, 101325, 'RelTol', self.rtol);

            self.phase.VH = {1.45, -1.45e7};
            self.verifyEqual(self.phase.V, 1.45, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.H, -1.45e7, 'RelTol', self.rtol);

            self.phase.TH = {400, -1.45e7};
            self.verifyEqual(self.phase.T, 400, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.H, -1.45e7, 'RelTol', self.rtol);

            self.phase.SH = {5000, -1.45e7};
            self.verifyEqual(self.phase.S, 5000, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.H, -1.45e7, 'RelTol', self.rtol);

            self.phase.ST = {5000, 400};
            self.verifyEqual(self.phase.S, 5000, 'RelTol', self.rtol);
            self.verifyEqual(self.phase.T, 400, 'RelTol', self.rtol);
        end

        function testSetMinMax(self)
            self.phase.TP = {self.phase.minTemp, 101325};
            self.verifyEqual(self.phase.T, self.phase.minTemp, ...
                             'RelTol', self.rtol);

            self.phase.TP = {self.phase.maxTemp, 101325};
            self.verifyEqual(self.phase.T, self.phase.maxTemp, ...
                            'RelTol', self.rtol);
        end

    end

end
