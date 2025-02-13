classdef ctTestKineticsRepeatability < matlab.unittest.TestCase

    properties
        phase
        X0
        X1
    end

    properties (SetAccess = immutable)
        T0 = 1200;
        T1 = 1300;
        rho0 = 2.4;
        rho1 = 3.1;
        P0 = 1.4e+05;
        P1 = 3.7e+06;
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
            copyfile('../data/pdep-test.yaml', ...
                     './pdep-test.yaml');
            copyfile('../data/sticking_coeff_check.yaml', ...
                     './sticking_coeff_check.yaml');
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            delete('./pdep-test.yaml');
            delete('./sticking_coeff_check.yaml');
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods

        function setup_phase(self, mech)
            clear self.phase
            self.phase = Solution(mech);
            self.X0 = 1 + sin(1:self.phase.nSpecies);
            self.X1 = 1 + sin(2:self.phase.nSpecies + 1);
        end

        function checkRatesX(self, mech)
            self.setup_phase(mech);

            self.phase.TDX = {self.T0, self.rho0, self.X0};
            w1 = self.phase.netProdRates;

            self.phase.TDX = {self.T1, self.rho1, self.X1};
            w2 = self.phase.netProdRates;

            self.phase.TDX = {self.T0, self.rho0, self.X1};
            w3 = self.phase.netProdRates;

            self.phase.TDX = {self.T0, self.rho0, self.X0};
            w4 = self.phase.netProdRates;

            self.verifyEqual(w1, w4, 'RelTol', self.rtol);
        end

        function checkRatesT1(self, mech)
            self.setup_phase(mech);

            self.phase.TDX = {self.T0, self.rho0, self.X0};
            w1 = self.phase.netProdRates;

            self.phase.TDX = {self.T1, self.rho1, self.X1};
            w2 = self.phase.netProdRates;

            self.phase.TDX = {self.T1, self.rho0, self.X0};
            w3 = self.phase.netProdRates;

            self.phase.TDX = {self.T0, self.rho0, self.X0};
            w4 = self.phase.netProdRates;

            self.verifyEqual(w1, w4, 'RelTol', self.rtol);
        end

        function checkRatesT2(self, mech)
            self.setup_phase(mech);

            self.phase.TPX = {self.T0, self.P0, self.X0};
            w1 = self.phase.netProdRates;

            self.phase.TPX = {self.T1, self.P1, self.X1};
            w2 = self.phase.netProdRates;

            self.phase.TPX = {self.T1, self.P0, self.X0};
            w3 = self.phase.netProdRates;

            self.phase.TPX = {self.T0, self.P0, self.X0};
            w4 = self.phase.netProdRates;

            self.verifyEqual(w1, w4, 'RelTol', self.rtol);
        end

        function checkRatesP(self, mech)
            self.setup_phase(mech);

            self.phase.TPX = {self.T0, self.P0, self.X0};
            w1 = self.phase.netProdRates;

            self.phase.TPX = {self.T1, self.P1, self.X1};
            w2 = self.phase.netProdRates;

            self.phase.TPX = {self.T0, self.P1, self.X0};
            w3 = self.phase.netProdRates;

            self.phase.TPX = {self.T0, self.P0, self.X0};
            w4 = self.phase.netProdRates;

            self.verifyEqual(w1, w4, 'RelTol', self.rtol);
        end

    end

    methods (Test)

        function testGRI30X(self)
            self.checkRatesX('gri30.yaml');
        end

        function testGRI30T(self)
            self.checkRatesT1('gri30.yaml');
            self.checkRatesT2('gri30.yaml');
        end

        function testGRI30P(self)
            self.checkRatesP('gri30.yaml');
        end

        function testH2O2X(self)
            self.checkRatesX('h2o2.yaml');
        end

        function testH2O2T(self)
            self.checkRatesT1('h2o2.yaml');
            self.checkRatesT2('h2o2.yaml');
        end

        function testH2O2P(self)
            self.checkRatesP('h2o2.yaml');
        end

        function testPdepX(self)
            self.checkRatesX('pdep-test.yaml');
        end

        function testPdepT(self)
            self.checkRatesT1('pdep-test.yaml');
            self.checkRatesT2('pdep-test.yaml');
        end

        function testPdepP(self)
            self.checkRatesP('pdep-test.yaml');
        end

    end

end