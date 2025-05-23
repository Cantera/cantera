classdef ctTestReactorSensitivities < matlab.unittest.TestCase

    properties
        gas1
        gas2
        r1
        r2
        interface
        net
        surf
    end

    properties (SetAccess = protected)
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

    methods (TestMethodTeardown)


        function deleteObjects(self)
            props = properties(self);
            for i = 1:length(props)
                prop = self.(props{i});
                if isa(prop, 'handle')
                    delete(prop)
                end
            end
        end

    end

    methods (Test)

        function testSensitivities1(self)
            self.net = ReactorNet();
            self.gas1 = Solution('gri30.yaml', '', 'none');
            self.gas1.TPX = {1300, 20 * OneAtm, 'CO:1.0, H2:0.1, CH4:0.1, H2O:0.5'};
            self.r1 = IdealGasReactor(self.gas1);
            self.net.addReactor(self.r1);

            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
            % self.verifyEqual(self.net.nSensitivityParams, 0);

            self.r1.addSensitivityReaction(41);
            self.r1.addSensitivityReaction(42);

            self.net.advance(0.1);

            self.verifyEqual(self.net.nSensitivityParams, 2);
            val = self.gas.nSpecies + self.r1.componentIndex(self.gas.speciesName(1));
            self.verifyEqual(self.net.nVars, val);

            S = self.net.sensitivities();
            val = size(S) == [self.net.nVars, self.net.nSensitivityParams];
            self.verifyTrue(all(val));
        end

        function testCoveragesRegression1(self)
            self.net = ReactorNet();
            self.interface = Interface('diamond.yaml', 'diamond_100');
            self.gas1 = self.interface.adjacent('gas');
            self.r1 = IdealGasReactor(self.gas1);
            self.net.addReactor(self.r1);
            self.net.setSensitivityTolerances(1e-8, 1e-10);

            self.gas2 = Solution('h2o2.yaml', '', 'none');
            self.gas2.TPX = {900, OneAtm, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'};
            self.r2 = IdealGasReactor(self.gas2);
            self.net.addReactor(self.r2);

            self.surf = ReactorSurface(self.interface, self.r1);
            self.surf.area = 1.5;

            self.assumeFail('Skipped until more methods are implemented for ReactorSurface');
        end

        function testParameterOrder1a(self)
            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
        end

        function testParameterOrder1b(self)
            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
        end

        function testParameterOrder2(self)
            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
        end

        function testParameterOrder3(self)
            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
        end

        function testIgnitionDelaySensitivity(self)
            self.assumeFail('Skipped until more methods are implemented for ReactorNet');
        end

    end

end
