classdef ctTestFlowReactor2 < matlab.unittest.TestCase

    properties
        gas
        reactor
        net
        surf
        rsurf
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

    methods

        function makeReactors(self)
            self.reactor = FlowReactor(self.gas);
            % self.reactor.area = 1e-4;
            % self.reactor.surfaceAreaToVolumeRatio = 5000;
            self.reactor.massFlowRate = 0.02;
            self.rsurf = ReactorSurface(self.surf, self.reactor);
            self.net = ReactorNet(self.reactor);
        end

        function getPhases(self)
            self.surf = Interface('SiF4_NH3_mec.yaml', 'SI3N4');
            self.gas = self.surf.adjacent('gas');
        end

    end

    methods (Test)

        function testAdvanceReverse(self)
            self.getPhases();
            self.gas.TPX = {1500, 4000, 'NH3:1.0, SiF4:0.4'};
            self.makeReactors();

            self.net.advance(0.1);

            try
                self.net.advance(0.09);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message,  'backwards in time');
            end
        end

        function testNoMassFlowRate(self)
            self.getPhases();
            self.makeReactors();

            self.assumeFail('Skipped until more Reactor methods are implemented');

            try
                self.net.initialize();
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message,  'mass flow rate');
            end
        end

        function testMixedReactorTypes(self)
            self.getPhases();
            r1 = FlowReactor(self.gas);
            r2 = IdealGasReactor(self.gas);

            try
                self.net = ReactorNet({r1, r2});
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message,  'Cannot mix Reactor types');
            end

            delete(r1);
            delete(r2);
        end

        function testUnrecoverableIntegratorError(self)
            self.assumeFail('Skipped until addReactions is implemented');
        end

        function testIntegratorErrorAdvance(self)
            self.assumeFail('Skipped until addReactions is implemented');
        end

        function testRecoverableIntegratorErrors(self)
            self.assumeFail('Skipped until addReactions is implemented');
        end

        function testMaxSteps(self)
            self.getPhases();
            self.gas.TPX = {1500, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.TP = self.gas.TP;
            self.makeReactors();

            self.assumeFail('Skipped until maxSteps is implemented');

            % self.net.maxSteps = 13;
            % self.verifyEqual(self.net.maxSteps, 13);

            % try
            %     self.net.advance(0.1);
            % catch ME
            %     self.verifySubstring(ME.identifier, 'Cantera:ctError');
            %     self.verifySubstring(ME.message,  'Maximum number of timesteps');
            % end
        end

        function testMaxTimeStep(self)
            self.getPhases();
            self.gas.TPX = {1500, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.TP = self.gas.TP;
            self.makeReactors();

            self.net.maxTimeStep = 0.002;
            self.net.advance(0.01);
            self.net.step();

            self.assumeFail('Skipped until ReactorNet.distance is implemented');

            % x1 = self.net.distance;
            % x2 = self.net.step();
            % dx_limit = 0.1 * (x2 - x1);

            % Setting a step size limit to take one additional step before it's
            % fully enforced.

            % self.net.step();

            % for i = 1:20
            %     tPrev = self.net.distance;
            %     tNow = self.net.step();
            %     self.verifyEqual(tNow - tPrev, 1.0001 * dx_limit, ...
            %                      'RelTol', self.rtol);
            % end
        end

        function testIndependentVariable(self)
            self.getPhases();
            self.makeReactors();

            try
                t = self.net.time;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message,  'independent variable');
            end
        end

        function testIterationLimits(self)
            self.getPhases();
            self.gas.TPX = {1700, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.TP = self.gas.TP;
            self.makeReactors();

            self.net.advance(0.1);

            self.assumeFail('Skipped until more ReactorNet methods are implemented');

            % try
            %     t = self.net.advance(0.2);
            % catch ME
            %     self.verifySubstring(ME.identifier, 'Cantera:ctError');
            %     self.verifySubstring(ME.message,  'corrector convergence');
            % end
        end

        function testSolverOrder(self)
            self.getPhases();
            self.gas.TPX = {1700, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.TP = self.gas.TP;
            self.makeReactors();

            self.assumeFail('Skipped until more ReactorNet methods are implemented');

            % try
            %     t = self.net.maxOrder = -1;
            % catch ME
            %     self.verifySubstring(ME.identifier, 'Cantera:ctError');
            %     self.verifySubstring(ME.message,  'IDA_ILL_INPUT');
            % end
        end

        function testReinitialization(self)
            self.getPhases();
            self.gas.TPX = {1700, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.TP = self.gas.TP;
            self.makeReactors();
            self.reactor.massFlowRate = 0.01;

            self.assumeFail('Skipped until syncState is implemented');
        end

        function testInitialConditionTolerances(self)
            self.getPhases();
            self.gas.TPX = {1700, 4000, 'NH3:1.0, SiF4:0.4'};
            self.surf.coverages = {ones(1, self.surf.nSpecies)};
            self.surf.TP = self.gas.TP;
            self.makeReactors();

            self.assumeFail('Skipped until initialize is implemented');
        end

    end

end
