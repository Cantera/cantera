classdef ctTestReactor < ctTestCase

    properties
        gas1
        gas2
        r1
        r2
        w
    end

    properties (SetAccess = protected)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods

        function makeReactors(self, arg)
            arguments
                self
                arg.nr (1,1) double {mustBeInteger} = 2
                arg.T1 (1,1) double {mustBeNumeric} = 300
                arg.P1 (1,1) double {mustBeNumeric} = 101325
                arg.X1 (1,:) char = 'O2:1.0';
                arg.T2 (1,1) double {mustBeNumeric} = 300
                arg.P2 (1,1) double {mustBeNumeric} = 101325
                arg.X2 (1,:) char = 'O2:1.0';
            end

            self.gas1 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {arg.T1, arg.P1, arg.X1};
            self.r1 = ct.zeroD.Reactor(self.gas1);

            if arg.nr == 2
                self.gas2 = ct.Solution('h2o2.yaml', '', 'none');
                self.gas2.TPX = {arg.T2, arg.P2, arg.X2};
                self.r2 = ct.zeroD.Reactor(self.gas2);
                self.r2.energyEnabled = true;
            end
        end

        function addWall(self, arg)
            arguments
                self
                arg.K (1,1) double {mustBeNumeric} = 0.0
                arg.U (1,1) double {mustBeNumeric} = 0.0
                arg.A (1,1) double {mustBeNumeric} = 1.0
            end

            self.w = ct.zeroD.Wall(self.r1, self.r2);
            self.w.area = arg.A;
            self.w.expansionRateCoeff = arg.K;
            self.w.heatTransferCoeff = arg.U;
        end

        function m = sensitivityReaction(self, gas)
            % 1-based index of the chain-branching reaction, to be used as a
            % sensitivity parameter.
            m = find(strcmp(gas.reactionEquations, 'H + O2 <=> O + OH'), 1);
            self.assertNotEmpty(m);
        end

        function net = makeSensitivityNet(self)
            % Single igniting reactor with one sensitivity parameter, advanced
            % past ignition so that the sensitivity coefficients are nonzero.
            self.gas1 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {1050, 5 * ct.OneAtm, 'H2:2.0, O2:1.0, AR:4.0'};
            self.r1 = ct.zeroD.IdealGasReactor(self.gas1);
            net = ct.zeroD.ReactorNet(self.r1);

            % A reactor has to be part of a network before a sensitivity
            % parameter can be added to it.
            self.r1.addSensitivityReaction(self.sensitivityReaction(self.gas1));
            net.setSensitivityTolerances(1e-6, 1e-6);
            net.advance(0.1);
        end

    end

    methods (Test)

        function testV(self)
            g = ct.Solution('h2o2.yaml', '', 'none');
            r = ct.zeroD.Reactor(g);
            self.verifyEqual(r.V, 1.0, 'AbsTol', self.atol);

            r.V = 9;
            self.verifyEqual(r.V, 9.0, 'AbsTol', self.atol);
        end

        function testTypes(self)
            self.makeReactors();
            self.verifyEqual(self.r1.type, 'Reactor');
        end

        function testIndependentVariable(self)
            self.makeReactors('nr', 1);
            net = ct.zeroD.ReactorNet(self.r1);
            self.verifyEqual(net.time, 0.0, 'AbsTol', self.atol);
        end

        function testDisjoint(self)
            self.makeReactors('T1', 300, 'P1', 101325, 'T2', 500, 'P2', 300000);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.advance(1.0);

            self.verifyEqual(self.gas1.T, 300, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.T, 500, 'AbsTol', self.atol);
            self.verifyEqual(self.gas1.P, 101325, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.P, 300000, 'AbsTol', self.atol);
        end

        function testTimeStepping(self)
            self.makeReactors();
            net = ct.zeroD.ReactorNet({self.r1, self.r2});

            tStart = 0.3;
            tEnd = 10.0;
            dtMax = 0.07;
            t = tStart;

            net.maxTimeStep = dtMax;
            net.time = tStart;
            self.verifyEqual(net.time, tStart);

            while t < tEnd
                tPrev = t;
                t = net.step;
                self.verifyLessThanOrEqual(t - tPrev, 1.0001 * dtMax);
                self.verifyEqual(t, net.time, 'AbsTol', self.atol);
            end
        end

        function testWall1(self)
            self.makeReactors('P1', 101325, 'P2', 300000);
            self.addWall('K', 0.1, 'A', 1.0);
            self.verifyEqual(self.w.type, 'Wall');
            self.w.name = 'name-of-wall';
            self.verifyEqual(self.w.name, 'name-of-wall');
        end

        function testWall2(self)
            self.makeReactors('nr', 1);
            res = ct.zeroD.Reservoir(self.gas1);
            w = ct.zeroD.Wall(self.r1, res);
            net = ct.zeroD.ReactorNet([self.r1]);
            self.verifyEqual(w.type, 'Wall');
        end

        function testWall3(self)
            self.makeReactors('nr', 1);
            res = ct.zeroD.Reservoir(self.gas1);
            w = ct.zeroD.Wall(res, self.r1);
            net = ct.zeroD.ReactorNet([self.r1]);
            self.verifyEqual(w.type, 'Wall');
        end

        function testEqualizePressure(self)
            self.makeReactors('P1', 101325, 'P2', 300000);
            self.addWall('K', 0.1, 'A', 1.0);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.advance(1.0);
            self.verifyEqual(net.time, 1.0, 'AbsTol', self.atol);
            self.verifyEqual(self.r1.phase.P, self.r2.phase.P, 'RelTol', self.rtol);
            self.verifyNotEqual(self.r1.T, self.r2.T);
        end

        function testHeatTransfer1(self)
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 500, 'A', 1.0);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});

            net.advance(10.0);
            self.verifyEqual(net.time, 10.0, 'RelTol', self.atol);
            self.verifyEqual(self.r1.T, self.r2.T, 'RelTol', 5e-7);
            self.verifyNotEqual(self.r1.P, self.r2.P);
        end

        function testHeatTransfer2(self)
            self.assumeFail('Skipped until heatTransferCoeff getter is implemented');
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 200, 'A', 1.0);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});

            net.advance(1.0);
            T1a = self.r1.T;
            T2a = self.r2.T;

            self.makeReactors('T1', 300, 'T2', 1000);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            self.r1.V = 0.25;
            self.r2.V = 0.25;
            self.addWall('U', 100, 'A', 0.5);

            Qdot1 = self.w.heatTransferCoeff * self.w.area * (self.r1.T - self.r2.T);
            self.verifyEqual(Qdot1, self.w.heatRate, 'RelTol', self.rtol);

            net.advance(1.0);
            Qdot2 = self.w.heatTransferCoeff * self.w.area * (self.r1.T - self.r2.T);
            self.verifyEqual(Qdot2, self.w.heatRate, 'RelTol', self.rtol);
            T1b = self.r1.T;
            T2b = self.r2.T;

            self.verifyEqual(T1a, T1b, 'RelTol', self.rtol);
            self.verifyEqual(T2a, T2b, 'RelTol', self.rtol);
        end

        function testTolerances(self)
            function n = integrate(atol, rtol)
                P0 = 10 * ct.OneAtm;
                T0 = 1100;
                X0 = 'H2:1.0, O2:0.5, AR:8.0';
                self.makeReactors('nr', 1, 'T1', T0, 'P1', P0, 'X1', X0);
                net = ct.zeroD.ReactorNet({self.r1});
                net.rtol = rtol;
                net.atol = atol;

                self.verifyEqual(net.rtol, rtol);
                self.verifyEqual(net.atol, atol);
                tEnd = 1.0;
                nSteps = 0;
                t = 0;

                while t < tEnd
                    t = net.step();
                    nSteps = nSteps + 1;
                end

                n = nSteps;
            end

            nbaseline = integrate(1e-10, 1e-20);
            nrtol = integrate(1e-8, 1e-20);
            natol = integrate(1e-10, 1e-5);

            self.verifyGreaterThan(nbaseline, nrtol);
            self.verifyGreaterThan(nbaseline, natol);
        end

        function testAdvanceReverse(self)
            self.makeReactors('nr', 1);
            net = ct.zeroD.ReactorNet(self.r1);
            net.advance(0.1);

            try
                net.advance(0.09);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'backwards in time');
            end
        end

        function testEquilibriumUV(self)
            P0 = 10 * ct.OneAtm;
            T0 = 1100;
            X0 = 'H2:1.0, O2:0.5, AR:8.0';
            self.makeReactors('nr', 1, 'T1', T0, 'P1', P0, 'X1', X0);
            net = ct.zeroD.ReactorNet({self.r1});

            net.advance(1.0);

            gas = ct.Solution('h2o2.yaml', '', 'none');
            gas.TPX = {T0, P0, X0};
            gas.equilibrate('UV');
            gas.basis = 'mass';

            self.verifyEqual(self.r1.T, gas.T, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.D, gas.D, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.P, gas.P, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.phase.X, gas.X, 'RelTol', self.rtol);
        end

        function testEquilibriumHP(self)
            P0 = 10 * ct.OneAtm;
            T0 = 1100;
            X0 = 'H2:1.0, O2:0.5, AR:8.0';

            self.gas1 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {T0, P0, X0};
            self.r1 = ct.zeroD.IdealGasConstPressureReactor(self.gas1);

            net = ct.zeroD.ReactorNet(self.r1);
            net.time = 0.0;
            net.advance(1.0);

            self.gas2 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas2.TPX = {T0, P0, X0};
            self.gas2.equilibrate('HP');
            self.gas2.basis = 'mass';

            self.verifyEqual(self.r1.T, self.gas2.T, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.D, self.gas2.D, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.P, P0, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.phase.X, self.gas2.X, 'RelTol', self.rtol);
        end

        function testWallVelocity(self)
            self.assumeFail('Skipped until velocity getter is implemented');

            self.makeReactors();
            A = 0.2;
            V1 = 2.0;
            V2 = 5.0;

            self.r1.V = V1;
            self.r2.V = V2;
            self.addWall('A', A);

            v = ct.Func1('tabulated-linear', [0.0, 1.0, 2.0], [0.0, 1.0, 0.0]);
            self.w.velocity = v;
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.advance(1.0);

            self.verifyEqual(self.w.velocity, v(1.0), 'AbsTol', self.atol);
            self.verifyEqual(self.w.expansionRate, 1.0 * A, 'AbsTol', self.atol);

            net.advance(2.0);

            self.verifyEqual(self.expansionRate, 0.0, 'AbsTol', self.atol);
            self.verifyEqual(self.r1.V, V1 + 1.0 * A, 'RelTol', self.rtol);
            self.verifyEqual(self.r2.V, V2 - 1.0 * A, 'RelTol', self.rtol);
        end

        function testDisableEnergy(self)
            self.makeReactors('T1', 500);
            self.r1.energyEnabled = false;
            self.addWall('A', 1.0, 'U', 2500);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.advance(11.0);

            self.verifyEqual(self.r1.T, 500, 'RelTol', self.rtol);
            self.verifyEqual(self.r2.T, 500, 'RelTol', self.rtol);
        end

        function testDisableChemistry(self)
            self.makeReactors('T1', 1000, 'nr', 1, 'X1', 'H2:2.0,O2:1.0');
            self.r1.chemistryEnabled = false;
            net = ct.zeroD.ReactorNet({self.r1});
            net.advance(11.0);
            x1 = self.r1.phase.X(self.r1.phase.speciesIndex('H2'));
            x2 = self.r1.phase.X(self.r1.phase.speciesIndex('O2'));

            self.verifyEqual(self.r1.T, 1000, 'RelTol', self.rtol);
            self.verifyEqual(x1, 2.0 / 3.0, 'RelTol', self.rtol);
            self.verifyEqual(x2, 1.0 / 3.0, 'RelTol', self.rtol);
        end

        function testHeatFluxFunc(self)
            self.assumeFail('Skipped until Func1.heatFlux getter is implemented');
            self.makeReactors('T1', 500, 'T2', 300);
            self.r1.V = 0.5;

            U1a = self.r1.V * self.r1.D * self.r1.phase.U;
            U2a = self.r2.V * self.r2.D * self.r2.phase.U;
            V1a = self.r1.V;
            V2a = self.r2.V;

            self.addWall('A', 0.3);
            f = ct.Func1('polynomial3', [-90000, 0, 90000]);
            self.w.heatFlux = f;
            Q = 0.3 * 60000;

            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.advance(1.1);
            self.verifyEqual(self.w.heatFlux, f(1.1), 'RelTol', self.rtol);
            U1b = self.r1.V * self.r1.D * self.r1.phase.U;
            U2b = self.r2.V * self.r2.D * self.r2.phase.U;

            self.verifyEqual(V1a, self.r1.V, 'RelTol', self.rtol);
            self.verifyEqual(V2a, self.r2.V, 'RelTol', self.rtol);
            self.verifyEqual(U1a - Q, U1b, 'RelTol', self.rtol);
            self.verifyEqual(U2a + Q, U2b, 'RelTol', self.rtol);
        end

        function testReinitialize(self)
            self.assumeFail('Skipped until Reactor.syncState is implemented');
        end

        function testPreconditionerUnsupported(self)
            self.assumeFail('Skipped until ReactorNet.preconditioner is implemented');
        end

        function testAddSensitivityReactionIndexing(self)
            self.gas1 = ct.Solution('h2o2.yaml', '', 'none');
            self.r1 = ct.zeroD.IdealGasReactor(self.gas1);
            net = ct.zeroD.ReactorNet(self.r1);

            % The last reaction is in range; 0 and one past the end are not.
            n = self.gas1.nReactions;
            self.r1.addSensitivityReaction(n);
            self.verifyError(@() self.r1.addSensitivityReaction(0), ...
                             'Cantera:ctError');
            self.verifyError(@() self.r1.addSensitivityReaction(n + 1), ...
                             'Cantera:ctError');
        end

        function testSensitivityComponentForms(self)
            net = self.makeSensitivityNet();

            s = net.sensitivity('temperature', 1);
            self.verifyTrue(isfinite(s));
            self.verifyNotEqual(s, 0);

            % A char array, a string, and the 1-based index of 'temperature' within the
            % global state vector all name the same component.
            self.verifyEqual(net.sensitivity("temperature", 1), s);
            k = net.globalComponentIndex('temperature');
            self.verifyEqual(net.sensitivity(k, 1), s);

            % The reactor may be omitted, given as an object, or given as its 1-based
            % position within the network.
            self.verifyEqual(net.sensitivity('temperature', 1, self.r1), s);
            self.verifyEqual(net.sensitivity('temperature', 1, 1), s);
        end

        function testSensitivitySpeciesComponent(self)
            net = self.makeSensitivityNet();

            s = net.sensitivity('OH', 1);
            self.verifyTrue(isfinite(s));
            self.verifyNotEqual(s, 0);

            k = net.globalComponentIndex('OH');
            self.verifyEqual(net.sensitivity(k, 1), s);
        end

        function testSensitivityReactorArgument(self)
            % In a network of two disjoint reactors, only the second one carries
            % a sensitivity parameter, so the component name has to be resolved
            % within the requested reactor to get the right answer.
            self.gas1 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {300, ct.OneAtm, 'H2:2.0, O2:1.0, AR:4.0'};
            self.r1 = ct.zeroD.IdealGasReactor(self.gas1);

            self.gas2 = ct.Solution('h2o2.yaml', '', 'none');
            self.gas2.TPX = {1050, 5 * ct.OneAtm, 'H2:2.0, O2:1.0, AR:4.0'};
            self.r2 = ct.zeroD.IdealGasReactor(self.gas2);

            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            self.r2.addSensitivityReaction(self.sensitivityReaction(self.gas2));
            net.setSensitivityTolerances(1e-6, 1e-6);
            net.advance(0.1);

            s2 = net.sensitivity('temperature', 1, self.r2);
            self.verifyTrue(isfinite(s2));
            self.verifyNotEqual(s2, 0);
            self.verifyEqual(net.sensitivity('temperature', 1, 2), s2);

            % The state of the second reactor follows that of the first.
            k2 = net.globalComponentIndex('temperature', self.r2);
            self.verifyEqual(net.sensitivity(k2, 1), s2);

            % The first reactor is unaffected by the parameter, and is what the
            % default reactor argument selects.
            self.verifyEqual(net.sensitivity('temperature', 1), 0, ...
                             'AbsTol', self.atol);
            k1 = net.globalComponentIndex('temperature');
            self.verifyEqual(net.sensitivity(k1, 1), 0, 'AbsTol', self.atol);
        end

        function testSensitivityInvalidComponent(self)
            net = self.makeSensitivityNet();

            try
                net.sensitivity({'temperature'}, 1);
                self.verifyFail('Expected an error for a cell array component.');
            catch ME
                self.verifySubstring(ME.message, 'must be a character array');
                self.verifySubstring(ME.message, 'cell');
            end

            try
                net.sensitivity('spam', 1);
                self.verifyFail('Expected an error for an unknown component.');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'spam');
            end
        end

        function testSensitivityInvalidReactor(self)
            net = self.makeSensitivityNet();

            gas = ct.Solution('h2o2.yaml', '', 'none');
            stranger = ct.zeroD.IdealGasReactor(gas);

            try
                net.sensitivity('temperature', 1, stranger);
                self.verifyFail('Expected an error for a foreign reactor.');
            catch ME
                self.verifySubstring(ME.message, 'not part of this network');
            end

            try
                net.sensitivity('temperature', 1, 'first');
                self.verifyFail('Expected an error for a char reactor argument.');
            catch ME
                self.verifySubstring(ME.message, 'Reactor must be a');
            end
        end

        function testSensitivityInvalidParameter(self)
            net = self.makeSensitivityNet();

            % Only one sensitivity parameter was registered.
            self.verifyError(@() net.sensitivity('temperature', 2), ...
                             'Cantera:ctError');
            k = net.globalComponentIndex('temperature');
            self.verifyError(@() net.sensitivity(k, 2), 'Cantera:ctError');
        end

        function testReactorComponents(self)
            gas = ct.Solution('h2o2.yaml', '', 'none');
            r = ct.zeroD.IdealGasReactor(gas);
            net = ct.zeroD.ReactorNet(r);

            % The state vector of an IdealGasReactor is
            % [mass, volume, temperature, mass fractions...]
            self.verifyEqual(r.nVars, 3 + gas.nSpecies);
            self.verifyEqual(net.nVars, r.nVars);
            self.verifyEqual(r.componentName(3), 'temperature');
            self.verifyEqual(r.componentIndex('temperature'), 3);
            self.verifyEqual(r.componentIndex('OH'), 3 + gas.speciesIndex('OH'));

            % Names and 1-based indices round-trip for every component
            for k = 1:r.nVars
                self.verifyEqual(r.componentIndex(r.componentName(k)), k);
            end

            % Indices outside the 1-based range are rejected
            self.verifyError(@() r.componentName(0), ...
                             'MATLAB:validators:mustBePositive');
            self.verifyError(@() r.componentName(r.nVars + 1), 'Cantera:ctError');

            % So is an unknown component name
            self.verifyError(@() r.componentIndex('spam'), 'Cantera:ctError');
        end

        function testReactorNetComponents(self)
            gas = ct.Solution('h2o2.yaml', '', 'none');
            first = ct.zeroD.IdealGasReactor(gas, 'first');
            second = ct.zeroD.IdealGasReactor(gas, 'second');
            net = ct.zeroD.ReactorNet({first, second});

            self.verifyEqual(net.nVars, first.nVars + second.nVars);

            % Components of the second reactor are offset by the length of the first
            k = net.globalComponentIndex('temperature', 2);
            self.verifyEqual(k, first.nVars + second.componentIndex('temperature'));

            % The reactor may be omitted, given as an object, or given as its
            % 1-based position within the network
            self.verifyEqual(net.globalComponentIndex('temperature', second), k);
            self.verifyEqual(net.globalComponentIndex('temperature'), ...
                             first.componentIndex('temperature'));

            % ReactorNet.componentName returns the reactor-prefixed name, while
            % globalComponentIndex takes the unprefixed one
            self.verifyEqual(net.componentName(k), 'second: temperature');
            self.verifyEqual(second.componentName( ...
                             second.componentIndex('temperature')), 'temperature');
            self.verifyError(@() net.globalComponentIndex('second: temperature', 2), ...
                             'Cantera:ctError');

            % Indices outside the 1-based range are rejected
            self.verifyError(@() net.componentName(0), ...
                             'MATLAB:validators:mustBePositive');
            self.verifyError(@() net.componentName(net.nVars + 1), 'Cantera:ctError');

            % So is a reactor that is not in the network, in either form
            try
                net.globalComponentIndex('temperature', 3);
                self.verifyFail('Expected an error for an out-of-range reactor.');
            catch ME
                self.verifySubstring(ME.message, 'between 1 and 2');
            end

            stranger = ct.zeroD.IdealGasReactor(gas);
            try
                net.globalComponentIndex('temperature', stranger);
                self.verifyFail('Expected an error for a foreign reactor.');
            catch ME
                self.verifySubstring(ME.message, 'not part of this network');
            end
        end

        function testReservoirComponentIndex(self)
            % ReactorBase::componentIndex and componentName are not implemented for
            % reactor types that have no state vector of their own.
            gas = ct.Solution('h2o2.yaml', '', 'none');
            res = ct.zeroD.Reservoir(gas);

            self.verifyEqual(res.nVars, 0);
            self.verifyError(@() res.componentIndex('temperature'), ...
                             'Cantera:ctError');
            self.verifyError(@() res.componentName(1), 'Cantera:ctError');
        end

        function testInvalidProperty(self)
            self.makeReactors();

            try
                self.r1.foobar = 3.14;
            catch ME
                self.verifySubstring(ME.identifier, 'MATLAB:');
                self.verifySubstring(ME.message, 'Unrecognized property');
            end
        end

    end

end
