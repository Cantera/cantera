classdef ctTestReactor < ctTestCase

    properties
        gas1
        gas2
        r1
        r2
        net
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

            self.gas1 = Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {arg.T1, arg.P1, arg.X1};
            self.r1 = Reactor(self.gas1);

            if arg.nr == 1
                self.net = ReactorNet(self.r1);
            elseif arg.nr >= 2
                self.gas2 = Solution('h2o2.yaml', '', 'none');
                self.gas2.TPX = {arg.T2, arg.P2, arg.X2};
                self.r2 = Reactor(self.gas2);
                self.r2.energy = 'on';
                self.net = ReactorNet({self.r1, self.r2});
            end
            self.verifyEqual(self.net.time, 0, 'AbsTol', self.atol);
        end

        function addWall(self, arg)
            arguments
                self
                arg.K (1,1) double {mustBeNumeric} = 0.0
                arg.U (1,1) double {mustBeNumeric} = 0.0
                arg.A (1,1) double {mustBeNumeric} = 1.0
            end

            self.w = Wall(self.r1, self.r2);
            self.w.area = arg.A;
            self.w.expansionRateCoeff = arg.K;
            self.w.heatTransferCoeff = arg.U;
        end

    end

    methods (Test)

        function testV(self)
            g = Solution('h2o2.yaml', '', 'none');
            r = Reactor(g);
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
            self.verifyEqual(self.net.time, 0.0, 'AbsTol', self.atol);
        end

        function testDisjoint(self)
            self.makeReactors('T1', 300, 'P1', 101325, 'T2', 500, 'P2', 300000);
            self.net.advance(1.0);

            self.verifyEqual(self.gas1.T, 300, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.T, 500, 'AbsTol', self.atol);
            self.verifyEqual(self.gas1.P, 101325, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.P, 300000, 'AbsTol', self.atol);
        end

        function testTimeStepping(self)
            self.makeReactors();

            tStart = 0.3;
            tEnd = 10.0;
            dtMax = 0.07;
            t = tStart;

            self.net.maxTimeStep = dtMax;
            self.net.time = tStart;
            self.verifyEqual(self.net.time, tStart);

            while t < tEnd
                tPrev = t;
                t = self.net.step;
                self.verifyLessThanOrEqual(t - tPrev, 1.0001 * dtMax);
                self.verifyEqual(t, self.net.time, 'AbsTol', self.atol);
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
            res = Reservoir(self.gas1);
            w = Wall(self.r1, res);
            net = ReactorNet([self.r1]);
            self.verifyEqual(w.type, 'Wall');
        end

        function testWall3(self)
            self.makeReactors('nr', 1);
            res = Reservoir(self.gas1);
            w = Wall(res, self.r1);
            net = ReactorNet([self.r1]);
            self.verifyEqual(w.type, 'Wall');
        end

        function testEqualizePressure(self)
            self.makeReactors('P1', 101325, 'P2', 300000);
            self.addWall('K', 0.1, 'A', 1.0);

            self.net.advance(1.0);
            self.verifyEqual(self.net.time, 1.0, 'AbsTol', self.atol);
            self.verifyEqual(self.r1.phase.P, self.r2.phase.P, 'RelTol', self.rtol);
            self.verifyNotEqual(self.r1.T, self.r2.T);
        end

        function testHeatTransfer1(self)
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 500, 'A', 1.0);

            self.net.advance(10.0);
            self.verifyEqual(self.net.time, 10.0, 'RelTol', self.atol);
            self.verifyEqual(self.r1.T, self.r2.T, 'RelTol', 5e-7);
            self.verifyNotEqual(self.r1.P, self.r2.P);
        end

        function testHeatTransfer2(self)
            self.assumeFail('Skipped until heatTransferCoeff getter is implemented');
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 200, 'A', 1.0);

            self.net.advance(1.0);
            T1a = self.r1.T;
            T2a = self.r2.T;

            self.makeReactors('T1', 300, 'T2', 1000);
            self.r1.V = 0.25;
            self.r2.V = 0.25;
            self.addWall('U', 100, 'A', 0.5);

            Qdot1 = self.w.heatTransferCoeff * self.w.area * (self.r1.T - self.r2.T);
            self.verifyEqual(Qdot1, self.w.heatRate, 'RelTol', self.rtol);

            self.net.advance(1.0);
            Qdot2 = self.w.heatTransferCoeff * self.w.area * (self.r1.T - self.r2.T);
            self.verifyEqual(Qdot2, self.w.heatRate, 'RelTol', self.rtol);
            T1b = self.r1.T;
            T2b = self.r2.T;

            self.verifyEqual(T1a, T1b, 'RelTol', self.rtol);
            self.verifyEqual(T2a, T2b, 'RelTol', self.rtol);
        end

        function testTolerances(self)
            function n = integrate(atol, rtol)
                P0 = 10 * OneAtm;
                T0 = 1100;
                X0 = 'H2:1.0, O2:0.5, AR:8.0';
                self.makeReactors('nr', 1, 'T1', T0, 'P1', P0, 'X1', X0);
                self.net.rtol = rtol;
                self.net.atol = atol;

                self.verifyEqual(self.net.rtol, rtol);
                self.verifyEqual(self.net.atol, atol);

                tEnd = 1.0;
                nSteps = 0;
                t = 0;

                while t < tEnd
                    t = self.net.step();
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
            self.net.advance(0.1);

            try
                self.net.advance(0.09);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'backwards in time');
            end
        end

        function testEquilibriumUV(self)
            P0 = 10 * OneAtm;
            T0 = 1100;
            X0 = 'H2:1.0, O2:0.5, AR:8.0';
            self.makeReactors('nr', 1, 'T1', T0, 'P1', P0, 'X1', X0);

            self.net.advance(1.0);

            gas = Solution('h2o2.yaml', '', 'none');
            gas.TPX = {T0, P0, X0};
            gas.equilibrate('UV');
            gas.basis = 'mass';

            self.verifyEqual(self.r1.T, gas.T, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.D, gas.D, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.P, gas.P, 'RelTol', self.rtol);
            self.verifyEqual(self.r1.phase.X, gas.X, 'RelTol', self.rtol);
        end

        function testEquilibriumHP(self)
            P0 = 10 * OneAtm;
            T0 = 1100;
            X0 = 'H2:1.0, O2:0.5, AR:8.0';

            self.gas1 = Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {T0, P0, X0};
            self.r1 = IdealGasConstPressureReactor(self.gas1);

            self.net = ReactorNet(self.r1);
            self.net.time = 0.0;
            self.net.advance(1.0);

            self.gas2 = Solution('h2o2.yaml', '', 'none');
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

            v = Func1('tabulated-linear', [0.0, 1.0, 2.0], [0.0, 1.0, 0.0]);
            self.w.velocity = v;
            self.net.advance(1.0);

            self.verifyEqual(self.w.velocity, v(1.0), 'AbsTol', self.atol);
            self.verifyEqual(self.w.expansionRate, 1.0 * A, 'AbsTol', self.atol);

            self.net.advance(2.0);

            self.verifyEqual(selfw.expansionRate, 0.0, 'AbsTol', self.atol);
            self.verifyEqual(self.r1.V, V1 + 1.0 * A, 'RelTol', self.rtol);
            self.verifyEqual(self.r2.V, V2 - 1.0 * A, 'RelTol', self.rtol);
        end

        function testDisableEnergy(self)
            self.makeReactors('T1', 500);
            self.r1.energy = 'off';
            self.addWall('A', 1.0, 'U', 2500);
            self.net.advance(11.0);

            self.verifyEqual(self.r1.T, 500, 'RelTol', self.rtol);
            self.verifyEqual(self.r2.T, 500, 'RelTol', self.rtol);
        end

        function testDisableChemistry(self)
            self.makeReactors('T1', 1000, 'nr', 1, 'X1', 'H2:2.0,O2:1.0');
            self.r1.chemistry = 'off';
            self.net.advance(11.0);
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
            f = Func1('polynomial3', [-90000, 0, 90000]);
            self.w.heatFlux = f;
            Q = 0.3 * 60000;

            self.net.advance(1.1);
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
