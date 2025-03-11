classdef ctTestReactor < matlab.unittest.TestCase

    properties
        gas1
        gas2
        r1
        r2
        net
        w
    end

    properties (SetAccess = immutable)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
            copyfile('../data/equilibrium.yaml', './equilibrium.yaml');
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            delete('./equilibrium.yaml');
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodTeardown)

        function deleteVariable(self)
            clear self.gas1 self.gas2 self.r1 self.r2 self.net self.w
        end

    end

    methods

        function makeReactors(self, varargin)
            param = inputParser;

            param.addParameter('independent', true, ...
                               @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
            param.addParameter('nr', 2, @(x) (isnumeric(x) && isscalar(x)));
            param.addParameter('T1', 300, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('P1', 101325, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('X1', 'O2:1.0', @(x) ischar(x));
            param.addParameter('T2', 300, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('P2', 101325, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('X2', 'O2:1.0', @(x) ischar(x));

            parse(param, varargin{:});

            independent = param.Results.independent;
            nr = param.Results.nr;
            T1 = param.Results.T1;
            P1 = param.Results.P1;
            X1 = param.Results.X1;
            T2 = param.Results.T2;
            P2 = param.Results.P2;
            X2 = param.Results.X2;

            self.net = ReactorNet();
            self.verifyEqual(self.net.time, 0, 'AbsTol', self.atol);

            self.gas1 = Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = {T1, P1, X1};
            self.r1 = Reactor(self.gas1);
            self.net.addReactor(self.r1);

            if independent
                self.gas2 = Solution('h2o2.yaml', '', 'none');
            else
                self.gas2 = self.gas1;
            end

            if nr >= 2
                self.gas2.TPX = {T2, P2, X2};
                self.r2 = Reactor(self.gas2);
                self.r2.energy = 'on';
                self.net.addReactor(self.r2);
            end
        end

        function addWall(self, varargin)
            param = inputParser;

            param.addParameter('K', 0.0, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('U', 0.0, @(x) isnumeric(x) && isscalar(x));
            param.addParameter('A', 1.0, @(x) isnumeric(x) && isscalar(x));

            parse(param, varargin{:});

            K = param.Results.K;
            U = param.Results.U;
            A = param.Results.A;

            self.w = Wall(self.r1, self.r2);
            self.w.area = A;
            self.w.expansionRateCoeff = K;
            self.w.heatTransferCoeff = U;
        end

    end

    methods (Test)

        function testVolume(self)
            g = Solution('h2o2.yaml', '', 'none');
            r = Reactor(g, '', '');
            self.verifyEqual(r.V, 1.0, 'AbsTol', self.atol);

            r.V = 9;
            self.verifyEqual(r.V, 9.0, 'AbsTol', self.atol);

            clear g r
        end

        function testTypes(self)
            self.makeReactors();
            self.verifyEqual(self.r1.type, 'Reactor');
        end

        function testIndependentVariable(self)
            self.makeReactors('independent', false, 'nr', 1);
            self.verifyEqual(self.net.time, 0.0, 'AbsTol', self.atol);
        end

        function testDisjoint1(self)
            self.makeReactors('T1', 300, 'P1', 101325, 'T2', 500, 'P2', 300000);
            self.net.advance(1.0);

            self.verifyEqual(self.gas1.T, 300, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.T, 500, 'AbsTol', self.atol);
            self.verifyEqual(self.gas1.P, 101325, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.P, 300000, 'AbsTol', self.atol);
        end

        function testDisjoint2(self)
            self.makeReactors('independent', false, ...
                              'T1', 300, 'P1', 101325, 'T2', 500, 'P2', 300000);
            self.net.advance(1.0);

            self.verifyEqual(self.r1.T, 300, 'AbsTol', self.atol);
            self.verifyEqual(self.r2.T, 500, 'AbsTol', self.atol);
            self.verifyEqual(self.r1.P, 101325, 'AbsTol', self.atol);
            self.verifyEqual(self.r2.P, 300000, 'AbsTol', self.atol);
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
                t = self.net.dt;
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

            clear res w net
        end

        function testWall3(self)
            self.makeReactors('nr', 1);
            res = Reservoir(self.gas1);
            w = Wall(res, self.r1);
            net = ReactorNet([self.r1]);
            self.verifyEqual(w.type, 'Wall');

            clear res w net
        end

        function testEqualizePressure(self)
            self.makeReactors('P1', 101325, 'P2', 300000);
            self.addWall('K', 0.1, 'A', 1.0);

            self.net.advance(1.0);
            self.verifyEqual(self.net.time, 1.0, 'AbsTol', self.atol);
            self.verifyEqual(self.gas1.P, self.gas2.P, 'RelTol', self.rtol);
            self.verifyNotEqual(self.r1.T, self.r2.T);
        end

        function testHeatTransfer1(self)
            self.assumeFail('Skipped since the two sides of the wall could not equilibrate in 10 s');
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 500, 'A', 1.0);

            self.net.advance(10.0);
            self.verifyEqual(self.net.time, 10.0, 'RelTol', self.atol);
            self.verifyEqual(self.r1.T, self.r2.T, 'RelTol', 5e-7);
            self.verifyNotEqual(self.r1.P, self.r2.P);
        end

        function testHeatTransfer2(self)
            self.assumeFail('Skipping until heatTransferCoeff getter is implemented');
            self.makeReactors('T1', 300, 'T2', 1000);
            self.addWall('U', 200, 'A', 1.0);

            self.net.advance(1.0);
            T1a = self.r1.T;
            T2a = self.r2.T;

            % Delete old objects to prevent illegal memory access
            delete(self.gas1);
            delete(self.gas2);
            delete(self.r1);
            delete(self.r2);
            delete(self.net);
            delete(self.w);

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

    end

end
