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

        function deleteSolution(self)
            clear self.gas1 self.gas2 self.r1 self.r2 self.net
        end

    end

    methods

        function makeReactors(self, varargin)
            p = inputParser;

            p.addOptional('independent', true, ...
                           @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
            p.addOptional('nr', 2, ...
                           @(x) (isnumeric(x) && isscalar(x)));
            p.addOptional('cond1', {300, 101325, 'O2:1.0'}, @(x) iscell(x));
            p.addOptional('cond2', {300, 101325, 'O2:1.0'}, @(x) iscell(x));

            parse(p, varargin{:});

            independent = p.Results.independent;
            nr = p.Results.nr;
            cond1 = p.Results.cond1;
            cond2 = p.Results.cond2;

            self.net = ReactorNet();
            self.verifyEqual(self.net.time, 0, 'AbsTol', self.atol);

            self.gas1 = Solution('h2o2.yaml', '', 'none');
            self.gas1.TPX = cond1;
            self.r1 = Reactor(self.gas1);
            self.net.addReactor(self.r1);

            if independent
                self.gas2 = Solution('h2o2.yaml', '', 'none');
            else
                self.gas2 = self.gas1;
            end

            if nr >= 2
                self.gas2.TPX = cond2;
                self.r2 = Reactor(self.gas2);
                self.net.addReactor(self.r2);
            end
        end

        function addWall(self, varargin)
            p = inputParser;

            p.addOptional('K', 0.0, @(x) isnumeric(x) && isscalar(x));
            p.addOptional('A', 1.0, @(x) isnumeric(x) && isscalar(x));

            parse(p, varargin{:});

            K = p.Results.K;
            A = p.Results.A;

            self.w = Wall(self.r1, self.r2);
            self.w.area = A;
            self.w.heatTransferCoeff = K;
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
            self.makeReactors(false, 1);
            self.verifyEqual(self.net.time, 0.0, 'AbsTol', self.atol);
        end

        function testDisjoint1(self)
            cond1 = {300, 101325, 'O2:1.0'};
            cond2 = {500, 300000, 'O2:1.0'};

            self.makeReactors(true, 2, cond1, cond2);
            self.net.advance(1.0);

            self.verifyEqual(self.gas1.T, 300, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.T, 500, 'AbsTol', self.atol);
            self.verifyEqual(self.gas1.P, 101325, 'AbsTol', self.atol);
            self.verifyEqual(self.gas2.P, 300000, 'AbsTol', self.atol);
        end

        function testDisjoint2(self)
            cond1 = {300, 101325, 'O2:1.0'};
            cond2 = {500, 300000, 'O2:1.0'};

            self.makeReactors(false, 2, cond1, cond2);
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
            cond1 = {300, 101325, 'O2:1.0'};
            cond2 = {300, 300000, 'O2:1.0'};

            self.makeReactors(true, 2, cond1, cond2);
            self.addWall(0.1, 1.0);
            self.verifyEqual(self.w.type, 'Wall');
            self.w.name = 'name-of-wall';
            self.verifyEqual(self.w.name, 'name-of-wall');
        end

    end

end
