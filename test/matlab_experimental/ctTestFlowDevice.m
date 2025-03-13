classdef ctTestFlowDevice < matlab.unittest.TestCase

    properties
        gas1
        gas2
        r1
        r2
        net
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

            param.parse(varargin{:});

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

    end

    methods (Test)

        function testMFC(self)
            self.assumeFail('Skipped until arbitrary functions can be set as Func1');
        end

        function testMFCType(self)
            self.makeReactors('nr', 2);
            mfc = MassFlowController(self.r1, self.r2);
            n1 = ReactorNet({self.r1, self.r2});

            self.verifyEqual(mfc.type, 'MassFlowController');
            self.verifyTrue(startsWith(mfc.name, 'MassFlowController_'));
            mfc.name = 'name-of-mfc';
            self.verifyEqual(mfc.name, 'name-of-mfc');

            clear n1 mfc
        end

        function testMFCErrors(self)
            self.assumeFail('Skipped until FlowDevice.pressureFunction is implemented');
        end

        function testValve1(self)
            self.makeReactors('P1', OneAtm, 'X1', 'AR:1.0', 'X2', 'O2:1.0');
            self.net.rtol = 1e-12;
            valve = Valve(self.r1, self.r2);
            k = 2e-5;
            valve.valveCoeff = k;
            self.net.time = 0;

            % Skipped until Reactor.inlets and Reactor.outlets are implemented
            % self.verifyEqual(self.r1.outlets, self.r2.inlets);

            % Skipped until Valve.valveCoeff getter is implemented
            % self.verifyEqual(valve.ValveCoeff, k);

            % self.verifyTrue(self.r1.energy);
            % self.verifyTrue(self.r2.energy);
            % self.verifyEqual((self.r1.P - self.r2.P) * k, valve.massFlowRate, ...
            %                  'RelTol', self.rtol);

            m1a = self.r1.D * self.r1.V;
            m2a = self.r2.D * self.r2.V;
            Y1a = self.r1.contents.Y;
            Y2a = self.r2.contents.Y;

            self.net.advance(0.1);

            m1b = self.r1.D * self.r1.V;
            m2b = self.r2.D * self.r2.V;
            Y1b = self.r1.contents.Y;
            Y2b = self.r2.contents.Y;

            % self.verifyEqual((self.r1.P - self.r2.P) * k, valve.massFlowRate, ...
            %                  'RelTol', self.rtol);
            self.verifyEqual(m1a + m2a, m1b + m2b, 'AbsTol', self.atol);
            self.verifyEqual(m1a .* Y1a + m2a .* Y2a, m1b .* Y1b + m2b .* Y2b, ...
                             'AbsTol', self.atol);
            self.verifyEqual(Y1a, Y1b, 'AbsTol', self.atol);

            clear valve
        end

        function testValve2(self)
            self.makeReactors('P1', OneAtm);
            self.net.rtol = 1e-11;
            self.r1.energy = 'off';
            self.r2.energy = 'off';
            valve = Valve(self.r1, self.r2);
            k = 2e-5;
            valve.valveCoeff = k;
            self.net.time = 0;

            % Skipped until Valve.valveCoeff getter is implemented
            % self.verifyEqual(valve.ValveCoeff, k);

            % self.verifyFalse(self.r1.energy);
            % self.verifyFalse(self.r2.energy);

            m1a = self.r1.D * self.r1.V;
            m2a = self.r2.D * self.r2.V;
            P1a = self.r1.P;
            P2a = self.r2.P;
            Y1 = self.r1.Y;
            A = k * P1a * (1 + m2a / m1a);
            B = k * (P1a / m1a + P2a / m2a);

            for t = linspace(1e-5, 0.5)
                self.net.advance(t);
                m1 = self.r1.D * self.r1.V;
                m2 = self.r2.D * self.r2.V;

                self.verifyEqual(m2, (m2a - A / B) * exp(-B * t) + A / B, ...
                                'AbsTol', self.atol);
                self.verifyEqual(m1a + m2a, m1 + m2, 'AbsTol', self.atol);
                self.verifyEqual(self.r1.Y, Y1, 'AbsTol', self.atol);
            end

            clear valve
        end

        function testValve3(self)
            self.assumeFail('Skipped until arbitrary functions can be set as Func1');
        end

        function testValveTiming(self)
            self.assumeFail('Skipped until time functions are implemented for Valve');
        end

        function testValveType1(self)
            self.makeReactors();
            res = Reservoir(self.gas1);
            v = Valve(self.r1, res);
            n1 = ReactorNet(self.r1);

            self.verifyTrue(startsWith(self.r1.name, sprintf('%s_', self.r1.type)));
            self.verifyTrue(startsWith(res.name, sprintf('%s_', res.type)));
            self.verifyEqual(v.type,'Valve');
            self.verifyTrue(startsWith(v.name, 'Valve_'));
            v.name = 'name-of-valve';
            self.verifyEqual(v.name, 'name-of-valve');

            clear res v n1
        end

        function testValveType2(self)
            self.makeReactors();
            res = Reservoir(self.gas1);
            v = Valve(res, self.r1);
            n1 = ReactorNet(self.r1);

            self.verifyTrue(startsWith(self.r1.name, sprintf('%s_', self.r1.type)));
            self.verifyTrue(startsWith(res.name, sprintf('%s_', res.type)));

            clear res v n1
        end

        function testValveErrors(self)
            self.makeReactors();
            v = Valve(self.r1, self.r2);

            try
                v.install(self.r2, self.r1);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Already installed');
            end

            clear v
        end

        function testPressureController(self)
            self.assumeFail('Skipped until PressureController is implemented');
        end

    end

end
