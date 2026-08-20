classdef ctTestFlowDevice < ctTestCase

    properties
        gas1
        gas2
        r1
        r2
    end

    properties (SetAccess = protected)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods

        function makeReactors(self, arg)
            arguments
                self
                arg.independent (1,1) logical = true
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

            if arg.independent
                self.gas2 = ct.Solution('h2o2.yaml', '', 'none');
            else
                self.gas2 = self.gas1;
            end

            if arg.nr >= 2
                self.gas2.TPX = {arg.T2, arg.P2, arg.X2};
                self.r2 = ct.zeroD.Reactor(self.gas2);
                self.r2.energyEnabled = true;
            end
        end

    end

    methods (Test)

        function testMFC(self)
            self.assumeFail('Skipped until arbitrary functions can be set as Func1');
        end

        function testMFCType(self)
            self.makeReactors('nr', 2);
            mfc = ct.zeroD.MassFlowController(self.r1, self.r2);
            net = ct.zeroD.ReactorNet({self.r1, self.r2});

            self.verifyEqual(mfc.type, 'MassFlowController');
            self.verifyTrue(startsWith(mfc.name, 'MassFlowController_'));
            mfc.name = 'name-of-mfc';
            self.verifyEqual(mfc.name, 'name-of-mfc');
        end

        function testMFCErrors(self)
            self.assumeFail('Skipped until FlowDevice.pressureFunction is implemented');
        end

        function testValve1(self)
            self.makeReactors('P1', ct.OneAtm, 'X1', 'AR:1.0', 'X2', 'O2:1.0');
            valve = ct.zeroD.Valve(self.r1, self.r2);
            k = 2e-5;
            valve.valveCoeff = k;
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.rtol = 1e-12;
            net.time = 0;

            % Skipped until Reactor.inlets and Reactor.outlets are implemented
            % self.verifyEqual(self.r1.outlets, self.r2.inlets);

            % Skipped until Valve.valveCoeff getter is implemented
            % self.verifyEqual(valve.ValveCoeff, k);

            % self.verifyTrue(self.r1.energyEnabled);
            % self.verifyTrue(self.r2.energyEnabled);
            % self.verifyEqual((self.r1.P - self.r2.P) * k, valve.massFlowRate, ...
            %                  'RelTol', self.rtol);

            m1a = self.r1.D * self.r1.V;
            m2a = self.r2.D * self.r2.V;
            Y1a = self.r1.phase.Y;
            Y2a = self.r2.phase.Y;

            net.advance(0.1);

            m1b = self.r1.D * self.r1.V;
            m2b = self.r2.D * self.r2.V;
            Y1b = self.r1.phase.Y;
            Y2b = self.r2.phase.Y;

            % self.verifyEqual((self.r1.P - self.r2.P) * k, valve.massFlowRate, ...
            %                  'RelTol', self.rtol);
            self.verifyEqual(m1a + m2a, m1b + m2b, 'AbsTol', self.atol);
            self.verifyEqual(m1a .* Y1a + m2a .* Y2a, m1b .* Y1b + m2b .* Y2b, ...
                             'AbsTol', self.atol);
            self.verifyEqual(Y1a, Y1b, 'AbsTol', self.atol);
        end

        function testValve2(self)
            self.makeReactors('P1', ct.OneAtm);
            net.rtol = 1e-11;
            self.r1.energyEnabled = false;
            self.r2.energyEnabled = false;
            valve = ct.zeroD.Valve(self.r1, self.r2);
            k = 2e-5;
            valve.valveCoeff = k;
            net = ct.zeroD.ReactorNet({self.r1, self.r2});
            net.time = 0;

            % Skipped until Valve.valveCoeff getter is implemented
            % self.verifyEqual(valve.ValveCoeff, k);

            % self.verifyFalse(self.r1.energyEnabled);
            % self.verifyFalse(self.r2.energyEnabled);

            m1a = self.r1.D * self.r1.V;
            m2a = self.r2.D * self.r2.V;
            P1a = self.r1.P;
            P2a = self.r2.P;
            Y1 = self.r1.Y;
            A = k * P1a * (1 + m2a / m1a);
            B = k * (P1a / m1a + P2a / m2a);

            for t = linspace(1e-5, 0.5)
                net.advance(t);
                m1 = self.r1.D * self.r1.V;
                m2 = self.r2.D * self.r2.V;

                self.verifyEqual(m2, (m2a - A / B) * exp(-B * t) + A / B, ...
                                'AbsTol', self.atol);
                self.verifyEqual(m1a + m2a, m1 + m2, 'AbsTol', self.atol);
                self.verifyEqual(self.r1.Y, Y1, 'AbsTol', self.atol);
            end
        end

        function testValve3(self)
            self.assumeFail('Skipped until arbitrary functions can be set as Func1');
        end

        function testValveTiming(self)
            self.assumeFail('Skipped until time functions are implemented for Valve');
        end

        function testValveType1(self)
            self.makeReactors();
            res = ct.zeroD.Reservoir(self.gas1);
            v = ct.zeroD.Valve(self.r1, res);
            net = ct.zeroD.ReactorNet(self.r1);

            self.verifyTrue(startsWith(self.r1.name, sprintf('%s_', self.r1.type)));
            self.verifyTrue(startsWith(res.name, sprintf('%s_', res.type)));
            self.verifyEqual(v.type,'Valve');
            self.verifyTrue(startsWith(v.name, 'Valve_'));
            v.name = 'name-of-valve';
            self.verifyEqual(v.name, 'name-of-valve');
        end

        function testValveType2(self)
            self.makeReactors();
            res = ct.zeroD.Reservoir(self.gas1);
            v = ct.zeroD.Valve(res, self.r1);
            net = ct.zeroD.ReactorNet(self.r1);

            self.verifyTrue(startsWith(self.r1.name, sprintf('%s_', self.r1.type)));
            self.verifyTrue(startsWith(res.name, sprintf('%s_', res.type)));
        end

        function testPressureController(self)
            self.makeReactors('nr', 1);
            g = ct.Solution('h2o2.yaml', '', 'none');
            g.TPX = {500, 2 * ct.OneAtm, 'H2:1.0'};
            inletReservoir = ct.zeroD.Reservoir(g);
            g.TP = {300, ct.OneAtm};
            outletReservoir = ct.zeroD.Reservoir(g);

            mdot = 0.1;
            mfc = ct.zeroD.MassFlowController(inletReservoir, self.r1);
            mfc.massFlowRate = mdot;

            pc = ct.zeroD.PressureController(self.r1, outletReservoir);
            pc.setPrimary(mfc);
            k = 2e-5;
            pc.pressureCoeff = k;
            self.verifyEqual(pc.pressureCoeff, k, 'RelTol', self.rtol);
            self.verifyEqual(pc.deviceCoefficient, k, 'RelTol', self.rtol);
            pc.pressureCoeff = 1e-5;
            k = 1e-5;
            self.verifyEqual(pc.pressureCoeff, k, 'RelTol', self.rtol);

            net = ct.zeroD.ReactorNet(self.r1);
            t = 0;

            while t < 1.0
                t = net.step();
                self.verifyEqual(mfc.massFlowRate, mdot, 'RelTol', self.rtol);
                dP = self.r1.P - outletReservoir.P;
                self.verifyEqual(pc.massFlowRate, max(mdot + k * dP, 0.0), ...
                                 'RelTol', self.rtol);
            end
        end

        function testPressureControllerType(self)
            self.makeReactors();
            res = ct.zeroD.Reservoir(self.gas1);
            mfc = ct.zeroD.MassFlowController(res, self.r1);
            mfc.massFlowRate = 0.6;
            pc = ct.zeroD.PressureController(self.r1, self.r2);
            pc.setPrimary(mfc);
            pc.pressureCoeff = 0.5;
            net = ct.zeroD.ReactorNet({self.r1, self.r2});

            self.verifyEqual(pc.type, 'PressureController');
            self.verifyTrue(startsWith(pc.name, 'PressureController_'));
            pc.name = 'name-of-pressure-controller';
            self.verifyEqual(pc.name, 'name-of-pressure-controller');
        end

        function testPressureControllerErrors(self)
            self.makeReactors();
            res = ct.zeroD.Reservoir(self.gas1);
            mfc = ct.zeroD.MassFlowController(res, self.r1);
            mfc.massFlowRate = 0.6;
            pc = ct.zeroD.PressureController(self.r1, self.r2);

            % Mass flow rate cannot be evaluated before the primary device is set
            self.verifyError(@() pc.massFlowRate, 'Cantera:ctError');

            % Only pressure controllers accept a primary device
            v = ct.zeroD.Valve(self.r1, self.r2);
            self.verifyError(@() v.setPrimary(mfc), ?MException);
        end

    end

end
