classdef ctTestWellStirredReactor < ctTestCase
    % Ignition (or not) of a well-stirred reactor

    properties
        gas
        fuel_in
        oxidizer_in
        exhaust
        combustor
        fuel_mfc
        oxidizer_mfc
        valve
        net
    end

    methods

        function makeReactors(self, T0, P0, mdot_fuel, mdot_ox)
            self.gas = Solution('../data/ch4_minimal.yaml', 'testWellStirredReactor');

            % fuel inlet
            self.gas.TPX = {T0, P0, 'CH4:1.0'};
            self.fuel_in = Reservoir(self.gas);

            % oxidizer inlet
            self.gas.TPX = {T0, P0, 'N2:3.76, O2:1.0'};
            self.oxidizer_in = Reservoir(self.gas);

            % reactor filled with N2
            self.gas.TPX = {T0, P0, 'N2:1.0'};
            self.combustor = IdealGasReactor(self.gas);
            self.combustor.V = 1.0;

            % outlet
            self.exhaust = Reservoir(self.gas);

            % connect the reactor to the reservois
            self.fuel_mfc = MassFlowController(self.fuel_in, self.combustor);
            self.fuel_mfc.massFlowRate = mdot_fuel;
            self.oxidizer_mfc = MassFlowController(self.oxidizer_in, self.combustor);
            self.oxidizer_mfc.massFlowRate = mdot_ox;
            self.valve = Valve(self.combustor, self.exhaust);
            self.valve.valveCoeff = 1.0;

            self.net = ReactorNet(self.combustor);
            % self.net.maxErrTestFails = 10;
        end

        function [times, T] = integrate(self, tf)
            t = 0.0;
            times = [];
            T = [];
            i = 0;

            while t < tf
                i = i + 1;
                t = self.net.step();
                times = [times, t];
                T = [T, self.combustor.T];
            end
        end

    end

    methods (Test)

        function testNonReacting(self)
            self.makeReactors(900.0, 10 * OneAtm, 1.0, 5.0);

            self.gas.setMultiplier(0.0);
            [t, T] = self.integrate(100.0);

            for i = 1:length(t)
                self.verifyEqual(T(i), 900, 'RelTol', 1e-5);
            end

            val1 = self.combustor.contents.Y...
                   (self.combustor.contents.speciesIndex('CH4'));
            val2 = (1.0 / 6.0);
            self.verifyEqual(val1, val2, 'RelTol', 1e-5);
        end

        function testIgnition1(self)
            self.makeReactors(900.0, 10 * OneAtm, 1.0, 5.0);

            [t, T] = self.integrate(10.0);

            self.verifyGreaterThan(T(end), 1200);
            for i = 1:length(t)
                if T(i) > 0.5 * (T(1) + T(end))
                    tIg = t(i);
                    break
                end
            end

            self.verifyEqual(tIg, 2.2249, 'RelTol', 1e-3);
        end

        function testIgnition2(self)
            self.makeReactors(900.0, 10 * OneAtm, 1.0, 20.0);

            [t, T] = self.integrate(10.0);

            self.verifyGreaterThan(T(end), 1200);
            for i = 1:length(t)
                if T(i) > 0.5 * (T(1) + T(end))
                    tIg = t(i);
                    break
                end
            end

            self.verifyEqual(tIg, 1.4856, 'RelTol', 1e-3);
        end

        function testIgnition3(self)
            self.makeReactors(900.0, 10 * OneAtm, 1.0, 80.0);
            self.net.maxTimeStep = 0.5;

            [t, T] = self.integrate(100.0);

            self.verifyLessThan(T(end), 910);
        end

        function testSteadyState(self)
            self.assumeFail('Skipped until advanceToSteadyState is implemented');
        end

    end

end
