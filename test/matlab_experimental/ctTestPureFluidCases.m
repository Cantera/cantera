classdef ctTestPureFluidCases < ctTestCase

    properties
        fluid
        states
        refState
        tol
        u0
        s0
    end

    methods

        function setupFluid(self, type)
            self.fluid = Solution('liquidvapor.yaml', type);
            self.fluid.basis = 'mass';
            self.fluid.TD = {self.refState.T, self.refState.D};
            self.u0 = self.fluid.U;
            self.s0 = self.fluid.S;
        end

        function val = a(self, T, rho)
            % Helmholtz free energy
            self.fluid.TD = {T, rho};
            val = self.fluid.U - T * self.fluid.S;
        end

        function testConsistencyTemperature(self)
            for state = self.states
                dT = 2e-5 * state.T;
                self.fluid.TD = {state.T - dT, state.D};
                s1 = self.fluid.S;
                u1 = self.fluid.U;
                self.fluid.TD = {state.T + dT, state.D};
                s2 = self.fluid.S;
                u2 = self.fluid.U;

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual((u2 - u1) / (s2 - s1), state.T, ...
                                 'RelTol', self.tol.dUdS, msg);
            end
        end

        function testConsistencyVolume(self)
            for state = self.states
                self.fluid.TD = {state.T, state.D};
                p1 = self.fluid.P;
                V = 1 / state.D;
                dV = 5e-6 * V;

                a1 = self.a(state.T, 1 / (V - 0.5 * dV));
                a2 = self.a(state.T, 1 / (V + 0.5 * dV));

                % dP/drho is high for liquids, so relax tolerances
                if strcmp(state.phase, 'liquid')
                    tol = 300 * self.tol.dAdV;
                else
                    tol = self.tol.dAdV;
                end

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual(-(a2 - a1) / dV, p1, 'RelTol', tol, msg);
            end
        end

        function testSaturation(self)
            for state = self.states
                if strcmp(state.phase, 'super')
                    continue
                end

                dT = 1e-6 * state.T;
                self.fluid.TQ = {state.T, 0};
                p1 = self.fluid.P;
                vf = 1.0 / self.fluid.D;
                hf = self.fluid.H;
                sf = self.fluid.S;

                self.fluid.TQ = {state.T + dT, 0};
                p2 = self.fluid.P;

                self.fluid.TQ = {state.T, 1};
                vg = 1.0 / self.fluid.D;
                hg = self.fluid.H;
                sg = self.fluid.S;

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual((p2 - p1) / dT, (hg - hf)/(state.T * (vg - vf)), ...
                                 'RelTol', self.tol.dPdT, msg);
                self.verifyEqual(hg - hf, state.T * (sg - sf), 'RelTol', self.tol.hTs);
            end
        end

        function testPressure(self)
            for state = self.states
                self.fluid.TD = {state.T, state.D};
                if strcmp(state.phase, 'liquid')
                    tol = 70 * self.tol.P;
                else
                    tol = self.tol.P;
                end
                tol = tol * state.tolMod;

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual(self.fluid.P, state.P, 'RelTol', tol, msg);
            end
        end

        function testInternalEnergy(self)
            for state = self.states
                self.fluid.TD = {state.T, state.D};
                tol = self.tol.U * state.tolMod;

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual(self.fluid.U - self.u0, state.U - self.refState.U, ...
                                 'RelTol', tol, msg);
            end
        end

        function testEntropy(self)
            for state = self.states
                self.fluid.TD = {state.T, state.D};
                tol = self.tol.S * state.tolMod;

                msg = sprintf('At state: T = %.2f, rho = %.2f', state.T, state.D);
                self.verifyEqual(self.fluid.U - self.u0, state.U - self.refState.U, ...
                                 'RelTol', tol, msg);
            end
        end

        function runAllTests(self)
            self.testConsistencyTemperature;
            self.testConsistencyVolume;
            self.testSaturation;
            self.testPressure;
            self.testInternalEnergy;
            self.testEntropy;
        end

    end

    methods (Test)

        function testHFC134a(self)
            self.states = [StateData('liquid', 175.0, 0.1, 'D', 1577.6239, ...
                           'U', 77.534586, 'S', 0.44788182),
                           StateData('liquid', 210.0, 0.1, 'D', 1483.2128, ...
                           'U', 119.48566, 'S', 0.66633877),
                           StateData('vapor', 250.0, 0.1, 'D', 5.1144317, ...
                           'U', 365.59424, 'S', 1.7577491),
                           StateData('vapor', 370.0, 0.1, 'D', 3.3472612, ...
                           'U', 459.82664, 'S', 2.0970769),
                           StateData('liquid', 290.0, 10, 'D', 1278.4700, ...
                           'U', 216.99119, 'S', 1.0613409),
                           StateData('super', 410.0, 10, 'D', 736.54666, ...
                           'U', 399.02258, 'S', 1.5972395),
                           StateData('super', 450.0, 40, 'D', 999.34087, ...
                           'U', 411.92422, 'S', 1.6108568)];
            self.states = self.states';
            self.refState = StateData('critical', 374.21, 4.05928, 'D', 511.900, ...
                                      'U', 381.70937, 'S', 1.5620991);
            self.tol = Tolerances();

            self.setupFluid('HFC-134a');
            self.runAllTests;

        end

        function testCarbonDioxide(self)
            self.states = [StateData('liquid', 230.0, 2.0, ...
                           'D', 1132.4, 'H', 28.25, 'S', 0.1208);
                           StateData('liquid', 270.0, 10.0, ...
                           'D', 989.97, 'H', 110.59, 'S', 0.4208);
                           StateData('vapor', 250.0, 1.788, ...
                           'V', 0.02140, 'H', 358.59, 'S', 1.4500, 'relax', true); % sat
                           StateData('vapor', 300.0, 2.0, ...
                           'V', 0.02535, 'H', 409.41, 'S', 1.6174);
                           StateData('super', 500.0, 1.0, ...
                           'V', 0.09376, 'H', 613.22, 'S', 2.2649);
                           StateData('super', 600.0, 20.0, ...
                           'V', 0.00554, 'H', 681.94, 'S', 1.8366)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 304.21, 7.3834, ...
                                      'D', 464.0, 'H', 257.31, 'S', 0.9312);
            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('carbon-dioxide');
            self.runAllTests;
        end

        function testHeptane(self)
            self.states = [StateData('liquid', 300.0, 0.006637, ...
                           'V', 0.001476, 'H', 0.0, 'S', 0.0, 'relax', true),
                           StateData('liquid', 400.0, 0.2175, ...
                           'V', 0.001712, 'H', 248.01, 'S', 0.709, 'relax', true),
                           StateData('vapor', 490.0, 1.282, ...
                           'V', 0.02222, 'H', 715.64, 'S', 1.7137, 'relax', true),
                           StateData('vapor', 480.0, 0.70, ...
                           'V', 0.04820, 'H', 713.04, 'S', 1.7477),
                           StateData('super', 600.0, 2.0, ...
                           'V', 0.01992, 'H', 1014.87, 'S', 2.2356),
                           StateData('super', 680.0, 0.2, ...
                           'V', 0.2790, 'H', 1289.29, 'S', 2.8450)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 537.68, 2.6199, ...
                                    'D', 197.60, 'H', 747.84, 'S', 1.7456);
            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('heptane');
            self.runAllTests;
        end

        function testHydrogen(self)
            self.states = [StateData('liquid', 18.0, 0.04807, ...
                           'V', 0.013660, 'H', 30.1, 'S', 1.856, 'relax', true),
                           StateData('liquid', 26.0, 0.4029, ...
                           'V', 0.015911, 'H', 121.2, 'S', 5.740, 'relax', true),
                           StateData('vapor', 30.0, 0.8214, ...
                           'V', 0.09207, 'H', 487.4, 'S', 17.859, 'relax', true),
                           StateData('super', 100.0, 0.20, ...
                           'V', 2.061, 'H', 1398.3, 'S', 39.869),
                           StateData('super', 200.0, 20.0, ...
                           'V', 0.04795, 'H', 3015.9, 'S', 31.274),
                           StateData('super', 300.0, 0.50, ...
                           'V', 2.482, 'H', 4511.6, 'S', 53.143),
                           StateData('super', 600.0, 1.00, ...
                           'V', 2.483, 'H', 8888.4, 'S', 60.398),
                           StateData('super', 800.0, 4.0, ...
                           'V', 0.8329, 'H', 11840.0, 'S', 58.890)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 32.938, 1.2838, ...
                                      'D', 31.36, 'H', 346.5, 'S', 12.536);

            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3, 'dUdS', 2e-4);

            self.setupFluid('hydrogen');
            self.runAllTests;
        end

        function testMethane(self)
            self.states = [StateData('liquid', 100.0, 0.50, ...
                           'D', 439.39, 'H', 31.65, 'S', 0.3206),
                           StateData('liquid', 140.0, 2.0, ...
                           'D', 379.51, 'H', 175.48, 'S', 1.4963),
                           StateData('vapor', 150.0, 0.20, ...
                           'V', 0.3772, 'H', 660.72, 'S', 5.5435),
                           StateData('vapor', 160.0, 1.594, ...
                           'V', 0.03932, 'H', 627.96, 'S', 4.3648, 'relax', true),
                           StateData('vapor', 175.0, 1.0, ...
                           'V', 0.08157, 'H', 692.55, 'S', 4.9558),
                           StateData('super', 200.0, 0.2, ...
                           'V', 0.5117, 'H', 767.37, 'S', 6.1574),
                           StateData('super', 300.0, 0.5, ...
                           'V', 0.3083, 'H', 980.87, 'S', 6.5513)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 190.555, 4.5988, ...
                                      'D', 160.43, 'H', 490.61, 'S', 3.2853);

            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('methane');
            self.runAllTests;
        end

        function testNitrogen(self)
            self.states = [StateData('liquid', 80.0, 0.1370, ...
                           'V', 0.001256, 'H', 33.50, 'S', 0.4668, 'relax', true),
                           StateData('vapor', 110.0, 1.467, ...
                           'V', 0.01602, 'H', 236.28, 'S', 2.3896, 'relax', true),
                           StateData('super', 200.0, 0.5, ...
                           'V', 0.1174, 'H', 355.05, 'S', 3.5019),
                           StateData('super', 300.0, 10.0, ...
                           'V', 0.00895, 'H', 441.78, 'S', 2.9797),
                           StateData('super', 500.0, 5.0, ...
                           'V', 0.03031, 'H', 668.48, 'S', 3.7722),
                           StateData('super', 600.0, 100.0, ...
                           'V', 0.00276, 'H', 827.54, 'S', 3.0208)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 126.200, 3.400, ...
                                      'D', 314.03, 'H', 180.78, 'S', 1.7903);

            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('nitrogen');
            self.runAllTests;
        end

        function testOxygen(self)
            self.states = [StateData('liquid', 80.0, 0.03009, ...
                           'V', 0.000840, 'H', 42.56, 'S', 0.6405, 'relax', true),
                           StateData('liquid', 125.0, 1.351, ...
                           'V', 0.001064, 'H', 123.24, 'S', 1.4236, 'relax', true),
                           StateData('vapor', 145.0, 3.448, ...
                           'V', 0.006458, 'H', 276.45, 'S', 2.4852, 'relax', true),
                           StateData('super', 200.0, 0.050, ...
                           'V', 1.038, 'H', 374.65, 'S', 4.1275),
                           StateData('super', 300.0, 1.0, ...
                           'V', 0.07749, 'H', 463.76, 'S', 3.7135),
                           StateData('super', 600.0, 0.20, ...
                           'V', 0.7798, 'H', 753.38, 'S', 4.7982),
                           StateData('super', 800.0, 5.0, ...
                           'V', 0.04204, 'H', 961.00, 'S', 4.2571)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 154.581, 5.0429, ...
                                      'D', 436.15, 'H', 226.53, 'S', 2.1080);

            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('oxygen');
            self.runAllTests;
        end

        function testWater(self)
            self.states = [StateData('liquid', 295.0, 0.002620, ...
                           'V', 0.0010025, 'H', 90.7, 'S', 0.3193, 'relax', true),
                           StateData('vapor', 315.0, 0.008143, ...
                           'V', 17.80, 'H', 2577.1, 'S', 8.2216, 'relax', true),
                           StateData('liquid', 440.0, 0.7332, ...
                           'V', 0.001110, 'H', 705.0, 'S', 2.0096, 'relax', true),
                           StateData('vapor', 510.0, 3.163, ...
                           'V', 0.06323, 'H', 2803.6, 'S', 6.1652, 'relax', true),
                           StateData('vapor', 400.0, 0.004, ...
                           'V', 46.13, 'H', 2738.8, 'S', 9.0035),
                           StateData('vapor', 500.0, 1.0, ...
                           'V', 0.2206, 'H', 2890.2, 'S', 6.8223),
                           StateData('super', 800.0, 0.01, ...
                           'V', 36.92, 'H', 3546.0, 'S', 9.9699),
                           StateData('super', 900.0, 0.70, ...
                           'V', 0.5917, 'H', 3759.4, 'S', 8.2621),
                           StateData('super', 1000.0, 30.0, ...
                           'V', 0.01421, 'H', 3821.6, 'S', 6.6373),
                           StateData('liquid', 500.0, 3.0, ...
                           'D', 832.04, 'H', 975.68, 'S', 2.58049)
                          ];
            self.states = self.states';
            self.refState = StateData('critical', 647.286, 22.089, ...
                                      'D', 317.0, 'H', 2098.8, 'S', 4.4289);

            self.tol = Tolerances('P', 2e-3, 'U', 2e-3, 'S', 2e-3);

            self.setupFluid('water');
            self.runAllTests;
        end

    end

end
