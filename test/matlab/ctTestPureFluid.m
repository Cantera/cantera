classdef ctTestPureFluid < ctTestCase

    properties
        fluid
    end

    properties (SetAccess = protected)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestMethodSetup)

        function createPhase(self)
            self.fluid = Water;
        end

    end

    methods

        function checkFdProperties(self, T1, P1, T2, P2, tol)
            self.fluid.TP = {T1, P1};
            self.fluid.basis = 'mass';

            h1a = self.fluid.H;
            cp1 = self.fluid.cp;
            cv1 = self.fluid.cv;
            k1 = self.fluid.isothermalCompressibility;
            alpha1 = self.fluid.thermalExpansionCoeff;
            h1b = self.fluid.H;

            self.fluid.TP = {T2, P2};
            h2a = self.fluid.H;
            cp2 = self.fluid.cp;
            cv2 = self.fluid.cv;
            k2 = self.fluid.isothermalCompressibility;
            alpha2 = self.fluid.thermalExpansionCoeff;
            h2b = self.fluid.H;

            self.verifyEqual(cp1, cp2, 'RelTol', tol);
            self.verifyEqual(cv1, cv2, 'RelTol', tol);
            self.verifyEqual(k1, k2, 'RelTol', tol);
            self.verifyEqual(alpha1, alpha2, 'RelTol', tol);

            self.verifyEqual(h1a, h1b, 'RelTol', 1e-9);
            self.verifyEqual(h2a, h2b, 'RelTol', 1e-9);
        end

    end

    methods (Test)

        function testTemperatureLimits(self)
            CO2 = CarbonDioxide;
            self.verifyEqual(CO2.minTemp, 216.54, 'RelTol', self.rtol);
            self.verifyEqual(CO2.maxTemp, 1500.0, 'RelTol', self.rtol);
        end

        function testCriticalProperties(self)
            self.verifyEqual(self.fluid.critPressure, 22.089e6, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.critTemperature, 647.286, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.critDensity, 317.0, 'RelTol', self.rtol);
        end

        function testSetState(self)
            self.fluid.PQ = {101325, 0.5};
            self.verifyEqual(self.fluid.P, 101325, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.Q, 0.5, 'RelTol', self.rtol);

            self.fluid.TQ = {500, 0.8};
            self.verifyEqual(self.fluid.T, 500, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.Q, 0.8, 'RelTol', self.rtol);
        end

        function testSubstanceSet(self)
            self.fluid.TV = {400, 1.45};
            self.verifyEqual(self.fluid.T, 400, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.V, 1.45, 'RelTol', self.rtol);
            try
                self.fluid.TV = {300, -1};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Negative specific volume');
            end

            self.fluid.PV = {101325, 1.45};
            self.verifyEqual(self.fluid.P, 101325, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.V, 1.45, 'RelTol', self.rtol);

            self.fluid.basis = 'mass';

            self.fluid.UP = {-1.45e7, 101325};
            self.verifyEqual(self.fluid.U, -1.45e7, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.P, 101325, 'RelTol', self.rtol);

            self.fluid.VH = {1.45, -1.45e7};
            self.verifyEqual(self.fluid.V, 1.45, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.H, -1.45e7, 'RelTol', self.rtol);

            self.fluid.TH = {400, -1.45e7};
            self.verifyEqual(self.fluid.T, 400, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.H, -1.45e7, 'RelTol', self.rtol);

            self.fluid.SH = {5000, -1.45e7};
            self.verifyEqual(self.fluid.S, 5000, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.H, -1.45e7, 'RelTol', self.rtol);

            self.fluid.ST = {5000, 400};
            self.verifyEqual(self.fluid.S, 5000, 'RelTol', self.rtol);
            self.verifyEqual(self.fluid.T, 400, 'RelTol', self.rtol);
        end

        function testSetMinMax(self)
            self.fluid.TP = {self.fluid.minTemp, 101325};
            self.verifyEqual(self.fluid.T, self.fluid.minTemp, ...
                             'RelTol', self.rtol);

            self.fluid.TP = {self.fluid.maxTemp, 101325};
            self.verifyEqual(self.fluid.T, self.fluid.maxTemp, ...
                            'RelTol', self.rtol);
        end

        function testPropertiesNearMin(self)
            self.checkFdProperties(self.fluid.minTemp * (1 + 1e-5), 101325, ...
                                   self.fluid.minTemp * (1 + 1e-4), 101325, 1e-2);
        end

        function testPropertiesNearMax(self)
            self.checkFdProperties(self.fluid.maxTemp * (1 - 1e-5), 101325, ...
                                   self.fluid.maxTemp * (1 - 1e-4), 101325, 1e-2);
        end

        function testPropertiesNearSat1(self)
            for T = [340, 390, 420]
                self.fluid.TQ = {T, 0.0};
                p = self.fluid.P;
                self.checkFdProperties(T, p + 0.01, T, p + 0.5, 1e-4);
            end
        end

        function testPropertiesNearSat2(self)
            for T = [340, 390, 420]
                self.fluid.TQ = {T, 0.0};
                p = self.fluid.P;
                self.checkFdProperties(T, p - 0.01, T, p - 0.5, 1e-4);
            end
        end

        function testIsothermalCompressibilityLowP(self)
            ref = Solution('h2o2.yaml', '', 'none');
            ref.TPX = {450, 12, 'H2O:1.0'};
            self.fluid.TP = {450, 12};

            self.verifyEqual(ref.isothermalCompressibility, ...
                             self.fluid.isothermalCompressibility, 'RelTol', 1e-5);
        end

        function testThermalExpansionCoeffLowP(self)
            ref = Solution('h2o2.yaml', '', 'none');
            ref.TPX = {450, 12, 'H2O:1.0'};
            self.fluid.TP = {450, 12};

            self.verifyEqual(ref.thermalExpansionCoeff, ...
                             self.fluid.thermalExpansionCoeff, 'RelTol', 1e-2);
        end

        function testThermalExpansionCoeffTD(self)
            self.fluid.basis = 'mass';
            for T = [440, 550, 660]
                self.fluid.TD = {T, 0.1};
                val = T * self.fluid.thermalExpansionCoeff;
                self.verifyEqual(val, 1.0, 'RelTol', 1e-2);
            end
        end

        function testPQSetterTripleCheck(self)
            self.fluid.PQ = {101325, 0.2};
            t1 = self.fluid.T;
            % Change T such that is would result in a Psat larger than P
            self.fluid.TP = {400, 101325};
            % Ensure that correct triple point pressure is recalculated
            self.fluid.PQ = {101325, 0.2};

            self.verifyEqual(self.fluid.T, t1, 'RelTol', self.rtol);

            try
                % minTemp is triple point temperature
                self.fluid.TP = {self.fluid.minTemp, 101325};
                % triple point pressure
                p1 = self.fluid.satPressure;
                self.fluid.PQ = {0.999 * p1, 0.2};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'below triple point');
            end
        end

        function testQualityExceptions(self)
            % Critical Point
            self.fluid.TP = {300, OneAtm};
            self.fluid.TQ = {self.fluid.critTemperature, 0.5};
            self.verifyEqual(self.fluid.P, self.fluid.critPressure, ...
                             'RelTol', self.rtol);

            self.fluid.TP = {300, OneAtm};
            self.fluid.PQ = {self.fluid.critPressure, 0.5};
            self.verifyEqual(self.fluid.T, self.fluid.critTemperature, ...
                            'RelTol', self.rtol);

            % Supercritical
            try
                self.fluid.TQ = {1.001 * self.fluid.critTemperature, 0.0};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'supercritical');
            end

            try
                self.fluid.PQ = {1.001 * self.fluid.critPressure, 0.0};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'supercritical');
            end

            % Q negative
            try
                self.fluid.TQ = {373.15, -0.001};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Invalid vapor fraction');
            end

            try
                self.fluid.PQ = {OneAtm, -0.001};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Invalid vapor fraction');
            end

            % Q larger than one
            try
                self.fluid.TQ = {373.15, 1.001};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Invalid vapor fraction');
            end

            try
                self.fluid.PQ = {OneAtm, 1.001};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Invalid vapor fraction');
            end
        end

        function testSaturatedMixture(self)
            w = Water;
            self.fluid.TP = {300, OneAtm};

            try
                self.fluid.TP = {300, self.fluid.satPressure};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Saturated mixture detected');
            end

            % Saturated vapor
            self.fluid.TQ = {373.15, 1.0};
            w.TP = {self.fluid.T, 0.999 * self.fluid.satPressure};

            self.verifyEqual(self.fluid.cp, w.cp, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.cv, w.cv, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.thermalExpansionCoeff, ...
                             w.thermalExpansionCoeff, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.isothermalCompressibility, ...
                             w.isothermalCompressibility, 'RelTol', 1e-3);

            % Saturated mixture
            self.fluid.TQ = {373.15, 0.5};

            self.verifyTrue(isinf(self.fluid.cp));
            self.verifyTrue(isnan(self.fluid.cv));
            self.verifyTrue(isinf(self.fluid.thermalExpansionCoeff));
            self.verifyTrue(isinf(self.fluid.isothermalCompressibility));

            % Saturated liquid
            self.fluid.TQ = {373.15, 0.0};
            w.TP = {self.fluid.T, 1.001 * self.fluid.satPressure};

            self.verifyEqual(self.fluid.cp, w.cp, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.cv, w.cv, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.thermalExpansionCoeff, ...
                             w.thermalExpansionCoeff, 'RelTol', 1e-3);
            self.verifyEqual(self.fluid.isothermalCompressibility, ...
                             w.isothermalCompressibility, 'RelTol', 1e-3);
        end

        function testSaturationNearLimits(self)
            % Low temperature limit (triple point)
            self.fluid.TP = {300, OneAtm};
            self.fluid.TP = {self.fluid.minTemp, OneAtm};
            psat = self.fluid.satPressure;

            self.fluid.TP = {300, OneAtm};
            self.fluid.TP = {300, psat};
            self.verifyEqual(self.fluid.satTemperature, self.fluid.minTemp, ...
                             'RelTol', self.rtol);

            % High temperature limit (critical point)
            self.fluid.TP = {300, OneAtm};
            self.fluid.TP = {self.fluid.critTemperature, self.fluid.critPressure};
            self.verifyEqual(self.fluid.satTemperature, self.fluid.critTemperature, ...
                             'RelTol', self.rtol);
            self.verifyEqual(self.fluid.satPressure, self.fluid.critPressure, ...
                             'RelTol', self.rtol);

            % Supercritical
            try
                self.fluid.TP = {1.001 * self.fluid.critTemperature, ...
                                 self.fluid.critPressure};
                p1 = self.fluid.satPressure;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Illegal temperature value');
            end

            try
                self.fluid.TP = {self.fluid.critTemperature, ...
                                 1.001 * self.fluid.critPressure};
                t1 = self.fluid.satTemperature;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Illegal pressure value');
            end

            % Below triple point
            try
                self.fluid.TP = {0.999 * self.fluid.critTemperature, OneAtm};
                p1 = self.fluid.satPressure;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Illegal temperature value');
            end

            % Test disabled pending fix of Github Issue #605
            try
                self.fluid.TP = {300, 0.999 * psat};
                t1 = self.fluid.satTemperature;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                % self.verifySubstring(ME.message, 'Illegal pressure value');
            end
        end

        function testTPQ(self)
            self.assumeFail('Skipping this test until TPQ setter is added');
        end

        function testPhaseOfMatter(self)
            self.assumeFail('Skipping this test until phase of matter is added');
        end

        function testWaterBackends(self)
            self.assumeFail('Skipping this test until thermo model is added');
        end

        function testWaterIPAWS(self)
            w = Water('IAPWS95');
            w.basis = 'mass';

            self.verifyEqual(w.critDensity, 322, 'RelTol', self.rtol);
            self.verifyEqual(w.critTemperature, 647.096, 'RelTol', self.rtol);
            self.verifyEqual(w.critPressure, 22064000, 'RelTol', self.rtol);

            % Test internal TP setters
            w.TP = {300, OneAtm};
            d1 = w.D;
            w.TP = {2000, OneAtm};
            % self.verifyEqual(w.phaseOfMatter, 'supercritical');
            w.TP = {300, OneAtm};
            self.verifyEqual(w.D, d1, 'RelTol', self.rtol);
            % self.verifyEqual(w.phaseOfMatter, 'liquid');

            % Test setters for critical conditions
            w.TP = {w.critTemperature, w.critPressure};
            self.verifyEqual(w.D, 322, 'RelTol', self.rtol);
            w.TP = {2000, OneAtm}; % Uses current density as initial guess
            w.TP = {273.16, OneAtm}; % Uses fixed density as initial guess
            self.verifyEqual(w.D, 999.84376, 'RelTol', self.rtol);
            % self.verifyEqual(w.phaseOfMatter, 'liquid');
            w.TP = {w.T, w.satPressure};
            % self.verifyEqual(w.phaseOfMatter, 'liquid');

            try
                w.TP = {273.15999999, OneAtm};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'assumes liquid phase');
            end

            try
                w.TP = {500, OneAtm};
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'assumes liquid phase');
            end

        end

    end

    properties (TestParameter)
        Tset1 = {273.161, 300.0, 350.0, 400.0, 500.0,
                 600.0, 640.0, 645.0, 646.0, 647.0,
                 647.1, 647.2, 647.22, 647.23, 647.25,
                 647.26, 647.27, 647.28, 647.282, 647.284,
                 647.285, 647.286, 647.287, 650.0, 800.0};
        Pset1 = {1234.0, 101325.0, 5e5, 22.0e6, 22.08e6, 22.09e6, 10001000.0};

        Uset1 = num2cell([0, 100, 200, 500, 1000, 1500, 2000] .* 1000 - 1.58581e7);
        Vset1 = {0.001, 0.002, 0.005, 0.010, 0.10, 0.5, 1.0, 1.5, 2.0};

        Hset1 = num2cell([0, 100, 200, 500, 1000, 1500, 2000] .* 1000 - 1.58581e7);
        Pset2 = {1234.0, 101325.0, 5e5, 22.0e6, 22.08e6, 22.09e6, 10001000.0};
    end

    methods (Test, ParameterCombination = 'exhaustive')

        function testTPConvergence(self, Tset1, Pset1)
            err = '';
            nerr = 0;

            try
                self.fluid.TP = {Tset1, Pset1};
                self.verifyEqual(self.fluid.T, Tset1, 'RelTol', self.rtol);
                self.verifyEqual(self.fluid.P, Pset1, 'RelTol', self.rtol);
            catch e
                err = [err, sprintf('Error at T=%.2f, P=%.2f:\n %s \n\n', ...
                                    Tset1, Pset1, e.message)];
                nerr = nerr + 1;
            end

            if err
                err = [err, sprintf('Total error count: %d', nerr)];
                error(err);
            end
        end

        function testUVConvergence(self, Uset1, Vset1)
            err = '';
            nerr = 0;
            self.fluid.basis = 'mass';

            try
                self.fluid.UV = {Uset1, Vset1};
                self.verifyEqual(self.fluid.U, Uset1, 'RelTol', self.rtol);
                self.verifyEqual(self.fluid.V, Vset1, 'RelTol', self.rtol);
            catch e
                err = [err, sprintf('Error at U=%.2f, V=%.2f:\n %s \n\n', ...
                                    Uset1, Vset1, e.message)];
                nerr = nerr + 1;
            end

            if err
                err = [err, sprintf('Total error count: %d', nerr)];
                error(err);
            end
        end

        function testHPConvergence(self, Hset1, Pset2)
            err = '';
            nerr = 0;
            self.fluid.basis = 'mass';

            try
                self.fluid.HP = {Hset1, Pset2};
                self.verifyEqual(self.fluid.H, Hset1, 'RelTol', self.rtol);
                self.verifyEqual(self.fluid.P, Pset2, 'RelTol', self.rtol);
            catch e
                err = [err, sprintf('Error at U=%.2f, V=%.2f:\n %s \n\n', ...
                                    Hset1, Pset2, e.message)];
                nerr = nerr + 1;
            end

            if err
                err = [err, sprintf('Total error count: %d', nerr)];
                error(err);
            end
        end

    end

end
