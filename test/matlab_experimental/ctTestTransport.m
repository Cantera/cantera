classdef ctTestTransport < matlab.unittest.TestCase

    properties
        phase
    end

    properties (SetAccess = immutable)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
            copyfile('../data/addReactions_err_test.yaml', ...
                     './addReactions_err_test.yaml');
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            delete('./addReactions_err_test.yaml');
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods
    
        function createPhase(self, src, tr)
            self.phase = Solutoin(src, '', tr);
            self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 1e-6, 0.4];
            self.phase.TP = {800, 2*OneAtm};
        end
    
    end

    methods (Test)

        function testScalarProperties(self)
            self.createPhase('h2o2.yaml', 'default');

            self.verifyGreaterThan(self.phase.viscosity, 0.0);
            self.verifyGreaterThan(self.phase.thermalConductivity, 0.0);
        end

        function testUnityLewis(self)
            self.createPhase('h2o2.yaml', 'unity-lewis-number');

            alpha = self.phase.thermalConductivity / (self.phase.D * self.phase.cp);
            Dkm_prime = phase.mixDiffCoeffs;

            self.verifyTrue(all(diff(Dkm_prime) < self.atol));
            self.verifyEqual(Dkm_prime(1), alpha, 'AbsTol', self.atol);
        end

        function testMixtureAveraged(self)
            self.createPhase('h2o2.yaml')

        end

        function testReactionEquations(self)
            self.verifyEqual(self.phase.nReactions, ...
                             length(self.phase.reactionEquations));
            s = strsplit(self.phase.reactionEquation(19), '<=>');
            r = s{1};
            p = s{2};
            self.verifySubstring(r, 'H');
            self.verifySubstring(r, 'H2O2');
            self.verifySubstring(p, 'HO2');
            self.verifySubstring(p, 'H2');

        end

        function testStoichCoeffs(self)
            nu_r = full(self.phase.reactantStoichCoeffs);
            nu_p = full(self.phase.productStoichCoeffs);

            function checkReactant(s, i, val)
                k = self.phase.kineticsSpeciesIndex(s);
                self.verifyEqual(self.phase.reactantStoichCoeffs(s, i), ...
                                 val);
                self.verifyEqual(self.phase.reactantStoichCoeffs(k, i), ...
                                 val);
                self.verifyEqual(nu_r(k, i), val);
            end

            function checkProduct(s, i, val)
                k = self.phase.kineticsSpeciesIndex(s);
                self.verifyEqual(self.phase.productStoichCoeffs(s, i), ...
                                 val);
                self.verifyEqual(self.phase.productStoichCoeffs(k, i), ...
                                 val);
                self.verifyEqual(nu_p(k, i), val);
            end

            % H + H2O2 <=> HO2 + H2
            checkReactant('H', 19, 1)
            checkReactant('H2O2', 19, 1)
            checkReactant('HO2', 19, 0)
            checkReactant('H2', 19, 0)

            checkProduct('H', 19, 0)
            checkProduct('H2O2', 19, 0)
            checkProduct('HO2', 19, 1)
            checkProduct('H2', 19, 1)

            % 2 O + M <=> O2 + M
            checkReactant('O', 1, 2)
            checkReactant('O2', 1, 0)
            checkProduct('O', 1, 0)
            checkProduct('O2', 1, 1)

        end

        function testRatesOfProgress(self)
            forROP = self.phase.forwardRatesOfProgress;
            revROP = self.phase.reverseRatesOfProgress;
            netROP = self.phase.netRatesOfProgress;

            l = length(netROP);
            tol = ones(1, l) .* self.atol;

            self.verifyEqual(l, self.phase.nReactions);
            self.verifyEqual(forROP - revROP, netROP, 'AbsTol', tol);
        end

        function testRateConstants(self)
            kf = self.phase.forwardRateConstants;
            kr = self.phase.reverseRateConstants;
            Keq = self.phase.equilibriumConstants;

            l = length(kf);
            self.verifyEqual(l, self.phase.nReactions);

            ix = find(kr ~= 0);
            tol = ones(1, l) .* self.rtol;
            self.verifyEqual(kf(ix) ./ kr(ix), Keq(ix), 'RelTol', tol);
        end

        function testSpeciesRates(self)
            nu_p = self.phase.productStoichCoeffs;
            nu_r = self.phase.reactantStoichCoeffs;
            forROP = self.phase.forwardRatesOfProgress;
            revROP = self.phase.reverseRatesOfProgress;
            cr = full(sum(nu_p .* forROP, 2) + sum(nu_r .* revROP, 2))';
            de = full(sum(nu_r .* forROP, 2) + sum(nu_p .* revROP, 2))';

            l = length(self.phase.nSpecies);
            tol = ones(1, l) .* self.rtol;

            self.verifyEqual(self.phase.creationRates, cr, 'RelTol', tol);
            self.verifyEqual(self.phase.destructionRates, de, 'RelTol', tol);
            self.verifyEqual(self.phase.netProdRates, cr - de, 'RelTol', tol);
        end

        function testReactionDeltas(self)
            H = self.phase.deltaEnthalpy;
            S = self.phase.deltaEntropy;
            G = self.phase.deltaGibbs;
            Hs = self.phase.deltaStandardEnthalpy;
            Ss = self.phase.deltaStandardEntropy;
            Gs = self.phase.deltaStandardGibbs;
            T = self.phase.T;

            l = length(H);
            tol = ones(1, l) .* self.rtol;

            self.verifyEqual(H - S .* T, G, 'RelTol', tol);
            self.verifyEqual(Hs - Ss .* T, Gs, 'RelTol', tol);
        end

        function testEmptyKinetics(self)
            try
                gas = Solution('air-no-reactions.yaml');
                arr = zeros(1, gas.nSpecies);
                tol = ones(1, gas.nSpecies) .* self.atol;

                self.verifyEqual(gas.nReactions, 0);
                self.verifyEqual(gas.creationRates, arr, 'AbsTol', tol);
                self.verifyEqual(gas.destructionRates, arr, 'AbsTol', tol);
                self.verifyEqual(gas.netProdRates, arr, 'AbsTol', tol);

                clear gas

            catch ME
                clear gas
                rethrow(ME);
            end
        end

        function testChemicallyActivated(self)
            try
                gas = Solution('chemically-activated-reaction.yaml');

                P = [2026.5, 202650.0, 10132500.0];
                Rf = [2.851022e-04, 2.775924e+00, 2.481792e+03];
                xx = [0.01, 0.01, 0.04, 0.10, 0.84];

                for i = 1:length(P)
                    gas.TPX = {900.0, P(i), xx};
                    self.verifyEqual(gas.forwardRatesOfProgress(1), Rf(i), 'RelTol', 2.0e-05);
                end

                clear gas

            catch ME
                clear gas
                rethrow(ME);
            end
        end

        function testPdepError(self)
            try
                gas = Solution('addReactions_err_test.yaml');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Invalid rate coefficient');
                clear gas
            end
        end

    end

end
