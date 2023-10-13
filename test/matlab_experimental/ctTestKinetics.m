classdef ctTestKinetics < matlab.unittest.TestCase

    properties
        phase
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodSetup)

        function createPhase(self)
            src = 'h2o2.yaml';
            id = 'ohmech';
            transport = 'none';
            self.phase = Solution(src, id, transport);
            self.phase.TPX = {800, 2 * OneAtm, ...
                              [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4, 0]};
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods (Test)

        function testCounts(self)
            self.verifyEqual(self.phase.nReactions, 29);
            self.verifyEqual(self.phase.nTotalSpecies, 10);
            self.verifyEqual(self.phase.nPhases, 1);
        end

        function testIsReversible(self)
            for i = 1:self.phase.nReactions
                self.verifyTrue(self.phase.isReversible(i));
            end

            diamond = Solution('diamond.yaml', 'diamond_100', 'none');
            self.verifyFalse(diamond.isReversible(20));
            clear diamond
        end

        function testMultipler(self)

            f0 = self.phase.forwardRatesOfProgress;
            r0 = self.phase.reverseRatesOfProgress;

            self.phase.setMultiplier(1, 2.0);
            self.phase.setMultiplier(7, 0.1);

            f1 = self.phase.forwardRatesOfProgress;
            r1 = self.phase.reverseRatesOfProgress;

            self.verifyEqual(2 .* f0(1), f1(1), 'AbsTol', self.atol);
            self.verifyEqual(0.1 .* f0(7), f1(7), 'AbsTol', self.atol);
            self.verifyEqual(2 .* r0(1), r1(1), 'AbsTol', self.atol);
            self.verifyEqual(0.1 .* r0(7), r1(7), 'AbsTol', self.atol);

            for i = 1:self.phase.nReactions

                if i ~= 1 && i ~= 7
                    self.verifyEqual(f0(i), f1(i), 'AbsTol', self.atol);
                    self.verifyEqual(r0(i), r1(i), 'AbsTol', self.atol);
                end

            end

            self.phase.setMultiplier(0.5);
            f2 = self.phase.forwardRatesOfProgress;
            r2 = self.phase.reverseRatesOfProgress;
            tol = ones(1, self.phase.nReactions) .* self.atol;
            self.verifyEqual(0.5 .* f0, f2, 'AbsTol', tol);
            self.verifyEqual(0.5 .* r0, r2, 'AbsTol', tol);

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

    end

end
