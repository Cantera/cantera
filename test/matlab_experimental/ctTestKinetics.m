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
            self.assumeEqual(self.phase.nReactions, 29);
            self.assumeEqual(self.phase.nTotalSpecies, 10);
            self.assumeEqual(self.phase.nPhases, 1);

            % Missing method
            % self.assumeEqual(self.phase.reactionPhaseIndex, 0);
        end

        function testIsReversible(self)
            % Fails because of incorrect indexing. 
            self.assumeFail();

            for i = 1:self.phase.nReactions
                self.assumeTrue(self.phase.isReversible(i));
            end

        end

        function testMultipler(self)

            f0 = self.phase.forwardRatesOfProgress;
            r0 = self.phase.reverseRatesOfProgress;

            self.phase.setMultiplier(1, 2.0);
            self.phase.setMultiplier(7, 0.1);

            f1 = self.phase.forwardRatesOfProgress;
            r1 = self.phase.reverseRatesOfProgress;

            self.assumeEqual(2 .* f0(1), f1(1), 'AbsTol', self.atol);
            self.assumeEqual(0.1 .* f0(7), f1(7), 'AbsTol', self.atol);
            self.assumeEqual(2 .* r0(1), r1(1), 'AbsTol', self.atol);
            self.assumeEqual(0.1 .* r0(7), r1(7), 'AbsTol', self.atol);

            for i = 1:self.phase.nReactions

                if i ~= 1 && i ~= 7
                    self.assumeEqual(f0(i), f1(i), 'AbsTol', self.atol);
                    self.assumeEqual(r0(i), r1(i), 'AbsTol', self.atol);
                end

            end

            self.phase.setMultiplier(0.5);
            f2 = self.phase.forwardRatesOfProgress;
            r2 = self.phase.reverseRatesOfProgress;
            tol = ones(1, self.phase.nReactions) .* self.atol;
            self.assumeEqual(0.5 .* f0, f2, 'AbsTol', tol);
            self.assumeEqual(0.5 .* r0, r2, 'AbsTol', tol);

        end

        function testReactionEquations(self)
            self.assumeEqual(self.phase.nReactions, ...
                             length(self.phase.reactionEquations));
            s = strsplit(self.phase.reactionEquation(19), '<=>');
            r = s{1};
            p = s{2};
            self.assumeSubstring(r, 'H');
            self.assumeSubstring(r, 'H2O2');
            self.assumeSubstring(p, 'HO2');
            self.assumeSubstring(p, 'H2');

        end

        function testStoichCoeffs(self)
            % Fails because StoichCoeffs methods does not convert species name to 
            % species index. 
            self.assumeFail();

            nu_r = self.phase.reactantStoichCoeffs;
            nu_p = self.phase.productStoichCoeffs;

            function checkReactnat(s, i, val)
                k = self.phase.kineticsSpeciesIndex(s);
                self.assumeEqual(self.phase.reactantStoichCoeffs(s, i), ...
                                 val);
                self.assumeEqual(self.phase.reactantStoichCoeffs(k, i), ...
                                 val);
                self.assumeEqual(nu_r(k, i), val);
            end

            function checkProduct(s, i, val)
                k = self.phase.kineticsSpeciesIndex(s);
                self.assumeEqual(self.phase.productStoichCoeffs(s, i), ...
                                 val);
                self.assumeEqual(self.phase.productStoichCoeffs(k, i), ...
                                 val);
                self.assumeEqual(nu_p(k, i), val);
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
