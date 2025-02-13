classdef ctTestUndeclared < matlab.unittest.TestCase

    properties
        phase
    end

    methods (TestClassSetup)

        function testSetUp(self)
            copyfile('../data/undeclared-tests.yaml', ...
                     './undeclared-tests.yaml');
            ctTestSetUp
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            delete('./undeclared-tests.yaml');
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods (Test)

        function testRaiseUndeclaredSpecies(self)
            try
                self.phase = Solution('undeclared-tests.yaml', 'A');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'contains undeclared species');
            end
        end

        function testRaiseUndeclaredThirdBodies(self)
            try
                self.phase = Solution('undeclared-tests.yaml', 'B');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'three-body reaction with undeclared');
            end
        end

        function testSkipUndeclaredThirdBodies1(self)
            self.phase = Solution('undeclared-tests.yaml', 'C');
            self.verifyEqual(self.phase.nReactions, 3);
        end

        function testSkipUndeclaredThirdBodies2(self)
            self.phase = Solution('undeclared-tests.yaml', 'D');
            rxns = self.phase.reactionEquations;
            self.verifyTrue(ismember('H + O2 + M <=> HO2 + M', rxns));
        end

        function testSkipUndeclaredOrders(self)
            self.phase = Solution('undeclared-tests.yaml', 'E');
            self.verifyEqual(self.phase.nReactions, 1);
        end

        function testRaiseNonreactantOrders(self)
            try
                self.phase = Solution('undeclared-tests.yaml', 'F');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Reaction order specified');
            end
        end

        function testRaiseUndeclaredOrders(self)
            try
                self.phase = Solution('undeclared-tests.yaml', 'G');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'reaction orders for undeclared');
            end
        end

        function testSkipUndeclaredSurfSpecies(self)
            try
                gas = Solution('undeclared-tests.yaml', 'gas');
                surf = Interface('undeclared-tests.yaml', 'Pt_surf', gas);
                self.verifyEqual(surf.nReactions, 14);

                clear gas
                clear surf

            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                clear gas
                clear surf
            end
        end

    end
end