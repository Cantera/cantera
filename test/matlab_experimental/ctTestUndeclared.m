classdef ctTestUndeclared < ctTestCase

    properties
        phase
    end

    properties (SetAccess = protected)
        inputfile = '../data/undeclared-tests.yaml';
    end

    methods (Test)

        function testRaiseUndeclaredSpecies(self)
            try
                self.phase = Solution(self.inputfile, 'A');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'contains undeclared species');
            end
        end

        function testRaiseUndeclaredThirdBodies(self)
            try
                self.phase = Solution(self.inputfile, 'B');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'three-body reaction with undeclared');
            end
        end

        function testSkipUndeclaredThirdBodies1(self)
            self.phase = Solution(self.inputfile, 'C');
            self.verifyEqual(self.phase.nReactions, 3);
        end

        function testSkipUndeclaredThirdBodies2(self)
            self.phase = Solution(self.inputfile, 'D');
            rxns = self.phase.reactionEquations;
            self.verifyTrue(ismember('H + O2 + M <=> HO2 + M', rxns));
        end

        function testSkipUndeclaredOrders(self)
            self.phase = Solution(self.inputfile, 'E');
            self.verifyEqual(self.phase.nReactions, 1);
        end

        function testRaiseNonreactantOrders(self)
            try
                self.phase = Solution(self.inputfile, 'F');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Reaction order specified');
            end
        end

        function testRaiseUndeclaredOrders(self)
            try
                self.phase = Solution(self.inputfile, 'G');
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'reaction orders for undeclared');
            end
        end

        function testSkipUndeclaredSurfSpecies(self)
            try
                gas = Solution(self.inputfile, 'gas');
                surf = Interface(self.inputfile, 'Pt_surf', gas);
                self.verifyEqual(surf.nReactions, 14);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
            end
        end

    end
end
