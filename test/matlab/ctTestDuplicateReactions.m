classdef ctTestDuplicateReactions < ctTestCase

    properties
        phase
    end

    methods

        function check(self, name)
            try
                self.phase = Solution('../data/duplicate-reactions.yaml', name);
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'duplicate reaction');
            end
        end

    end

    methods (Test)

        function testForwardMultiple(self)
            self.check('A');
        end

        function testOpposite1(self)
            self.check('B');
        end

        function testOpposite2(self)
            self.check('C');
        end

        function testOpposite3(self)
            self.check('D');
        end

        function testOpposite4(self)
            self.check('E');
            self.verifyEqual(self.phase.nReactions, 2);
        end

        function testCommonEfficiencies(self)
            self.check('F');
        end

        function testDisjointEfficiencies(self)
            self.check('G');
            self.verifyEqual(self.phase.nReactions, 2);
        end

        function testDifferentType(self)
            self.check('H');
            self.verifyEqual(self.phase.nReactions, 2);
        end

        function testDeclaredDupliate(self)
            self.check('I');
            self.verifyEqual(self.phase.nReactions, 2);
        end

        function testUnmatchedDuplicate(self)
            self.check('J');
        end

        function testNonreactingSpecies(self)
            self.check('K');
            self.verifyEqual(self.phase.nReactions, 3);
        end

    end

end
