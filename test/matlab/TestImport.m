classdef TestImport < TestCase
    methods
        function self = TestImport(name)
            self = self@TestCase(name);
        end

        function testImportXML(self)
            gas = Solution('gri30.xml', 'gri30_mix');
            dkm = mixDiffCoeffs(gas);
            assertEqual(length(dkm), nSpecies(gas))
        end

        function testImportCTI(self)
            gas = Solution('h2o2.cti');
            assertEqual(temperature(gas), 300)
        end
    end
end
