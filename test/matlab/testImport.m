function test_suite = testImport
initTestSuite;

function testImportXML
gas = Solution('gri30.xml', 'gri30_mix');
dkm = mixDiffCoeffs(gas);
assertEqual(length(dkm), nSpecies(gas))

function testImportCTI
gas = Solution('h2o2.cti');
assertEqual(temperature(gas), 300)

cleanup()
