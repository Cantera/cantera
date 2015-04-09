function test_suite = testImport
initTestSuite;

function testImportXML
gas = importPhase('gri30.xml', 'gri30_mix');
dkm = mixDiffCoeffs(gas);
assertEqual(length(dkm), nSpecies(gas))

function testImportCTI
gas = importPhase('h2o2.cti');
assertEqual(temperature(gas), 300)

cleanup()
