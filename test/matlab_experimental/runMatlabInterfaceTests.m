import matlab.unittest.TestRunner

setenv('CANTERA_ROOT', pwd);
testFolder = fullfile(pwd, 'test', 'matlab_experimental');
addpath(genpath(testFolder));
suite = testsuite(testFolder);
runner = TestRunner.withTextOutput;

results = runner.run(suite);
disp(results);
