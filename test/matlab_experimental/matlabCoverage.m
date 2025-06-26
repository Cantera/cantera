% generateCodeCoverageReport.m
% Generate a Cobertura-style code coverage report for the experimental MATLAB interface.

import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat

sourceFolder = fullfile(pwd, 'interfaces', 'matlab_experimental');
testFolder = fullfile(pwd, 'test', 'matlab_experimental');

suite = testsuite(testFolder);
runner = TestRunner.withTextOutput;

reportFormat = CoberturaFormat(fullfile(testFolder, 'matlabCoverage.xml'));
coveragePlugin = CodeCoveragePlugin.forFolder(sourceFolder, ...
                                              'IncludingSubfolders', true, ...
                                              'Producing', reportFormat);
runner.addPlugin(coveragePlugin);

results = runner.run(suite);
disp(results);
