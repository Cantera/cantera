import matlab.unittest.TestRunner

thisFile = mfilename('fullpath');
testFolder = fileparts(thisFile);
addpath(genpath(testFolder));
suite = testsuite(testFolder);
% Sample tests are tagged as 'Slow' due to the extra headroom
% required when Cantera is loaded in 'outofprocess' mode.
isSlow = cellfun(@(tags) any(strcmp(tags, 'Slow')), {suite.Tags});
suite = suite(~isSlow);

runner = TestRunner.withTextOutput;

results = runner.run(suite);
disp(results);
