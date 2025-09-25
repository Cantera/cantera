import matlab.unittest.TestRunner

testFolder = fullfile(pwd, 'test', 'matlab_experimental');
addpath(genpath(testFolder));
suite = testsuite(testFolder);
isSlow = cellfun(@(tags) any(strcmp(tags, 'Slow')), {suite.Tags});
suite = suite(~isSlow);

runner = TestRunner.withTextOutput;

results = runner.run(suite);
disp(results);
