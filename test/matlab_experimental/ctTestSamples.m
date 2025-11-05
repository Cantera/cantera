classdef (TestTags = {'Slow'}) ctTestSamples < ctTestCase

    methods (TestMethodSetup)

        function addSampleFolder(self)
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab_experimental');
            addpath(samplePath);
        end

    end

    methods (TestMethodTeardown)

        function removeSampleFolder(self)
            close all
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab_experimental');
            rmpath(samplePath);
        end

    end

    properties (TestParameter)

            % `ignite`, `flame2`, and `diff_flame` will solve without crashing,
            % but produce results that differ from expected. These will be
            % re-enabled when the `Flame` classes are re-implemeneted in MATLAB.
            ExampleScript = {
                'equil', 'isentropic', 'reactor1', 'reactor2', 'surf_reactor', ...
                'periodic_cstr', 'lithium_ion_battery', ...
                'rankine', 'prandtl1', 'prandtl2', 'catcomb', 'diamond_cvd', ...
                'flame1', ...
                % 'ignite_hp', 'ignite_uv', 'plug_flow_reactor', % disabled due to excessive run times; see GH issue #2034
                % 'flame2', 'diff_flame', % broken; fix included in GH PR #2031
                % 'ignite', % disabled as it is broken; see GH issue #2033
            };
    end

    methods (Test, ParameterCombination = 'sequential')

        function testExample(self, ExampleScript)
            try
                feval(ExampleScript);
            catch ME
                self.assertFail(sprintf('Example %s failed: %s', ...
                                ExampleScript, ME.message));
            end
        end

    end

end
