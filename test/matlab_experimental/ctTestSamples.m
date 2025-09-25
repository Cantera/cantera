classdef ctTestSamples < ctTestCase

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
            ExampleScript = {
                'equil', 'isentropic', 'reactor1', 'reactor2', 'surf_reactor', ...
                'periodic_cstr', 'plug_flow_reactor', 'lithium_ion_battery', ...
                'rankine', 'prandtl1', 'prandtl2', 'catcomb','ignite', ...
                'ignite_hp', 'ignite_uv', ...
                % 'flame1', 'flame2', 'diff_flame', 'diamond_cvd'
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
