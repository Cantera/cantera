classdef (TestTags = {'Slow'}) ctTestSamples < ctTestCase

    methods (TestMethodSetup)

        function addSampleFolder(self)
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab');
            addpath(samplePath);
        end

    end

    methods (TestMethodTeardown)

        function removeSampleFolder(self)
            close all
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab');
            rmpath(samplePath);
        end

    end

    properties (TestParameter)

            ExampleScript = {
                'equil', 'isentropic', 'reactor1', 'reactor2', 'surf_reactor', ...
                'periodic_cstr', 'lithium_ion_battery', ...
                'rankine', 'prandtl1', 'prandtl2', 'catalytic_combustion', ...
                'burner_flame', 'diffusion_flame', ...
                'diamond_cvd', ...
                % 'ignite_hp', 'ignite_uv', 'plug_flow_reactor', % disabled due to excessive run times; see GH issue #2034
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
