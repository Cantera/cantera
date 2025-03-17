classdef ctTestTransport < matlab.unittest.TestCase

    properties
        phase
    end

    properties (SetAccess = protected)
        rtol_min = 1e-8;
        rtol_max = 1e-2;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (TestMethodSetup)

        function createPhase(self)
            src = 'h2o2.yaml';
            self.phase = Solution(src);
            self.phase.TPX = {800, 2*OneAtm, ...
                              [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 1e-6, 0.4]};
        end

    end

    methods (TestMethodTeardown)

        function deleteObjects(self)
            props = properties(self);
            for i = 1:length(props)
                prop = self.(props{i});
                if isa(prop, 'handle')
                    delete(prop)
                end
            end
        end

    end

    methods (Test)

        function testScalarProperties(self)
            self.verifyGreaterThan(self.phase.viscosity, 0.0);
            self.verifyGreaterThan(self.phase.thermalConductivity, 0.0);
        end

        function testUnityLewis(self)
            self.phase.transportModel = 'unity-Lewis-number';

            alpha = self.phase.thermalConductivity / (self.phase.D * self.phase.cp);
            Dkm_prime = self.phase.mixDiffCoeffs;

            self.verifyTrue(all(diff(Dkm_prime) < self.atol));
            self.verifyEqual(Dkm_prime(1), alpha, 'AbsTol', self.atol);
        end

        function testMixtureAveraged(self)
            self.verifyEqual(self.phase.transportModel, 'mixture-averaged');

            Dkm1 = self.phase.mixDiffCoeffs;
            Dbin1 = self.phase.binDiffCoeffs;

            self.phase.transportModel = 'multicomponent';
            Dkm2 = self.phase.mixDiffCoeffs;
            Dbin2 = self.phase.binDiffCoeffs;

            self.verifyEqual(Dkm1, Dkm2, 'AbsTol', self.atol);
            self.verifyEqual(Dbin1, Dbin2, 'AbsTol', self.atol);

        end

        function testMixDiffCoeffsChange(self)
            Dkm1 = self.phase.mixDiffCoeffs;
            self.phase.TP = {self.phase.T + 1, self.phase.P};
            Dkm2 = self.phase.mixDiffCoeffs;
            self.verifyTrue(all(Dkm2 > Dkm1));
        end

        function testCKMode(self)
            mu_ct = self.phase.viscosity;
            cond_ct = self.phase.thermalConductivity;
            diff_ct = self.phase.binDiffCoeffs;

            self.phase.transportModel = 'mixture-averaged-CK';
            self.verifyEqual(self.phase.transportModel, 'mixture-averaged-CK');

            mu_ck = self.phase.viscosity;
            cond_ck = self.phase.thermalConductivity;
            diff_ck = self.phase.binDiffCoeffs;

            self.verifyEqual(mu_ct, mu_ck, 'RelTol', self.rtol_max);
            self.verifyNotEqual(mu_ct, mu_ck);
            self.verifyEqual(cond_ct, cond_ck, 'RelTol', self.rtol_max);
            self.verifyNotEqual(cond_ct, cond_ck);
            self.verifyEqual(diff_ct, diff_ck, 'RelTol', self.rtol_max);
            self.verifyNotEqual(diff_ct, diff_ck);
        end

        function testMultiComponent(self)
            try
                a = self.phase.multiDiffCoeffs
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'NotImplementedError');
            end

            self.verifyEqual(self.phase.thermalDiffCoeffs, ...
                             zeros(1, self.phase.nSpecies), 'AbsTol', self.atol);

            self.phase.transportModel = 'multicomponent';
            self.verifyTrue(all(self.phase.multiDiffCoeffs(:) >= 0.0));
            self.verifyTrue(all(self.phase.thermalDiffCoeffs(:) ~= 0.0));
        end

    end

end
