classdef ctTestUtility < ctTestCase

    methods (Test)

        function testDataDirectories(self)
            d = ct.dataDirectories();
            self.verifyClass(d, 'cell');
            self.verifyGreaterThan(numel(d), 0);
            for i = 1:numel(d)
                self.verifyClass(d{i}, 'char');
                self.verifyFalse(contains(d{i}, pathsep));
            end
            self.verifyTrue(any(cellfun(@isfolder, d)));
        end

        function testReleaseAfterUnload(self)
            % An object that outlives ct.unload must be released without error.
            % ct.cleanUp only clears the base workspace, so a local variable
            % stands in for any reference it cannot reach.
            self.assumeTrue(string(ct.executionMode) == "outofprocess", ...
                            'Unloading requires out-of-process execution.');

            gas = ct.Solution('h2o2.yaml', '', 'none'); %#ok<NASGU>
            ct.unload;
            ct.load;

            lastwarn('');
            clear gas
            [msg, id] = lastwarn();
            self.verifyNotEqual(id, 'MATLAB:class:DestructorError', msg);
        end

        function testStaleHandleAfterReload(self)
            % An object created before ct.unload must not release an unrelated
            % object that reuses its index after the library is reloaded.
            self.assumeTrue(string(ct.executionMode) == "outofprocess", ...
                            'Unloading requires out-of-process execution.');

            % Start from an empty library so both objects get the same index.
            ct.unload;
            ct.load;
            stale = ct.Solution('h2o2.yaml', '', 'none');
            ct.unload;
            ct.load;
            fresh = ct.Solution('h2o2.yaml', '', 'none');
            self.assertEqual(fresh.solnID, stale.solnID);

            clear stale
            try
                ct.zeroD.IdealGasReactor(fresh);
            catch ME
                self.verifyFail(sprintf(['Releasing a stale object deleted ' ...
                                         'an unrelated one: %s'], ME.message));
            end
        end

    end

end
