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

    end

end
