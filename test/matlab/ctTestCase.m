classdef ctTestCase < matlab.unittest.TestCase

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

end
