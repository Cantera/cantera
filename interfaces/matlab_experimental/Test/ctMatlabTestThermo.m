classdef ctMatlabTestThermo < matlab.unittest.TestCase

    properties
        phase
        rtol
    end

    methods (TestClassSetup)
        function ctLoad(testCase)
            % Load Cantera
            ctLoad
            testCase.rtol = 1e-6;
        end
    end

    methods (TestClassTeardown)
        function ctUnload(testCase)
            % Unload Cantera
            ctUnload
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
        function createSolution(testCase)
            src = 'gri30.yaml';
            id = 'gri30';
            trans = 'mixture-averaged';
            testCase.phase = Solution(src, id, trans);
        end

    end

    methods (TestMethodTeardown)
        % Destroy object
        function deleteSolution(testCase)
            clear testCase.phase;
        end
    end

    methods (Test)
        % Test methods

        function testSource(testCase)
            
        end

        function temperatureTest(testCase)
            val = testCase.phase.T;
            exp = 300;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);
        end

        function pressureTest(testCase)
            val = testCase.phase.P;
            exp = 101325;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);
        end

        function temperatureSetTest(testCase)
            setPoint = 500;
            testCase.phase.TP = {setPoint, testCase.phase.P};
            val = testCase.phase.T;
            diff = abs(val - setPoint)/setPoint;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            setPoint = -1;
            errMessage = ?MException;
            function setTnegative(sp)
                testCase.phase.TP = {sp, testCase.phase.P};
            end

            testCase.verifyError(@() setTnegative(setPoint), ... 
                                errMessage);        
    end

    end

end