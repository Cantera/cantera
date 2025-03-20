classdef ctTestSurfaceKinetics < matlab.unittest.TestCase

    properties
        gas
        r1
        r2
        net
        surf1
        surf2
    end

    properties (SetAccess = protected)
        rtol = 1e-6;
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

    methods

        function makeReactors(self)
            self.net = ReactorNet();
            self.interface = Interface('diamond.yaml', 'diamond_100');
            self.gas = self.interface.adjacent('gas');
            self.gas.TPX = {self.gas.T, 1.0e3, 'H:0.002, H2:1, CH4:0.01, CH3:0.0002'};
            self.r1 = IdealGasReactor(self.gas);
            self.r1.V = 0.01;
            self.net.addReactor(self.r1);
            self.r2 = IdealGasReactor(self.gas);
            self.r2.V = 0.01;
            self.net.addReactor(self.r2);
        end

    end

    methods (Test)

        function testCoverages(self)
            self.assumeFail('Skipped until ReactorSurface.coverages is implemented');
        end

        function testCoveragesRegression1(self)
            self.assumeFail('Skipped until ReactorSurface.coverages is implemented');
        end

        function testCoveragesRegression2(self)
            self.assumeFail('Skipped until ReactorSurface.coverages is implemented');
        end

    end

end
