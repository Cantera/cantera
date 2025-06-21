classdef ctTestConstPressureReactor < matlab.unittest.TestCase

    properties
        gas
        gas1
        gas2
        r1
        r2
        solid
        resGas
        env
        mfc1
        mfc2
        interface1
        interface2
        surf1
        surf2
        net1
        net2
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

        function makeReactors(self, arg)
            arguments
                self
                arg.addQ (1,1) logical = false
                arg.addMdot (1,1) logical = false
                arg.addSurf (1,1) logical = false
            end

            src = '../data/testConstPressureReactor.yaml';

            self.gas = Solution(src);
            self.gas.TPX = {900, 25 * OneAtm, 'CO:0.5, H2O:0.2'};

            self.gas1 = Solution(src);
            self.gas2 = Solution(src);

            self.resGas = Solution(src);
            self.solid = Solution('diamond.yaml', 'diamond');

            T0 = 1200;
            P0 = 25 * OneAtm;
            X0 = 'CH4:0.5, H2O:0.2, CO:0.3';

            self.gas1.TPX = {T0, P0, X0};
            self.gas2.TPX = {T0, P0, X0};

            self.r1 = IdealGasReactor(self.gas1);
            self.r2 = Reactor(self.gas2);

            self.r1.V = 0.2;
            self.r2.V = 0.2;

            self.resGas.TP = {T0 - 300, P0};
            self.env = Reservoir(self.resGas);

            U = 300 * arg.addQ;

            if arg.addMdot
                self.mfc1 = MassFlowController(self.env, self.r1);
                self.mfc1.massFlowRate = 0.05;
                self.mfc2 = MassFlowController(self.env, self.r2);
                self.mfc2.massFlowRate = 0.05;
            end

            if arg.addSurf
                self.interface1 = Interface('diamond.yaml', 'diamond_100', ...
                                            self.gas1, self.solid);
                self.interface2 = Interface('diamond.yaml', 'diamond_100', ...
                                            self.gas2, self.solid);

                C = zeros(1, self.interface1.nSpecies);
                C(1) = 0.3;
                C(5) = 0.7;
                C = {C};

                self.surf1 = ReactorSurface(self.interface1, self.r1);
                self.surf1.area = 0.2;
                self.interface1.coverages = C;
                self.surf2 = ReactorSurface(self.interface2, self.r2);
                self.surf2.area = 0.2;
                self.interface2.coverages = C;
            end

            self.net1 = ReactorNet(self.r1);
            self.net2 = ReactorNet(self.r2);
            self.net1.maxTimeStep = 0.05;
            self.net2.maxTimeStep = 0.05;
            % self.net2.maxErrTestFails = 10;
        end

        function assertReactorStatesEqual(self, surf)
            arguments
                self
                surf (1,1) logical = false
            end

            self.verifyEqual(self.r1.contents.Y, self.r2.contents.Y, ...
                            'RelTol', 5e-4, 'AbsTol', 1e-6);
            self.verifyEqual(self.r1.T, self.r2.T, 'RelTol', 5e-5);
            self.verifyEqual(self.r1.P, self.r2.P, 'RelTol', 1e-6);

            if surf
                self.verifyEqual(self.interface1.coverages, ...
                                self.interface2.coverages, ...
                                'RelTol', 1e-4, 'AbsTol', 1e-8);
            end
        end


        function integrate(self, surf)
            arguments
                self
                surf (1,1) logical = false
            end

            for t = 0.5:1.0:50
                self.net1.advance(t);
                self.net2.advance(t);
                self.assertReactorStatesEqual(surf)
            end
        end

    end

    methods (Test)

        function testReactorSurfaceType(self)
            self.makeReactors('addSurf', true);
            self.verifyTrue(isa(self.surf1, 'ReactorSurface'));
            self.verifyTrue(startsWith(self.surf1.name, 'ReactorSurface_'));
            self.surf1.name = 'name-of-reactor-surface';
            self.verifyEqual(self.surf1.name, 'name-of-reactor-surface');
        end

        function testClosed(self)
            self.makeReactors();
            self.integrate();
        end

        function testWithHeatTransfer(self)
            self.makeReactors('addQ', true);
            self.integrate();
        end

        function testWithMdot(self)
            self.makeReactors('addMdot', true);
            self.integrate();
        end

        function testWithSurfaceReactions(self)
            self.makeReactors('addSurf', true);
            self.net1.atol = 1e-18;
            self.net2.atol = 1e-18;
            self.net1.rtol = 1e-9;
            self.net2.rtol = 1e-9;
            self.integrate();
        end

    end

end
