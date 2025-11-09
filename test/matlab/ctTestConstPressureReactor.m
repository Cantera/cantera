classdef ctTestConstPressureReactor < ctTestCase
    % The constant pressure reactor should give essentially the same results as
    % as a regular "Reactor" with a wall with a very high expansion rate
    % coefficient. Replicates the tests done in the Python test suite in
    % class `TestConstPressureReactor`.

    properties
        gas
        gas1
        gas2
        r1
        r2
        solid
        resGas
        env
        w1
        w2
        mfc1
        mfc2
        interface1
        interface2
        surf1
        surf2
        net1
        net2
    end

    methods

        function makeReactors(self, arg)
            arguments
                self
                arg.addQ (1,1) logical = false
                arg.addMdot (1,1) logical = false
                arg.addSurf (1,1) logical = false
            end

            src = '../data/ch4_minimal.yaml';

            self.gas = ct.Solution(src, 'testConstPressureReactor');
            self.gas.TPX = {900, 25 * ct.OneAtm, 'CO:0.5, H2O:0.2'};

            self.gas1 = ct.Solution(src, 'testConstPressureReactor');
            self.gas2 = ct.Solution(src, 'testConstPressureReactor');

            self.resGas = ct.Solution(src, 'testConstPressureReactor');
            self.solid = ct.Solution('diamond.yaml', 'diamond');

            T0 = 1200;
            P0 = 25 * ct.OneAtm;
            X0 = 'CH4:0.5, H2O:0.2, CO:0.3';

            self.gas1.TPX = {T0, P0, X0};
            self.gas2.TPX = {T0, P0, X0};

            self.r1 = ct.zeroD.IdealGasReactor(self.gas1);
            self.r2 = ct.zeroD.ConstPressureReactor(self.gas2);

            self.r1.V = 0.2;
            self.r2.V = 0.2;

            self.resGas.TP = {T0 - 300, P0};
            self.env = ct.zeroD.Reservoir(self.resGas);

            U = 300 * arg.addQ;

            self.w1 = ct.zeroD.Wall(self.r1, self.env);
            self.w1.area = 0.1;
            self.w1.expansionRateCoeff = 1e3;
            self.w1.heatTransferCoeff = U;

            self.w2 = ct.zeroD.Wall(self.r2, self.env);
            self.w2.area = 0.1;
            self.w2.heatTransferCoeff = U;

            if arg.addMdot
                self.mfc1 = ct.zeroD.MassFlowController(self.env, self.r1);
                self.mfc1.massFlowRate = 0.05;
                self.mfc2 = ct.zeroD.MassFlowController(self.env, self.r2);
                self.mfc2.massFlowRate = 0.05;
            end

            if arg.addSurf
                self.interface1 = ct.Interface('diamond.yaml', 'diamond_100', ...
                                            self.gas1, self.solid);
                self.interface2 = ct.Interface('diamond.yaml', 'diamond_100', ...
                                            self.gas2, self.solid);

                C = zeros(1, self.interface1.nSpecies);
                C(1) = 0.3;
                C(5) = 0.7;

                self.surf1 = ct.zeroD.ReactorSurface(self.interface1, self.r1);
                self.surf1.area = 0.2;
                self.interface1.coverages = C;
                self.surf2 = ct.zeroD.ReactorSurface(self.interface2, self.r2);
                self.surf2.area = 0.2;
                self.interface2.coverages = C;
            end

            self.net1 = ct.zeroD.ReactorNet(self.r1);
            self.net2 = ct.zeroD.ReactorNet(self.r2);
            self.net1.maxTimeStep = 0.05;
            self.net2.maxTimeStep = 0.05;
        end

        function assertReactorStatesEqual(self, surf)
            arguments
                self
                surf (1,1) logical = false
            end

            self.verifyEqual(self.r1.phase.Y, self.r2.phase.Y, ...
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
            self.verifyTrue(isa(self.surf1, 'ct.zeroD.ReactorSurface'));
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
