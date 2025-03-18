classdef ctTestFlowReactor1 < matlab.unittest.TestCase

    properties
        gas
        reactor
        net
        surf
        rsurf
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

    methods (Test)

        function testNonReacting(self)
            self.gas = Solution('testConstPressureReactor.yaml');
            self.gas.TPX = {300, OneAtm, 'O2:1.0'};
            self.reactor = FlowReactor(self.gas);
            self.reactor.massFlowRate = 10;

            self.net = ReactorNet();
            self.net.addReactor(self.reactor);

            self.assumeFail('Skipped until Reactor.speed is implemeneted');
            % x = 0;
            % v0 = self.reactor.speed;
            % self.verifyEqual(v0, 10 / r.D, 'RelTol', self.rtol);
            % while x < 10.0
            %     self.net.step();
            %     self.verifyEqual(v0, self.reactor.speed, 'RelTol', self.rtol);
            % end
        end

        function testReacting(self)
            self.gas = Solution('testConstPressureReactor.yaml');
            self.gas.TPX = {1400, 20 * OneAtm, 'CO:1.0, H2O:1.0'};
            self.reactor = FlowReactor(self.gas);
            self.reactor.massFlowRate = 10;

            self.net = ReactorNet();
            self.net.addReactor(self.reactor);

            self.assumeFail('Skipped until Reactor.speed is implemeneted');

            % i = 0
            % while self.net.distance < 1.0
            %     self.net.step();
            %     i = i + 1;
            %     self.verifyEqual(self.reactor.speed * self.reactor.D * self.reactor.area, ...
            %                      10, 'RelTol', self.rtol);
            % end

            % stats = self.net.solverStats;
            % self.verifyEqual(stats('step'), i);

            % x_now = self.net.distance;
            % self.net.advance(x_now);
            % sefl.verifyEqual(self.net.solverStats('steps'), i);
        end

        function testCatalyticSurface(self)
            T0 = 1073.15;
            P0 = OneAtm;
            X0 = 'CH4:1, O2:1.5, AR:0.1';

            self.surf = Interface('methane_pox_on_pt.yaml', 'Pt_surf');
            self.surf.TP = {T0, P0};

            self.assumeFail('Skipped until more Reactor methods are implemented');

            self.gas = self.surf.adjacent('gas');
            self.gas.TPX = {T0, P0, X0};

            % self.reactor = FlowReactor(self.gas);
            % self.reactor.area = 1e-4;
            % porosity = 0.3;
            % velocity = 0.4 / 60;
            % mdot = velocity * self.gas.density * self.reactor.area * porosity;
            % self.reactor.surfaceAreaToVolumeRatio = porosity * 1e5;
            % self.reactor.massFlowRate = mdot;
            % self.reactor.energy = 'off';

            % self.rsurf = ReactorSurface(self.surf, self.reactor);

            % self.net = ReactorNet(self.reactor);
            % kCH4 = self.gas.speciesIndex('CH4');
            % kH2 = self.gas.speciesIndex('H2');
            % kCO = self.gas.speciesIndex('CO');
        end

        function testComponentNames(self)
            self.assumeFail('Skipped until component names are implemented');
        end

    end

end
