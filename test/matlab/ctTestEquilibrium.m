classdef ctTestEquilibrium < ctTestCase

    properties
        phase
        mix
    end

    properties (SetAccess = protected)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods
        function checkval(self, names, moles)
            ntot = sum(moles);
            xx = moles ./ ntot;
            for i = 1:length(names)
                val = self.phase.X(self.phase.speciesIndex(names{i}));
                self.verifyEqual(val, xx(i), 'AbsTol', self.atol);
            end
        end
    end

    methods (Test)

        function testEquilCompleteStoichiometric(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'complete');
            self.phase.TPX = {298, 1.0e6, 'CH4:1.0, O2:2.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 0, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilCompleteLean(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'complete');
            self.phase.TPX = {298, 1.0e6, 'CH4:1.0, O2:3.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 1, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilIncompleteStoichiometric(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'incomplete');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:2.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 0, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilIncompleteLean(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'incomplete');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:3.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 1, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilGriStoichiometric(self)
            self.phase = ct.Solution('gri30.yaml', '', 'none');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:2.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 0, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilGriLean(self)
            self.phase = ct.Solution('gri30.yaml', '', 'none');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:3.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 1, 2, 1];
            self.checkval(names, moles);
        end

        function testEquilOverconstrained1(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'overconstrained-1');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:1.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2'};
            moles = [1, 1];
            self.checkval(names, moles);
        end

        function testEquilOverconstrained2(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', 'overconstrained-2');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:1.0'};
            self.phase.equilibrate('TP');

            names = {'CH4', 'O2'};
            moles = [1, 1];
            self.checkval(names, moles);
        end

        function testEquilGriStoichiometricGibbs(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', '', 'none');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:2.0'};
            self.phase.equilibrate('TP', 'gibbs');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 0, 2, 1];

            % self.assumeFail('Skipping multi-phase equilibrium test');
            self.checkval(names, moles);
        end

        function testEquilGriLeanGibbs(self)
            self.phase = ct.Solution('../data/equilibrium.yaml', '', 'none');
            self.phase.TPX = {301, 1.0e6, 'CH4:1.0, O2:3.0'};
            self.phase.equilibrate('TP', 'gibbs');

            names = {'CH4', 'O2', 'H2O', 'CO2'};
            moles = [0, 1, 2, 1];

            % self.assumeFail('Skipping multi-phase equilibrium test');
            self.checkval(names, moles);
        end

        function testKOHEquilTP(self)
            phasenames = {'K_solid', 'K_liquid', ...
                          'KOH_a', 'KOH_b', 'KOH_liquid', ...
                          'K2O2_solid', 'K2O_solid', 'KO2_solid', ...
                          'ice', 'liquid_water', 'KOH_plasma'};
            phases = ImportPhases('KOH.yaml', phasenames);
            self.mix = ct.Mixture(phases);

            temp = 350:300:5000;
            data = zeros(length(temp), self.mix.nSpecies + 1);
            data(:, 1) = temp;

            for i  = 1:length(temp)
                self.mix.T = temp(i);
                self.mix.P = ct.OneAtm;
                self.mix.setSpeciesMoles('K:1.03, H2:2.12, O2:0.9');
                self.mix.equilibrate('TP');

                data(i, 2:end) = self.mix.speciesMoles;
            end

            refData = readmatrix('../data/koh-equil-TP.csv');
            self.verifySize(data, size(refData), ...
                            'Generated data does not match reference size');
            self.verifyEqual(data, refData, 'AbsTol', self.atol);
        end

        function testKOHEquilHP(self)
            phasenames = {'K_solid', 'K_liquid', ...
                          'KOH_a', 'KOH_b', 'KOH_liquid', ...
                          'K2O2_solid', 'K2O_solid', 'KO2_solid', ...
                          'ice', 'liquid_water', 'KOH_plasma'};
            phases = ImportPhases('KOH.yaml', phasenames);
            self.mix = ct.Mixture(phases);

            temp = 350:300:5000;
            dT = 1;
            data = zeros(length(temp), self.mix.nSpecies + 2);
            data(:, 1) = temp;

            self.mix.P = ct.OneAtm;
            for i  = 1:length(temp)
                self.mix.setSpeciesMoles('K:1.03, H2:2.12, O2:0.9');
                self.mix.T = temp(i) - dT;
                self.mix.equilibrate('TP');
                self.mix.T = temp(i);
                self.mix.equilibrate('HP');

                data(i, 2) = self.mix.T;
                data(i, 3:end) = self.mix.speciesMoles;
            end

            refData = readmatrix('../data/koh-equil-HP.csv');
            self.verifySize(data, size(refData), ...
                            'Generated data does not match reference size');
            self.verifyEqual(data, refData, 'AbsTol', self.atol);
        end

        function testIdealSolidSolnPhaseEquil(self)
            self.phase = ct.Solution('../data/IdealSolidSolnPhaseExample.yaml');
            self.phase.TPX = {500, ct.OneAtm, 'C2H2-graph: 1.0'};
            self.phase.equilibrate('TP', 'element_potential');

            self.verifyEqual(self.phase.X(self.phase.speciesIndex('C-graph')), ...
                             2.0/3.0, 'AbsTol', self.atol);
            self.verifyEqual(self.phase.X(self.phase.speciesIndex('H2-solute')), ...
                             1.0/3.0, 'AbsTol', self.atol);
        end

    end

end
