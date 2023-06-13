classdef ctTestThermo < matlab.unittest.TestCase

    properties
        phase
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
            % Clean up Cantera
            ctCleanUp
            ctTestTearDown
        end
    end

    methods (TestMethodSetup)

        function createPhase(self)
            src = 'h2o2.yaml';
            id = 'ohmech';
            transport = 'none';
            self.phase = Solution(src, id, transport);
        end

    end

    methods (TestMethodTeardown)

        function deleteSolution(self)
            clear self.phase;
        end

    end

    methods

        % Generic function to set invalid values to attribute and assume errors
        function setInvalidValue(self, attr, val, errMessage)
            try 
                self.phase.(attr) = val;
            catch ME
                self.assumeSubstring(ME.message, errMessage);
            end
        end

        % Generic function to get invalid values of an attribute and assume errors
        function val = getInvalidValue(self, attr, args, errMessage)
            try
                if nargin == 3
                    val = self.phase.(attr);
                else 
                    val = self.phase.(attr)(args{:});
                end
            catch ME
                self.assumeSubstring(ME.message, errMessage);
            end                
        end
        
        % Generic function to get an invalid property
        function a = getInvalidProperty(self)
            a = self.phase.foobar;
        end
        
        % Generic function to set an invalid property
        function setInvalidProperty(self, val)
            self.phase.foobar = val;
        end

        % Check state
        function checkState(self, T, D, Y)
            self.assumeEqual(self.phase.T, T, 'RelTol', self.rtol);
            self.assumeEqual(self.phase.D, D, 'RelTol', self.rtol);            
            self.assumeEqual(self.phase.Y, Y, 'AbsTol', self.atol);            
        end

        % Check multi properties
        function checkMultiProperties(self, str)
            val = self.phase.(str);
            for i = 1:length(str)
                attr = str(i);
                exp = self.phase.(attr);
                self.assumeEqual(val{i}, exp, 'RelTol', self.rtol);
            end
        end

        % Check getter
        function checkGetters(self)
            self.checkMultiProperties('TD');
            self.checkMultiProperties('TDX');
            self.checkMultiProperties('TDY');

            self.checkMultiProperties('TP');
            self.checkMultiProperties('TPX');
            self.checkMultiProperties('TPY');

            self.checkMultiProperties('HP');
            self.checkMultiProperties('HPX');
            self.checkMultiProperties('HPY'); 
            
            self.checkMultiProperties('UV');
            self.checkMultiProperties('UVX');
            self.checkMultiProperties('UVY');
            
            self.checkMultiProperties('SP');
            self.checkMultiProperties('SPX');
            self.checkMultiProperties('SPY');
            
            self.checkMultiProperties('SV');
            self.checkMultiProperties('SVX');
            self.checkMultiProperties('SVY'); 

            self.checkMultiProperties('DP');
            self.checkMultiProperties('DPX');
            self.checkMultiProperties('DPY'); 
        end

        % Check setter
        function checkSetters(self, T1, D1, Y1)
            val = self.phase.TDY;
            T0 = val{1};
            D0 = val{2};
            Y0 = val{3};

            self.phase.TDY = {T1, D1, Y1};
            X1 = self.phase.X;
            P1 = self.phase.P;
            H1 = self.phase.H;
            S1 = self.phase.S;
            U1 = self.phase.U;
            V1 = self.phase.V;

            self.phase.TDY = {T0, D0, Y0};
            self.phase.TPY = {T1, P1, Y1};
            self.checkState(T1, D1, Y1);

            self.phase.TDY = {T0, D0, Y0};
            self.phase.UVY = {U1, V1, Y1};
            self.checkState(T1, D1, Y1); 
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.HPY = {H1, P1, Y1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.SPY = {S1, P1, Y1};
            self.checkState(T1, D1, Y1);

            self.phase.TDY = {T0, D0, Y0};
            self.phase.TPX = {T1, P1, X1};
            self.checkState(T1, D1, Y1); 
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.UVX = {U1, V1, X1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.HPX = {H1, P1, X1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.SPX = {S1, P1, X1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.SVX = {S1, V1, X1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.SVY = {S1, V1, Y1};
            self.checkState(T1, D1, Y1); 

            self.phase.TDY = {T0, D0, Y0};
            self.phase.DPX = {D1, P1, X1};
            self.checkState(T1, D1, Y1);
            
            self.phase.TDY = {T0, D0, Y0};
            self.phase.DPY = {D1, P1, Y1};
            self.checkState(T1, D1, Y1); 
        end
    end

    methods (Test)
        % Test methods

        function testBaseAttributes(self)
            self.assumeInstanceOf(self.phase.solnName, ...
                                      'char');
            
            % self.assumeInstanceOf(self.phase.phaseName, ...
            %                           'char');

            % self.phase.phaseName = 'spam';
            % self.assumeMatches(self.phase.phaseName, 'spam');

            self.assumeGreaterThanOrEqual(self.phase.tpID, 0);

            self.assumeMatches(self.phase.basis, 'molar');
            self.phase.basis = 'mass';
            self.assumeMatches(self.phase.basis, 'mass');            
        end

        function testPhases(self)
            % Note: ThermoPhase.phaseOfMatter is not implemented in Clib
            self.assumeEqual(self.phase.nPhases, 1);
        end

        function testSpecies(self)
            self.assumeEqual(self.phase.nSpecies, 10);

            names = self.phase.speciesNames;
            for i = 1:10
                n = self.phase.speciesName(i);
                m = self.phase.speciesIndex(n{:});

                self.assumeMatches(n{:}, names{i});
                self.assumeEqual(i, m);
            end

            self.getInvalidValue('speciesNames', {11}, 'must not exceed');
        end

        function testElements(self)
            % Note: ThermoPhase.elementName is not implemented in Clib
            self.assumeEqual(self.phase.nElements, 4);
        end

        function testNAtoms(self)
            data = {{1, 'O', 'O'}, {2, 'O', 'O2'}, {1, 'H', 'OH'},...
                    {2, 'H', 'H2O'}, {2, 'O', 'H2O2'}, {1, 'Ar', 'AR'},...
                    {0, 'O', 'H'}, {0, 'H', 'AR'}, {0, 'Ar', 'HO2'}};
            for i = 1:length(data)
                n = data{i}{1};
                element = data{i}{2};
                species = data{i}{3};
                mElem = self.phase.elementIndex(element);
                kSpec = self.phase.speciesIndex(species);
                n1 = self.phase.nAtoms(species, element);
                n2 = self.phase.nAtoms(kSpec, mElem);

                self.assumeEqual(n1, n);
                self.assumeEqual(n2, n);

                self.getInvalidValue('nAtoms', {'C', 'H2'}, 'no such species');
                self.getInvalidValue('nAtoms', {'H', 'CH4'}, 'no such element');
            end
        end

        function testElementalMassFraction(self)
            self.phase.Y = 'H2O:0.5, O2:0.5';
            Zo = self.phase.elementalMassFraction('O');
            Zh = self.phase.elementalMassFraction('H');
            Zar = self.phase.elementalMassFraction('Ar');

            exp1 = 0.5 + 0.5 * (15.999 / 18.015);
            exp2 = 0.5 * (2.016 / 18.015);
            exp3 = 0.0;

            self.assumeEqual(Zo, exp1, 'AbsTol', self.atol);
            self.assumeEqual(Zh, exp2, 'AbsTol', self.atol);
            self.assumeEqual(Zar, exp3, 'AbsTol', self.atol);

            self.getInvalidValue('elementalMassFraction', {'C'}, 'No such element');
            self.getInvalidValue('elementalMassFraction', {5}, 'No such element');
        end

        function testWeights(self)
            aw = self.phase.atomicMasses;
            mw = self.phase.molecularWeights;

            self.assumeEqual(length(aw), self.phase.nElements);
            self.assumeEqual(length(mw), self.phase.nSpecies);

            for i = 1:length(mw)
                testWeight = 0.0;
                for j = 1:length(aw)
                    testWeight = testWeight + ...
                                 aw(j) * self.phase.nAtoms(i, j);
                end
                self.assumeEqual(testWeight, mw(i), 'RelTol', self.rtol);
            end
        end

        function testCharges(self)
            chargePhase = Solution('gri30_ion.yaml', 'gas');
            charges = chargePhase.charges;
            test = {{'E',-1}, {'N2',0}, {'H3O+',1}};

            for i = 1:length(test)
                species = test{i}{1};
                charge = test{i}{2};

                flag = sum(ismember(chargePhase.speciesNames, species));
                self.assumeGreaterThan(flag, 0);
                
                idx = chargePhase.speciesIndex(species);
                self.assumeEqual(charges(idx), charge);
            end
            clear chargePhase
        end

        % function testReport(self)
        %     str = self.phase.report;

        %     self.assumeSubstring(str, self.phase.phaseName);
        %     self.assumeSubstring(str, 'temperature');

        %     for i = 1:self.phase.nSpecies
        %         name = self.phase.speciesName(i);
        %         self.assumeSubstring(str, name{:});
        %     end
        % end

        % function testRefInfo(self)
        %     self.assumeEqual(self.phase.refPressure, OneAtm, 'RelTol', self.rtol);
        %     self.assumeEqual(self.phase.minTemp, 300, 'RelTol', self.rtol);
        %     self.assumeEqual(self.phase.maxTemp, 3500, 'RelTol', self.rtol);
        % end

        function testSingleGetters(self)
            val = self.phase.T;
            exp = 300;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);
            self.assumeGreaterThan(self.phase.maxTemp, 0);
            self.assumeGreaterThan(self.phase.minTemp, 0);

            val = self.phase.P;
            exp = OneAtm;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

            val = self.phase.D;
            exp = self.phase.P * self.phase.meanMolecularWeight / ...
                  (GasConstant * self.phase.T);
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

            self.phase.basis = 'mass';
            val = self.phase.V;
            exp = 1/exp;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);
            self.phase.basis = 'molar';
            val = self.phase.V;
            exp = exp * self.phase.meanMolecularWeight;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

            val = self.phase.molarDensity;
            exp = self.phase.D/self.phase.meanMolecularWeight;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

            self.assumeMatches(self.phase.eosType, 'ideal-gas');
            self.assumeTrue(self.phase.isIdealGas);

            val = self.phase.X;
            exp = zeros(1, 10);
            exp(1) = 1.0;
            tol = ones(1, 10).*self.rtol;
            self.assumeEqual(val, exp, 'RelTol', tol);

            val = self.phase.Y;
            self.assumeEqual(val, exp, 'RelTol', tol);

            val1 = [self.phase.H, self.phase.S, ...
                    self.phase.U, self.phase.G, ...
                    self.phase.cp, self.phase.cv];
            self.phase.basis = 'mass';
            val2 = [self.phase.H, self.phase.S, ...
                    self.phase.U, self.phase.G, ...
                    self.phase.cp, self.phase.cv];
            exp = val2.*self.phase.meanMolecularWeight;
            tol = ones(1, 9).*self.rtol;
            self.assumeEqual(val1, exp, 'RelTol', tol);

            val = self.phase.isothermalCompressibility;
            exp = 1.0 / self.phase.P;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

            val = self.phase.thermalExpansionCoeff;
            exp = 1.0 / self.phase.T;
            self.assumeEqual(val, exp, 'RelTol', self.rtol);

        end

        function testGetStateMole(self)
            self.phase.TDX = {350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'};
            self.checkGetters;
        end

        function testGetStateMass(self)
            self.phase.basis = 'mass';
            self.phase.TDY = {350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'};
            self.checkGetters;
        end

        function testSetComposition(self)
            xx = zeros(1, self.phase.nSpecies);
            xx(3) = 1.0;
            self.phase.X = xx;
            yy = self.phase.Y;
            tol = ones(1, 10).*self.atol;

            self.assumeEqual(xx, yy, 'AbsTol', tol)
        end

        function testSetCompositionBadLength(self)
            xx = zeros(1, 5);
            self.setInvalidValue('X', [], 'cannot be empty');
            self.setInvalidValue('X', xx, 'must be equal');
            self.setInvalidValue('Y', xx, 'must be equal');
        end

        function testSetCompositionString(self)
            self.phase.X = 'h2:1.0, o2:1.0';
            xx = self.phase.X;

            self.assumeEqual(xx(1), 0.5, 'AbsTol', self.atol);
            self.assumeEqual(xx(4), 0.5, 'AbsTol', self.atol);
        end

        function testSetCompositionStringBad(self)
            self.setInvalidValue('X', 'H2:1.0,CO2:1.5', 'Unknown species');
            self.setInvalidValue('X', 'H2:1.0,O2:asdf', 'Trouble processing');
            self.setInvalidValue('X', 'H2:1.e-x4', 'Trouble processing');
            self.setInvalidValue('X', 'H2:1e-1.4', 'decimal point in exponent');
            self.setInvalidValue('X', 'H2:0.5,O2:1.0,H2:0.1', 'Duplicate key');
            self.setInvalidValue('X', '', 'cannot be empty');
        end

        function testSetStateMole(self)
            self.checkSetters(750, 0.07, [0.2, 0.1, 0.0, 0.3, 0.1, ...
                                  0.0, 0.0, 0.2, 0.1, 0.0]);
        end

        % function testSetStateMass(self)
        %     self.phase.basis = 'mass';
        %     self.checkSetters(500, 1.5, [0.1, 0.0, 0.0, 0.1, 0.4, ...
        %                           0.2, 0.0, 0.0, 0.2, 0.0]);            
        % end

        function testSetterErrors(self)
            self.setInvalidValue('TD', 400, 'not supported');
            self.setInvalidValue('TP', {300, 101325, 'CH4:1.0'}, ...
                                     'incorrect number');
            self.setInvalidValue('HPY', {1.2e6, 101325}, ...
                                     'incorrect number'); 
            self.setInvalidValue('UVX', {-4e5, 4.4, 'H2:1.0', -1}, ...
                                     'incorrect number');           
        end

        function testInvalidProperty(self)
            self.assumeError(@() self.getInvalidProperty(),...
                                'MATLAB:noSuchMethodOrField');
            self.assumeError(@() self.setInvalidProperty(300),...
                                'MATLAB:noPublicFieldForClass');            
        end

    end

end
