classdef ctMatlabTestThermo < matlab.unittest.TestCase

    properties
        phase
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassTeardown)
        function testTearDown(testCase)
            % Clean up Cantera
            ctCleanUp
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
        function createPhase(testCase)
            src = 'h2o2.yaml';
            id = 'ohmech';
            transport = 'none';
            testCase.phase = Solution(src, id, transport);
        end

    end

    methods (TestMethodTeardown)
        % Destroy object
        function deleteSolution(testCase)
            clear testCase.phase;
        end
    end

    methods
        % Generic function to check whether a value is near expected value (relative).
        function verifyNearRelative(testCase, val, exp)
            diff = max(abs((val - exp) ./ exp));
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);
        end

        % Generic function to check whether a value is near expected value (absolute).
        function verifyNearAbsolute(testCase, val, exp)
            diff = max(abs(val - exp));
            testCase.verifyLessThanOrEqual(diff, testCase.atol);
        end

        % Generic function to set invalid values to attribute and verify errors
        function setInvalidValue(testCase, attr, val, errMessage)
            try 
                testCase.phase.(attr) = val;
            catch ME
                testCase.verifySubstring(ME.message, errMessage);
            end
        end

        % Generic function to get invalid values of an attribute and verify errors
        function val = getInvalidValue(testCase, attr, args, errMessage)
            try
                if nargin == 3
                    val = testCase.phase.(attr);
                else 
                    val = testCase.phase.(attr)(args{:});
                end
            catch ME
                testCase.verifySubstring(ME.message, errMessage);
            end                
        end
        
        % Generic function to get an invalid property
        function a = getInvalidProperty(testCase)
            a = testCase.phase.foobar;
        end
        
        % Generic function to set an invalid property
        function setInvalidProperty(testCase, val)
            testCase.phase.foobar = val;
        end

        % Check state
        function checkState(testCase, T, D, Y)
            testCase.verifyNearRelative(testCase.phase.T, T);
            testCase.verifyNearRelative(testCase.phase.D, D);            
            testCase.verifyNearAbsolute(testCase.phase.Y, Y);            
        end

        % Check multi properties
        function checkMultiProperties(testCase, str)
            val = testCase.phase.(str);
            for i = 1:length(str)
                attr = str(i);
                exp = testCase.phase.(attr);
                testCase.verifyNearAbsolute(val{i}, exp);
            end
        end

        % Check getter
        function checkGetters(testCase)
            testCase.checkMultiProperties('TD');
            testCase.checkMultiProperties('TDX');
            testCase.checkMultiProperties('TDY');

            testCase.checkMultiProperties('TP');
            testCase.checkMultiProperties('TPX');
            testCase.checkMultiProperties('TPY');

            testCase.checkMultiProperties('HP');
            testCase.checkMultiProperties('HPX');
            testCase.checkMultiProperties('HPY'); 
            
            testCase.checkMultiProperties('UV');
            testCase.checkMultiProperties('UVX');
            testCase.checkMultiProperties('UVY');
            
            testCase.checkMultiProperties('SP');
            testCase.checkMultiProperties('SPX');
            testCase.checkMultiProperties('SPY');
            
            testCase.checkMultiProperties('SV');
            testCase.checkMultiProperties('SVX');
            testCase.checkMultiProperties('SVY'); 

            testCase.checkMultiProperties('DP');
            testCase.checkMultiProperties('DPX');
            testCase.checkMultiProperties('DPY'); 
        end

        % Check setter
        function checkSetters(testCase, T1, D1, Y1)
            val = testCase.phase.TDY;
            T0 = val{1};
            D0 = val{2};
            Y0 = val{3};

            testCase.phase.TDY = {T1, D1, Y1};
            X1 = testCase.phase.X;
            P1 = testCase.phase.P;
            H1 = testCase.phase.H;
            S1 = testCase.phase.S;
            U1 = testCase.phase.U;
            V1 = testCase.phase.V;

            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.TPY = {T1, P1, Y1};
            testCase.checkState(T1, D1, Y1);

            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.UVY = {U1, V1, Y1};
            testCase.checkState(T1, D1, Y1); 
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.HPY = {H1, P1, Y1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.SPY = {S1, P1, Y1};
            testCase.checkState(T1, D1, Y1);

            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.TPX = {T1, P1, X1};
            testCase.checkState(T1, D1, Y1); 
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.UVX = {U1, V1, X1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.HPX = {H1, P1, X1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.SPX = {S1, P1, X1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.SVX = {S1, V1, X1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.SVY = {S1, V1, Y1};
            testCase.checkState(T1, D1, Y1); 

            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.DPX = {D1, P1, X1};
            testCase.checkState(T1, D1, Y1);
            
            testCase.phase.TDY = {T0, D0, Y0};
            testCase.phase.DPY = {D1, P1, Y1};
            testCase.checkState(T1, D1, Y1); 
        end
    end

    methods (Test)
        % Test methods

        function testBaseAttributes(testCase)
            testCase.verifyInstanceOf(testCase.phase.solnName, ...
                                      'char');
            
            testCase.verifyInstanceOf(testCase.phase.phaseName, ...
                                      'char');

            testCase.phase.phaseName = 'spam';
            testCase.verifyMatches(testCase.phase.phaseName, 'spam');

            testCase.verifyGreaterThanOrEqual(testCase.phase.tpID, 0);

            testCase.verifyMatches(testCase.phase.basis, 'molar');
            testCase.phase.basis = 'mass';
            testCase.verifyMatches(testCase.phase.basis, 'mass');            
        end

        function testPhases(testCase)
            % Note: ThermoPhase.phaseOfMatter is not implemented in Clib
            testCase.verifyEqual(testCase.phase.nPhases, 1);
        end

        function testSpecies(testCase)
            testCase.verifyEqual(testCase.phase.nSpecies, 10);

            names = testCase.phase.speciesNames;
            for i = 1:10
                n = testCase.phase.speciesName(i);
                m = testCase.phase.speciesIndex(n{:});

                testCase.verifyMatches(n{:}, names{i});
                testCase.verifyEqual(i, m);
            end

            testCase.getInvalidValue('speciesNames', {11}, 'must not exceed');
        end

        function testElements(testCase)
            % Note: ThermoPhase.elementName is not implemented in Clib
            testCase.verifyEqual(testCase.phase.nElements, 4);
        end

        function testNAtoms(testCase)
            data = {{1, 'O', 'O'}, {2, 'O', 'O2'}, {1, 'H', 'OH'},...
                    {2, 'H', 'H2O'}, {2, 'O', 'H2O2'}, {1, 'Ar', 'AR'},...
                    {0, 'O', 'H'}, {0, 'H', 'AR'}, {0, 'Ar', 'HO2'}};
            for i = 1:length(data)
                n = data{i}{1};
                element = data{i}{2};
                species = data{i}{3};
                mElem = testCase.phase.elementIndex(element);
                kSpec = testCase.phase.speciesIndex(species);
                n1 = testCase.phase.nAtoms(species, element);
                n2 = testCase.phase.nAtoms(kSpec, mElem);

                testCase.verifyEqual(n1, n);
                testCase.verifyEqual(n2, n);

                testCase.getInvalidValue('nAtoms', {'C', 'H2'}, 'no such species');
                testCase.getInvalidValue('nAtoms', {'H', 'CH4'}, 'no such element');
            end
        end

        function testElementalMassFraction(testCase)
            testCase.phase.Y = 'H2O:0.5, O2:0.5';
            Zo = testCase.phase.elementalMassFraction('O');
            Zh = testCase.phase.elementalMassFraction('H');
            Zar = testCase.phase.elementalMassFraction('Ar');

            exp1 = 0.5 + 0.5 * (15.999 / 18.015);
            exp2 = 0.5 * (2.016 / 18.015);
            exp3 = 0.0;

            testCase.verifyNearRelative(Zo, exp1);
            testCase.verifyNearRelative(Zh, exp2);
            testCase.verifyNearAbsolute(Zar, exp3);

            testCase.getInvalidValue('elementalMassFraction', {'C'}, 'No such element');
            testCase.getInvalidValue('elementalMassFraction', {5}, 'No such element');
        end

        function testWeights(testCase)
            aw = testCase.phase.atomicMasses;
            mw = testCase.phase.molecularWeights;

            testCase.verifyEqual(length(aw), testCase.phase.nElements);
            testCase.verifyEqual(length(mw), testCase.phase.nSpecies);

            for i = 1:length(mw)
                testWeight = 0.0;
                for j = 1:length(aw)
                    testWeight = testWeight + ...
                                 aw(j) * testCase.phase.nAtoms(i, j);
                end
                testCase.verifyNearRelative(testWeight, mw(i));
            end
        end

        function testCharges(testCase)
            chargePhase = Solution('gri30_ion.yaml', 'gas');
            charges = chargePhase.charges;
            test = {{'E',-1}, {'N2',0}, {'H3O+',1}};

            for i = 1:length(test)
                species = test{i}{1};
                charge = test{i}{2};

                flag = sum(ismember(chargePhase.speciesNames, species));
                testCase.verifyGreaterThan(flag, 0);
                
                idx = chargePhase.speciesIndex(species);
                testCase.verifyEqual(charges(idx), charge);
            end
            clear chargePhase
        end

        function testReport(testCase)
            str = testCase.phase.report;

            testCase.verifySubstring(str, testCase.phase.phaseName);
            testCase.verifySubstring(str, 'temperature');

            for i = 1:testCase.phase.nSpecies
                name = testCase.phase.speciesName(i);
                testCase.verifySubstring(str, name{:});
            end
        end

        function testRefInfo(testCase)
            testCase.verifyNearRelative(testCase.phase.refPressure, OneAtm);
            testCase.verifyNearRelative(testCase.phase.minTemp, 300);
            testCase.verifyNearRelative(testCase.phase.maxTemp, 3500);
        end

        function testSingleGetters(testCase)
            val = testCase.phase.T;
            exp = 300;
            testCase.verifyNearRelative(val, exp);
            testCase.verifyGreaterThan(testCase.phase.maxTemp, 0);
            testCase.verifyGreaterThan(testCase.phase.minTemp, 0);

            val = testCase.phase.P;
            exp = OneAtm;
            testCase.verifyNearRelative(val, exp);

            val = testCase.phase.D;
            exp = testCase.phase.P * testCase.phase.meanMolecularWeight / ...
                  (GasConstant * testCase.phase.T);
            testCase.verifyNearRelative(val, exp);

            testCase.phase.basis = 'mass';
            val = testCase.phase.V;
            exp = 1/exp;
            testCase.verifyNearRelative(val, exp);
            testCase.phase.basis = 'molar';
            val = testCase.phase.V;
            exp = exp * testCase.phase.meanMolecularWeight;
            testCase.verifyNearRelative(val, exp);

            val = testCase.phase.molarDensity;
            exp = testCase.phase.D/testCase.phase.meanMolecularWeight;
            testCase.verifyNearRelative(val, exp);

            testCase.verifyMatches(testCase.phase.eosType, 'ideal-gas');
            testCase.verifyTrue(testCase.phase.isIdealGas);

            val = testCase.phase.X;
            exp = zeros(1, 10);
            exp(1) = 1.0;
            testCase.verifyNearAbsolute(val, exp);

            val = testCase.phase.Y;
            testCase.verifyNearAbsolute(val, exp);

            val1 = [testCase.phase.H, testCase.phase.S, ...
                    testCase.phase.U, testCase.phase.G, ...
                    testCase.phase.cp, testCase.phase.cv];
            testCase.phase.basis = 'mass';
            val2 = [testCase.phase.H, testCase.phase.S, ...
                    testCase.phase.U, testCase.phase.G, ...
                    testCase.phase.cp, testCase.phase.cv];
            exp = val2.*testCase.phase.meanMolecularWeight;
            testCase.verifyNearRelative(val1, exp);

            val = testCase.phase.isothermalCompressibility;
            exp = 1.0 / testCase.phase.P;
            testCase.verifyNearRelative(val, exp);

            val = testCase.phase.thermalExpansionCoeff;
            exp = 1.0 / testCase.phase.T;
            testCase.verifyNearRelative(val, exp);

        end

        function testGetStateMole(testCase)
            testCase.phase.TDX = {350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'};
            testCase.checkGetters;
        end

        function testGetStateMass(testCase)
            testCase.phase.basis = 'mass';
            testCase.phase.TDY = {350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'};
            testCase.checkGetters;
        end

        function testSetComposition(testCase)
            xx = zeros(1, testCase.phase.nSpecies);
            xx(3) = 1.0;
            testCase.phase.X = xx;
            yy = testCase.phase.Y;

            testCase.verifyNearAbsolute(xx, yy)
        end

        function testSetCompositionBadLength(testCase)
            xx = zeros(1, 5);
            testCase.setInvalidValue('X', [], 'cannot be empty');
            testCase.setInvalidValue('X', xx, 'must be equal');
            testCase.setInvalidValue('Y', xx, 'must be equal');
        end

        function testSetCompositionString(testCase)
            testCase.phase.X = 'h2:1.0, o2:1.0';
            xx = testCase.phase.X;

            diff1 = abs(xx(1) - 0.5)/0.5;
            diff2 = abs(xx(4) - 0.5)/0.5;

            testCase.verifyLessThanOrEqual(diff1, testCase.rtol);
            testCase.verifyLessThanOrEqual(diff2, testCase.rtol);
        end

        function testSetCompositionStringBad(testCase)
            testCase.setInvalidValue('X', 'H2:1.0,CO2:1.5', 'Unknown species');
            testCase.setInvalidValue('X', 'H2:1.0,O2:asdf', 'Trouble processing');
            testCase.setInvalidValue('X', 'H2:1.e-x4', 'Trouble processing');
            testCase.setInvalidValue('X', 'H2:1e-1.4', 'decimal point in exponent');
            testCase.setInvalidValue('X', 'H2:0.5,O2:1.0,H2:0.1', 'Duplicate key');
            testCase.setInvalidValue('X', '', 'cannot be empty');
        end

        function testSetStateMole(testCase)
            testCase.checkSetters(750, 0.07, [0.2, 0.1, 0.0, 0.3, 0.1, ...
                                  0.0, 0.0, 0.2, 0.1, 0.0]);
        end

        function testSetStateMass(testCase)
            testCase.phase.basis = 'mass';
            testCase.checkSetters(500, 1.5, [0.1, 0.0, 0.0, 0.1, 0.4, ...
                                  0.2, 0.0, 0.0, 0.2, 0.0]);            
        end

        function testSetterErrors(testCase)
            testCase.setInvalidValue('TD', 400, 'not supported');
            testCase.setInvalidValue('TP', {300, 101325, 'CH4:1.0'}, ...
                                     'incorrect number');
            testCase.setInvalidValue('HPY', {1.2e6, 101325}, ...
                                     'incorrect number'); 
            testCase.setInvalidValue('UVX', {-4e5, 4.4, 'H2:1.0', -1}, ...
                                     'incorrect number');           
        end

        function testInvalidProperty(testCase)
            testCase.verifyError(@() testCase.getInvalidProperty(),...
                                'MATLAB:noSuchMethodOrField');
            testCase.verifyError(@() testCase.setInvalidProperty(300),...
                                'MATLAB:noPublicFieldForClass');            
        end

    end

end