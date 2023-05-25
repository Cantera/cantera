classdef ctMatlabTestThermo < matlab.unittest.TestCase

    properties
        phase
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)
        function ctLoad(testCase)
            % Load Cantera
            ctLoad
            import matlab.unittest.constraints.Throws;
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
        % Generic function to set invalid values to attribute
        function setInvalidValue(testCase, attr, val)
            testCase.phase.(attr) = val;
        end
        % Generic function to get an invalid property
        function a = getInvalidProperty(testCase)
            a = testCase.phase.foobar;
        end
        % Generic function to set an invalid property
        function setInvalidProperty(testCase, val)
            testCase.phase.foobar = val;
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
            testCase.verifyEqual(testCase.phase.phaseName, ...
                                 'spam');

            testCase.verifyGreaterThanOrEqual(testCase.phase.tpID, 0);

            testCase.verifyEqual(testCase.phase.basis, ...
                                 'molar');
            testCase.phase.basis = 'mass';
            testCase.verifyEqual(testCase.phase.basis, ...
                                 'mass');            
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

                testCase.verifyEqual(n{:}, names{i});
                testCase.verifyEqual(i, m);
            end

            testCase.verifyError(@() testCase.phase.speciesName(11), 'Cantera:ctError');
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

                testCase.verifyError(@() testCase.phase.nAtoms('C', 'H2'),...
                                    'Cantera:ctError');
                testCase.verifyError(@() testCase.phase.nAtoms('H', 'CH4'),...
                                    'Cantera:ctError');
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

            diff1 = abs(Zo - exp1)/exp1;
            diff2 = abs(Zh - exp2)/exp2;
            diff3 = abs(Zar - exp3);

            testCase.verifyLessThanOrEqual(diff1, testCase.rtol);
            testCase.verifyLessThanOrEqual(diff2, testCase.rtol);
            testCase.verifyLessThanOrEqual(diff3, testCase.atol);

            testCase.verifyError(@() testCase.phase.elementalMassFraction('C'),...
                                'Cantera:ctError');        
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
                diff = (testWeight - mw(i))/mw(i);

                testCase.verifyLessThanOrEqual(diff, testCase.rtol);
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

        function testSetComposition(testCase)
            xx = zeros(1, testCase.phase.nSpecies);
            xx(3) = 1.0;
            testCase.phase.X = xx;
            yy = testCase.phase.Y;

            testCase.verifyEqual(xx, yy)
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
            testCase.verifyError(@()...
                                testCase.setInvalidValue('X','H2:1.0,CO2:1.5'),...
                                'Cantera:ctError'); 
            testCase.verifyError(@()...
                                testCase.setInvalidValue('X','H2:1.0,O2:asdf'),...
                                'Cantera:ctError'); 
            testCase.verifyError(@()...
                                testCase.setInvalidValue('X','H2:1.e-x4'),...
                                'Cantera:ctError'); 
            testCase.verifyError(@()...
                                testCase.setInvalidValue('X','H2:1e-1.4'),...
                                'Cantera:ctError');  
            testCase.verifyError(@()...
                                testCase.setInvalidValue('X','H2:0.5,O2:1.0,H2:0.1'),...
                                'Cantera:ctError');                 
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

        function checkGetters(testCase)
            val = testCase.phase.T;
            exp = 300;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);
            testCase.verifyGreaterThan(testCase.phase.maxTemp, 0);
            testCase.verifyGreaterThan(testCase.phase.minTemp, 0);

            val = testCase.phase.P;
            exp = 101325;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            val = testCase.phase.D;
            exp = 0.081894;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            val = testCase.phase.V;
            exp = 1/0.081894;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            val = testCase.phase.molarDensity;
            exp = testCase.phase.D/testCase.phase.meanMolecularWeight;
            diff = abs(val - exp)/exp;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            testCase.verifyEqual(testCase.phase.eosType, 'ideal-gas');
            testCase.verifyTrue(testCase.phase.isIdealGas);

            val = testCase.phase.X;
            exp = zeros(1, 10);
            exp(1) = 1.0;
            testCase.verifyEqual(val, exp);

            val = testCase.phase.Y;
            testCase.verifyEqual(val, exp);

            val1 = [testCase.phase.H, testCase.phase.S, ...
                    testCase.phase.U, testCase.phase.G, ...
                    testCase.phase.cp,testCase.phase.cv];
            testCase.phase.basis = 'mass';
            val2 = [testCase.phase.H, testCase.phase.S, ...
                    testCase.phase.U, testCase.phase.G, ...
                    testCase.phase.cp,testCase.phase.cv];
            exp = val2.*testCase.phase.meanMolecularWeight;
            testCase.verifyEqual(val1, exp);

        end

        function checkSetters(testCase)
            
        end

        function testInvalidProperty(testCase)
            testCase.verifyError(@() testCase.getInvalidProperty(),...
                                'MATLAB:noSuchMethodOrField');
            testCase.verifyError(@() testCase.setInvalidProperty(300),...
                                'MATLAB:noPublicFieldForClass');            
        end

        function temperatureSetTest(testCase)
            setPoint = 500;
            testCase.phase.TP = {setPoint, testCase.phase.P};
            val = testCase.phase.T;
            diff = abs(val - setPoint)/setPoint;
            testCase.verifyLessThanOrEqual(diff, testCase.rtol);

            setPoint = -1;
            errMessage = 'Cantera:ctError';

            % testCase.verifyError(@() testCase.setInvalidValue('TP',...
            %                     {setPoint, testCase.phase.P}), ... 
            %                     errMessage);    
            testCase.verifyThat(@() testCase.setInvalidValue('TP',...
                                {setPoint, testCase.phase.P}), ... 
                                Throws(errMessage));   
        end

    end

end