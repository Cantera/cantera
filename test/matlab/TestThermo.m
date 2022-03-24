classdef TestThermo < TestCase
    properties
        gas
    end

    methods
        function self = TestThermo(name)
            self = self@TestCase(name);
        end

        function setUp(self)
            global staticTestThermoGas
            if isempty(staticTestThermoGas)
                staticTestThermoGas = Solution('../data/steam-reforming.yaml', 'full');
            end
            self.gas = staticTestThermoGas;
            set(self.gas, 'T', 300, 'P', oneatm, 'Y', [0.5, 0, 0.5, 0, 0, 0, 0]);
        end

%        function tearDown(self)
%        end

        function testCounts(self)
            assertEqual(nElements(self.gas), 4)
            assertEqual(nSpecies(self.gas), 7)
        end

        function testElements(self)
            for i = 1:nElements(self.gas)
                name = elementName(self.gas, i);
                assertEqual(i, elementIndex(self.gas, name))
            end
        end

        function testSpecies(self)
            for i = 1:nSpecies(self.gas)
                name = speciesName(self.gas, i);
                assertEqual(i, speciesIndex(self.gas, name))
            end
        end

        function test_nAtoms(self)
           assertEqual(nAtoms(self.gas, 1, 3), 4)
           assertEqual(nAtoms(self.gas, 1, 4), 0)
           assertEqual(nAtoms(self.gas, 3, 2), 2)
           assertExceptionThrown(@() nAtoms(self.gas, 2, 5), '')
           assertExceptionThrown(@() nAtoms(self.gas, 8, 1), '')
        end

        function testSetState(self)
            u0 = intEnergy_mass(self.gas);
            h0 = enthalpy_mass(self.gas);
            s0 = entropy_mass(self.gas);
            v0 = 1/density(self.gas);
            T0 = temperature(self.gas);
            P0 = pressure(self.gas);

            set(self.gas, 'T', 400, 'P', 5*oneatm);
            assertAlmostEqual(temperature(self.gas), 400)
            assertAlmostEqual(pressure(self.gas), 5*oneatm)

            set(self.gas, 'H', h0, 'P', P0);
            assertAlmostEqual(temperature(self.gas), T0, 1e-8)
            assertAlmostEqual(entropy_mass(self.gas), s0, 1e-8)

            set(self.gas, 'T', 400, 'P', 5*oneatm);
            set(self.gas, 'U', u0, 'V', v0);
            assertAlmostEqual(pressure(self.gas), P0, 1e-8)
            assertAlmostEqual(enthalpy_mass(self.gas), h0, 1e-8)
        end
    end
end
