// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera.Tests;

public class ThermoPhaseTest
{
    [Fact]
    public void ThermoPhase_SpeciesRetrieved()
    {
        using var thermo = Application.CreateThermoPhase("gri30.yaml");

        Assert.NotEmpty(thermo.Species);
        Assert.NotEmpty(thermo.Species[0].Name);
        Assert.Equal(thermo.Species.Count, thermo.Species.MoleFractions.Length);
        Assert.Equal(thermo.Species.Count, thermo.Species.MassFractions.Length);
    }

    [Fact]
    public void ThermoPhase_PairsSet()
    {
        using var thermo = Application.CreateThermoPhase("gri30.yaml");

        Assert.Equal(101325, thermo.Pressure, 6);
        Assert.Equal(300, thermo.Temperature, 6);
        Assert.Equal(0.081894, thermo.Density, 6);

        thermo.SetPair(ThermoPair.TemperaturePressure, 350, 113000);

        Assert.Equal(113000, thermo.Pressure, 6);
        Assert.Equal(350, thermo.Temperature, 6);
        Assert.Equal(0.078283, thermo.Density, 6);

        thermo.SetPair(ThermoPair.DensityPressure, .094, 110000);

        Assert.Equal(110000, thermo.Pressure, 6);
        Assert.Equal(283.740398, thermo.Temperature, 6);
        Assert.Equal(0.094, thermo.Density, 6);
    }

    [Fact]
    public void ThermoPhase_ThrowsOnBadPair()
    {
        using var thermo = Application.CreateThermoPhase("gri30.yaml");

        Assert.Throws<ArgumentOutOfRangeException>(() => thermo.SetPair((ThermoPair)(-1), 0, 0));
    }
}
