using Xunit;

namespace Cantera.Tests;

public class ThermoPhaseTest
{
    static readonly string InputFile = Path.GetFullPath("gri30.yaml");

    [Fact]
    public void ThermoPhase_SpeciesRetrieved()
    {
        using var thermo = new ThermoPhase(InputFile);

        Assert.NotEmpty(thermo.Species);
        Assert.NotEmpty(thermo.Species[0].Name);
        Assert.Equal(thermo.Species.Count, thermo.Species.MoleFractions.Length);
        Assert.Equal(thermo.Species.Count, thermo.Species.MassFractions.Length);
    }
} 