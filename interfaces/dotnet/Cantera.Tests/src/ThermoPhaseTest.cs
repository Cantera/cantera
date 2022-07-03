using Xunit;

namespace Cantera.Tests;

[Collection("Application")]
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
}
