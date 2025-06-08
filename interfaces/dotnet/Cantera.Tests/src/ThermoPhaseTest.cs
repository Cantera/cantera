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
}
