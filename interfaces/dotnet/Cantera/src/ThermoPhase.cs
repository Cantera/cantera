using Cantera.Interop;

namespace Cantera;

public partial class ThermoPhase
{
    readonly Lazy<SpeciesCollection> _species;

    public SpeciesCollection Species => _species.Value;

    public ThermoPhase(string filename, string? phasename = null)
    {
        _handle = LibCantera.thermo_newFromFile(filename, phasename ?? "");
        _handle.EnsureValid();

        _species = new(() => new SpeciesCollection(_handle));
    }
}