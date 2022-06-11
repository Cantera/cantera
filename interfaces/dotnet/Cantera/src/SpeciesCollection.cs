using System.Collections;
using Cantera.Interop;

namespace Cantera;

public class SpeciesCollection : IReadOnlyList<Species>
{
    readonly List<Species> _species;

    public Species this[int index] => _species[index];

    public int Count => _species.Count;

    internal unsafe SpeciesCollection(ThermoPhaseHandle handle)
    {
        var count = InteropUtil.CheckReturn(LibCantera.thermo_nSpecies(handle));

        var molecularWeights = InteropUtil.GetDoubles(handle, count, LibCantera.thermo_getMolecularWeights);
        var moleFractions = InteropUtil.GetDoubles(handle, count, LibCantera.thermo_getMoleFractions);
        var massFractions = InteropUtil.GetDoubles(handle, count, LibCantera.thermo_getMassFractions);
        var charges = InteropUtil.GetDoubles(handle, count, LibCantera.thermo_getCharges);

        _species = new((int) count);

        for (var i = (nuint) 0; i < count; i++)
        {
            var name = InteropUtil.GetString(10, (length, buffer) =>
                LibCantera.thermo_getSpeciesName(handle, i, (nuint) length, buffer));

            _species.Add(new
            (
                name,
                molecularWeights[i],
                moleFractions[i],
                massFractions[i],
                charges[i]
            ));
        }
    }

    public IEnumerator<Species> GetEnumerator() =>
        _species.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() =>
        _species.GetEnumerator();
}