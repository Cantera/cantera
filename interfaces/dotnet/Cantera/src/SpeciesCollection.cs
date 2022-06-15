using System.Buffers;
using System.Collections;
using Cantera.Interop;

namespace Cantera;

public class SpeciesCollection : IReadOnlyList<Species>
{
    readonly ThermoPhaseHandle _handle;
    readonly List<Species> _species;

    public Species this[int index] => _species[index];

    public int Count => _species.Count;

    public unsafe double[] MassFractions
    {
        get => InteropUtil.GetDoubles(_handle, _species.Count, LibCantera.thermo_getMassFractions);
        set => InteropUtil.CheckReturn(
            LibCantera.thermo_setMassFractions(_handle, (nuint) _species.Count, value, InteropConsts.True));
    }

    public unsafe double[] MoleFractions
    {
        get => InteropUtil.GetDoubles(_handle, _species.Count, LibCantera.thermo_getMoleFractions);
        set => InteropUtil.CheckReturn(
            LibCantera.thermo_setMoleFractions(_handle, (nuint) _species.Count, value, InteropConsts.True));
    }

    internal unsafe SpeciesCollection(ThermoPhaseHandle handle)
    {
        _handle = handle;

        var count = (int) InteropUtil.CheckReturn(LibCantera.thermo_nSpecies(handle));

        var pool = MemoryPool<double>.Shared;

        using var __1 = pool.Rent(count, out var molecularWeights);
        using var __2 = pool.Rent(count, out var charges);

        InteropUtil.GetDoubles(handle, molecularWeights, LibCantera.thermo_getMolecularWeights);
        InteropUtil.GetDoubles(handle, charges, LibCantera.thermo_getCharges);

        _species = new(count);

        for (var i = 0; i < count; i++)
        {
            var name = InteropUtil.GetString(10, (length, buffer) =>
                LibCantera.thermo_getSpeciesName(handle, (nuint) i, (nuint) length, buffer));

            _species.Add(new
            (
                name,
                molecularWeights[i],
                charges[i]
            ));
        }
    }

    public IEnumerator<Species> GetEnumerator() =>
        _species.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() =>
        _species.GetEnumerator();

    public void SetUnormalizedMassFractions(double[] fractions) =>
        LibCantera.thermo_setMassFractions(_handle, (nuint) _species.Count, fractions, InteropConsts.False);

    public void SetUnormalizedMoleFractions(double[] fractions) =>
        LibCantera.thermo_setMoleFractions(_handle, (nuint) _species.Count, fractions, InteropConsts.False);
}