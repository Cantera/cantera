using System.Buffers;
using System.Collections;
using Cantera.Interop;

namespace Cantera;

public class SpeciesCollection : IReadOnlyList<Species>
{
    readonly ThermoPhaseHandle _handle;

    // this collection should be eargerly-initialized because it depends on interop
    readonly List<Species> _species;

    // this collection can be lazy-initialzied because it relies only on other elements
    // from the above list
    readonly Lazy<Dictionary<string, Species>> _speciesByName;

    public Species this[int index] => _species[index];

    public Species this[string name] => _speciesByName.Value[name];

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

        using (pool.Rent(count, out var molecularWeights))
        using (pool.Rent(count, out var charges))
        {
            InteropUtil.GetDoubles(handle, molecularWeights, LibCantera.thermo_getMolecularWeights);
            InteropUtil.GetDoubles(handle, charges, LibCantera.thermo_getCharges);

            _species = new(count);

            for (var i = 0; i < count; i++)
            {
                var name = InteropUtil.GetString(10, (length, buffer) =>
                    LibCantera.thermo_getSpeciesName(handle, (nuint)i, (nuint)length, buffer));

                _species.Add(new
                (
                    name,
                    molecularWeights[i],
                    charges[i]
                ));
            }
        }

        _speciesByName = new(() => _species.ToDictionary(s => s.Name, s => s));
    }

    public IEnumerator<Species> GetEnumerator() =>
        _species.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() =>
        _species.GetEnumerator();

    public int IndexOf(string name) =>
        _species.FindIndex(s => s.Name == name);

    public int IndexOf(Species species) =>
        _species.IndexOf(species);

    public bool Contains(string name) =>
        _speciesByName.Value.ContainsKey(name);

    public bool Contains(Species species) =>
        _species.Contains(species);

    public void SetMassFractions(params (string speciesName, double fraction)[] tuples) =>
        SetMassFractions((IEnumerable<(string, double)>) tuples);

    public void SetMassFractions(IEnumerable<(string speciesName, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var tuple in tuples)
            fractions[EnsureIndexOf(tuple.speciesName)] = tuple.fraction;

        MassFractions = fractions;
    }

    public void SetMassFractions(params (Species species, double fraction)[] tuples) =>
        SetMassFractions((IEnumerable<(Species, double)>) tuples);

    public void SetMassFractions(IEnumerable<(Species species, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var tuple in tuples)
            fractions[EnsureIndexOf(tuple.species)] = tuple.fraction;

        MassFractions = fractions;
    }

    public void SetMoleFractions(params (string speciesName, double fraction)[] tuples) =>
        SetMoleFractions((IEnumerable<(string, double)>) tuples);

    public void SetMoleFractions(IEnumerable<(string speciesName, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var tuple in tuples)
            fractions[EnsureIndexOf(tuple.speciesName)] = tuple.fraction;

        MoleFractions = fractions;
    }

    public void SetMoleFractions(params (Species species, double fraction)[] tuples) =>
        SetMoleFractions((IEnumerable<(Species, double)>) tuples);

    public void SetMoleFractions(IEnumerable<(Species species, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var tuple in tuples)
            fractions[EnsureIndexOf(tuple.species)] = tuple.fraction;

        MoleFractions = fractions;
    }

    int EnsureIndexOf(Species species)
    {
        var index = IndexOf(species);

        if (index == -1)
        {
            throw new InvalidOperationException(
                $"Unknown species! {species.Name} is not present in this collection.");
        }

        return index;
    }

    int EnsureIndexOf(string name)
    {
        var index = IndexOf(name);

        if (index == -1)
        {
            throw new InvalidOperationException(
                $"Unknown species! {name} is not present in this collection.");
        }

        return index;
    }

    public void SetUnormalizedMassFractions(double[] fractions) =>
        LibCantera.thermo_setMassFractions(_handle, (nuint) _species.Count, fractions, InteropConsts.False);

    public void SetUnormalizedMoleFractions(double[] fractions) =>
        LibCantera.thermo_setMoleFractions(_handle, (nuint) _species.Count, fractions, InteropConsts.False);
}
