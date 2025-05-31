// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using System.Buffers;
using System.Collections;
using Cantera.Interop;

namespace Cantera;

/// <summary>
/// The collection of species in a thermodynamic phase.
/// </summary>
public class SpeciesCollection : IReadOnlyList<Species>
{
    readonly ThermoPhaseHandle _handle;

    // this collection should be eagerly-initialized because it depends on interop
    readonly List<Species> _species;

    // this collection can be lazy-initialized because it relies only on other elements
    // from the above list
    readonly Lazy<Dictionary<string, int>> _speciesIndexByName;

    /// <summary>
    /// Gets a particular <see cref="Species" /> contained in this collection
    /// by it’s index.
    /// </summary>
    public Species this[int index] => _species[index];

    /// <summary>
    /// Gets a particular <see cref="Species" /> contained in this collection
    /// by it’s name.
    /// </summary>
    public Species this[string name] => _species[_speciesIndexByName.Value[name]];

    /// <summary>
    /// The number of <see cref="Species" /> in this collection.
    /// </summary>
    public int Count => _species.Count;

    /// <summary>
    /// Gets or sets the mass fractions for all the species at once,
    /// in the order in which they appear in this collection. When setting,
    /// normalizes the mass fractions so they sum to 1.
    /// </summary>
    public double[] MassFractions
    {
        get => InteropUtil.GetDoubles(_handle, _species.Count,
            LibCantera.thermo_getMassFractions);

        set => LibCantera.thermo_setMassFractions(_handle,  value.Length, value);
    }

    /// <summary>
    /// Gets or sets the mole fractions for all the species at once,
    /// in the order in which they appear in this collection. When setting,
    /// normalizes the mole fractions so they sum to 1.
    /// </summary>
    public double[] MoleFractions
    {
        get => InteropUtil.GetDoubles(_handle, _species.Count,
            LibCantera.thermo_getMoleFractions);

        set => LibCantera.thermo_setMoleFractions(_handle, value.Length, value);
    }

    internal SpeciesCollection(ThermoPhaseHandle handle)
    {
        _handle = handle;

        var count = LibCantera.thermo_nSpecies(handle);

        var pool = MemoryPool<double>.Shared;

        using (pool.Rent(count, out var molecularWeights))
        using (pool.Rent(count, out var charges))
        {
            InteropUtil.GetDoubles(handle, molecularWeights,
                LibCantera.thermo_getMolecularWeights);
            InteropUtil.GetDoubles(handle, charges, LibCantera.thermo_getCharges);

            _species = new(count);

            for (var i = 0; i < count; i++)
            {
                int getName(int length, Span<byte> buffer) => LibCantera
                    .thermo_speciesName(handle, i, length, buffer);

                var name = InteropUtil.GetString(10, getName);

                _species.Add(new
                (
                    name,
                    molecularWeights[i],
                    charges[i]
                ));
            }
        }

        _speciesIndexByName = new(() => _species
            .Select((s, i) => (species: s, index: i))
            .ToDictionary(s => s.species.Name, s => s.index,
                StringComparer.OrdinalIgnoreCase));
    }

    /// <inheritdoc />
    public IEnumerator<Species> GetEnumerator() =>
        _species.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() =>
        _species.GetEnumerator();

    /// <summary>
    /// Determines the index of the given <see cref="Species" /> in this collection,
    /// returning -1 if the species is not found.
    /// </summary>
    public int IndexOf(Species species) =>
        IndexOf(species.Name);

    /// <summary>
    /// Determines the index of a <see cref="Species" /> with the given name in
    /// this collection, returning -1 if the no species with the is not found.
    /// </summary>
    public int IndexOf(string name)
    {
        if (_speciesIndexByName.Value.TryGetValue(name, out var index))
        {
            return index;
        }

        return -1;
    }

    /// <summary>
    /// Reports whether the given <see cref="Species" /> is found in this collection.
    /// </summary>
    public bool Contains(Species species) =>
        Contains(species.Name);


    /// <summary>
    /// Reports whether a <see cref="Species" /> with the given name is
    /// found in this collection.
    /// </summary>
    public bool Contains(string name) =>
        _speciesIndexByName.Value.ContainsKey(name);

    /// <summary>
    /// Sets the mass fractions for a subset of the species, given as pairs of
    /// their name and fraction. The fractions of all the species are then normalized.
    /// </summary>
    public void SetMassFractions(params (string speciesName, double fraction)[] tuples)
    {
        SetMassFractions((IEnumerable<(string, double)>)tuples);
    }

    /// <inheritdoc cref="SetMassFractions(ValueTuple{string, double}[])" />
    public void SetMassFractions(
                              IEnumerable<(string speciesName, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var (name, fraction) in tuples)
            fractions[EnsureIndexOf(name)] = fraction;

        MassFractions = fractions;
    }

    /// <summary>
    /// Sets the mass fractions for a subset of the species, given as pairs of
    /// the <see cref="Species" /> and its fraction. The fractions of all the species
    /// are then normalized.
    /// </summary>
    public void SetMassFractions(params (Species species, double fraction)[] tuples) =>
        SetMassFractions((IEnumerable<(Species, double)>) tuples);

    /// <inheritdoc cref="SetMassFractions(ValueTuple{Species, double}[])" />
    public void SetMassFractions(IEnumerable<(Species species, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var (species, fraction) in tuples)
            fractions[EnsureIndexOf(species)] = fraction;

        MassFractions = fractions;
    }

    /// <summary>
    /// Sets the mole fractions for a subset of the species, given as pairs of
    /// their name and fraction. The fractions of all the species are then normalized.
    /// </summary>
    public void SetMoleFractions(params (string speciesName, double fraction)[] tuples)
    {
        SetMoleFractions((IEnumerable<(string, double)>)tuples);
    }

    /// <inheritdoc cref="SetMoleFractions(ValueTuple{string, double}[])" />
    public void SetMoleFractions(
                              IEnumerable<(string speciesName, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var (speciesName, fraction) in tuples)
            fractions[EnsureIndexOf(speciesName)] = fraction;

        MoleFractions = fractions;
    }

    /// <summary>
    /// Sets the mole fractions for a subset of the species, given as pairs of
    /// the <see cref="Species" /> and its fraction. The fractions of all the species
    /// are then normalized.
    /// </summary>
    public void SetMoleFractions(params (Species species, double fraction)[] tuples) =>
        SetMoleFractions((IEnumerable<(Species, double)>) tuples);

    /// <inheritdoc cref="SetMoleFractions(ValueTuple{Species, double}[])" />
    public void SetMoleFractions(IEnumerable<(Species species, double fraction)> tuples)
    {
        var fractions = new double[_species.Count];

        foreach(var (species, fraction) in tuples)
            fractions[EnsureIndexOf(species)] = fraction;

        MoleFractions = fractions;
    }

    int EnsureIndexOf(Species species) =>
        EnsureIndexOf(species.Name);

    int EnsureIndexOf(string name)
    {
        if (_speciesIndexByName.Value.TryGetValue(name, out var index))
        {
            return index;
        }

        throw new InvalidOperationException(
            $"Unknown species! {name} is not present in this collection.");
    }
}
