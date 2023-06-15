// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

namespace Cantera;

/// <summary>
/// Represents a chemical species.
/// </summary>
public class Species : IEquatable<Species>
{
#pragma warning disable CS1591
    public string Name { get; }
    public double MolecularWeight { get; }
    public double Charge { get; }
#pragma warning restore CS1591

    internal Species(string name, double molecularWeight, double charge)
    {
        Name = name;
        MolecularWeight = molecularWeight;
        Charge = charge;
    }

    /// <inheritdoc />
    public bool Equals(Species? other) =>
        StringComparer.OrdinalIgnoreCase.Equals(Name, other?.Name);

    /// <inheritdoc />
    public override bool Equals(object? obj) =>
        obj is Species other && Equals(other);

    /// <inheritdoc />
    public override int GetHashCode() =>
        StringComparer.OrdinalIgnoreCase.GetHashCode(Name);
}
