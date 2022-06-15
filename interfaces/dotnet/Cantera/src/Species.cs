namespace Cantera;

public class Species : IEquatable<Species>
{
    public string Name { get; }
    public double MolecularWeight { get; }
    public double Charge { get; }

    internal Species(string name, double molecularWeight, double charge)
    {
        Name = name;
        MolecularWeight = molecularWeight;
        Charge = charge;
    }

    public bool Equals(Species? other) =>
        Name == other?.Name;

    public override bool Equals(object? obj) =>
        obj is Species other && Equals(other);

    public override int GetHashCode() =>
        Name.GetHashCode();
}
