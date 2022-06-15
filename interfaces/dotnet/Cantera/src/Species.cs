namespace Cantera;

public class Species
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
}