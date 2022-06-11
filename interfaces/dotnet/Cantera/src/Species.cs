namespace Cantera;

public class Species
{
    public string Name { get; }
    public double MolecularWeight { get; }
    public double MoleFraction { get; }
    public double MassFraction { get; }
    public double Charge { get; }

    internal Species(string name, double molecularWeight, double moleFraction, double massFraction, double charge)
    {
        Name = name;
        MolecularWeight = molecularWeight;
        MoleFraction = moleFraction;
        MassFraction = massFraction;
        Charge = charge;
    }
}