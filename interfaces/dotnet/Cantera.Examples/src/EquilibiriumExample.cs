namespace Cantera.Examples;

class EquilibriumExample : IExample
{
    public void Run()
    {
        Application.DataDirectories.AddAssemblyDirectory();

        Application.AddConsoleLogging();

        var phase = Application.CreateThermoPhase("gri30.yaml");
        var species = phase.Species;

        species.SetMoleFractions(("CH4", 1.00), ("O2", 2.20), ("N2", 7.52));

        for (var n = 0; n < 27; n++)
        {
            var T = 300.0 + 100.0 * n;
            equilSoundSpeeds(phase, T);
        }
    }

    void equilSoundSpeeds(ThermoPhase gas, double T, double rtol = 1.0e-6, int maxIter=5000)
    {
        gas.SetPair(ThermoPair.TemperaturePressure, T, Consts.OneAtm);

        gas.Equilibrate(ThermoPair.TemperaturePressure, tolerance: rtol, maxIterations: maxIter);

        var s0 = gas.MassEntropy;
        var p0 = gas.Pressure;
        var r0 = gas.Density;

        var p1 = p0 * 1.0001;

        gas.SetPair(ThermoPair.EntropyPressure, s0, p1);

        var aFrozen = Math.Sqrt((p1 - p0)/(gas.Density - r0));

        gas.Equilibrate(ThermoPair.EntropyPressure, tolerance: rtol, maxIterations: maxIter);

        var aEquil = Math.Sqrt((p1 - p0)/(gas.Density - r0));

        var gamma = gas.MassCp / gas.MassCv;

        var aFrozen2 = Math.Sqrt(gamma * Consts.GasConstant * gas.Temperature / gas.MeanMolecularWeight);

        Console.WriteLine($"T: {T}, aEquil: {aEquil}, aFrozen: {aFrozen}, aFrozen2: {aFrozen2}");
    }
}
