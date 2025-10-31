/*
Cantera C# Examples
Application example

Computes the “equilibrium” and “frozen” sound speeds for a gas.
*/

using Cantera;

Application.DataDirectories.AddAssemblyDirectory();

Application.AddConsoleLogging();

var phase = Application.CreateThermoPhase("gri30.yaml");
var species = phase.Species;

species.SetMoleFractions(("CH4", 1.00), ("O2", 2.00), ("N2", 7.52));

for (var n = 0; n < 27; n++)
{
    var T = 300.0 + 100.0 * n;
    var (aEquil, aFrozen, aFrozen2) = EquilSoundSpeeds(phase, T);

    Console.WriteLine(
        $"T: {T}, aEquil: {aEquil}, aFrozen: {aFrozen}, aFrozen2: {aFrozen2}");
}

/// <summary>
/// Returns a tuple containing the equilibrium and frozen sound speeds for a gas
/// with an equilibrium composition. The gas is first set to an equilibrium state
/// at the temperature and pressure of the gas, because the equilibrium sound speed
/// is otherwise undefined.
/// </summary>
static (double equil, double frozen, double frozen2)
    EquilSoundSpeeds(ThermoPhase gas, double T, double rTol = 1.0e-6,
                     int maxIter = 5000)
{
    gas.SetPair(ThermoPair.TemperaturePressure, T, Consts.OneAtm);

    // Set the gas to equilibrium at its current T and P.
    gas.Equilibrate(
        ThermoPair.TemperaturePressure, tolerance: rTol, maxIterations: maxIter);

    // Save properties.
    var s0 = gas.MassEntropy;
    var p0 = gas.Pressure;
    var r0 = gas.Density;

    // Perturb the pressure.
    var p1 = p0 * 1.0001;

    // Set the gas to a state with the same entropy and composition but
    // the perturbed pressure.
    gas.SetPair(ThermoPair.EntropyPressure, s0, p1);

    // frozen sound speed
    var aFrozen = Math.Sqrt((p1 - p0)/(gas.Density - r0));

    // Now equilibrate the gas holding S and P constant.
    gas.Equilibrate(
        ThermoPair.EntropyPressure, tolerance: rTol, maxIterations: maxIter);

    // equilibrium sound speed
    var aEquil = Math.Sqrt((p1 - p0)/(gas.Density - r0));

    // Compute the frozen sound speed using the ideal gas expression as a check.
    var gamma = gas.MassCp / gas.MassCv;
    var aFrozen2 = Math.Sqrt(
        gamma * Consts.GasConstant * gas.Temperature / gas.MeanMolecularWeight);

    return (aEquil, aFrozen, aFrozen2);
}
