// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;

namespace Cantera;

public partial class ThermoPhase
{
    /// <summary>
    /// Represents a func the sets a pair of thermo variables using a pointer
    /// to a pair of doubles to stand in for a stack-allocated array with two elements.
    /// </summary>
    unsafe delegate int SetPairFunc(ThermoPhaseHandle n, (double, double)* vals);

    /// <summary>
    /// Using reflection and the fact the CLIB follows a naming convention for
    /// the functions that set the pairs of thermodynamic variables simultaneously
    /// </summary>
    static readonly Lazy<Dictionary<ThermoPair, SetPairFunc>> PairSetters;

    static ThermoPhase()
    {
        PairSetters = new(() =>
        {
            var methods = typeof(LibCantera).GetMethods();

            var pairs = ThermoPairExtensions.GetThermoPairEnumFieldsWithTwoCharName()
                .Select(f =>
                (
                    pair: f,
                    method: methods
                        .SingleOrDefault(m => m.Name == "thermo_set_" + f.Name)
                ))
                .Where(t => t.method is not null)
                .ToDictionary(
                    t => (ThermoPair) t.pair.GetValue(null)!,
                    t => (SetPairFunc) t.method!.CreateDelegate(typeof(SetPairFunc)));

            return pairs;
        });
    }

    readonly Lazy<SpeciesCollection> _species;

    public SpeciesCollection Species => _species.Value;

    internal ThermoPhase(string filename, string? phasename)
    {
        _handle = LibCantera.thermo_newFromFile(filename, phasename ?? "");
        _handle.EnsureValid();

        _species = new(() => new SpeciesCollection(_handle));
    }

    public void Equilibrate(ThermoPair thermoPair,
                            EquilibriumSolver solver = EquilibriumSolver.Auto,
                            double tolerance = 1e-9, int maxSteps = 1000,
                            int maxIterations = 100, bool logProgress = false)
    {
        var interopString = thermoPair.ToInteropString();
        var interopLogProgress = InteropUtil.GetInteropBool(logProgress);

        var retVal = LibCantera.thermo_equilibrate(_handle, interopString, (int) solver,
            tolerance, maxSteps, maxIterations, interopLogProgress);

        InteropUtil.CheckReturn(retVal);
    }

    public unsafe void SetPair(ThermoPair pair, double first, double second)
    {
        if (!PairSetters.Value.TryGetValue(pair, out var setter))
        {
            throw new InvalidOperationException($"Cannot set thermo pair {pair}!");
        }

        var tuple = (first, second);

        InteropUtil.CheckReturn(setter(_handle, &tuple));
    }
}
