// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Represents a thermodynamic phase.
/// </summary>
public partial class ThermoPhase
{
    readonly SolutionHandle _sol;

    /// <summary>
    /// Represents a func that sets a pair of thermo variables using a pointer
    /// to a pair of doubles to stand in for a stack-allocated array with two elements.
    /// </summary>
    unsafe delegate int SetPairFunc(ThermoPhaseHandle n, (double, double)* values);

    /// <summary>
    /// Using reflection and the fact that CLib follows a naming convention for
    /// the functions that set the pairs of thermodynamic variables simultaneously
    /// </summary>
    static readonly Lazy<Dictionary<ThermoPair, SetPairFunc>> s_pairSetters;

    static ThermoPhase()
    {
        s_pairSetters = new(() =>
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

    /// <summary>
    /// The collection of species that make up this phase.
    /// </summary>
    public SpeciesCollection Species => _species.Value;

    internal ThermoPhase(string filename, string? phaseName)
    {
        _sol = LibCantera.soln_newSolution(filename, phaseName ?? "", "none");
        _handle = LibCantera.soln_thermo(_sol);
        _handle.EnsureValid();

        _species = new(() => new SpeciesCollection(_handle));
    }

    /// <summary>
    /// Simulates bringing the phase to thermodynamic equilibrium by holding the
    /// specified <see cref="ThermoPair" /> constant and using the algorithm(s)
    /// identified by the solver string.
    /// </summary>
    public void Equilibrate(ThermoPair thermoPair,
                            string solver = "auto",
                            double tolerance = 1e-9, int maxSteps = 1000,
                            int maxIterations = 100, int logVerbosity = 0)
    {
        var interopString = thermoPair.ToInteropString();

        var retVal = LibCantera.thermo_equilibrate(_handle, interopString, solver,
            tolerance, maxSteps, maxIterations, logVerbosity);

        InteropUtil.CheckReturn(retVal);
    }

    /// <summary>
    /// Sets the given pair of thermodynamic properties for this phase together.
    /// </summary>
    public unsafe void SetPair(ThermoPair pair, double first, double second)
    {
        if (!s_pairSetters.Value.TryGetValue(pair, out var setter))
        {
            throw new InvalidOperationException($"Cannot set thermo pair {pair}!");
        }

        var tuple = (first, second);

        InteropUtil.CheckReturn(setter(_handle, &tuple));
    }
}
