// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

using Cantera.Interop;

namespace Cantera;

/// <summary>
/// Represents a thermodynamic phase.
/// </summary>
public partial class ThermoPhase
{
    readonly SolutionHandle _solutionHandle;

    /// <summary>
    /// Represents a func that sets a pair of thermo variables.
    /// </summary>
    delegate int SetPairFunc(ThermoPhaseHandle n, double val0, double val1);

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
                        .SingleOrDefault(m => m.Name == "thermo_setState_" + f.Name)
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

    ThermoPhase(SolutionHandle solutionHandle)
    {
        _solutionHandle = solutionHandle;
        _handle = LibCantera.sol_thermo(_solutionHandle);

        _species = new(() => new SpeciesCollection(_handle));
    }

    /// <summary>
    /// Returns a new <see cref="ThermoPhase" /> object by loading and parsing the
    /// given configuration file. Optionally chooses the phase to load by
    /// looking up the given name.
    /// </summary>
    public static ThermoPhase Load(string filename, string? phaseName = null) =>
        new(LibCantera.sol_newSolution(filename, phaseName ?? "", "none"));

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

        LibCantera.thermo_equilibrate(_handle, interopString, solver,
            tolerance, maxSteps, maxIterations, logVerbosity);
    }

    /// <summary>
    /// Sets the given pair of thermodynamic properties for this phase together.
    /// </summary>
    public void SetPair(ThermoPair pair, double first, double second)
    {
        if (!s_pairSetters.Value.TryGetValue(pair, out var setter))
        {
            throw new InvalidOperationException($"Cannot set thermo pair {pair}!");
        }

        setter(_handle, first, second);
    }

    partial void ExtraDispose()
    {
        _solutionHandle.Dispose();
    }
}
