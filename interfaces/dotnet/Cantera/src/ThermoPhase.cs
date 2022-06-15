using System.Reflection;
using System.Text.RegularExpressions;
using Cantera.Interop;

namespace Cantera;

public partial class ThermoPhase
{
    /// <summary>
    /// Using reflection and the fact the CLIB follows a naming convention for the functions that
    /// set the pairs of thermodynamic variables simultaneously
    /// </summary>
    static readonly Lazy<Dictionary<ThermoPair, Func<ThermoPhaseHandle, double[], int>>> PairSetters;

    static ThermoPhase()
    {
        PairSetters = new(() =>
        {
            var methods = typeof(LibCantera).GetMethods();

            var pairs = typeof(ThermoPair)
                .GetFields(BindingFlags.Static | BindingFlags.Public)
                .Where(f => Regex.IsMatch(f.Name, "^[A-Z]{2}$")) // match exactly two uppercase
                .Select(f => (pair: f, method: methods.SingleOrDefault(m => m.Name == "thermo_set_" + f.Name)))
                .Where(t => t.method is not null)
                .ToDictionary(
                    t => (ThermoPair) t.pair.GetValue(null)!,
                    t => (Func<ThermoPhaseHandle, double[], int>)
                        t.method!.CreateDelegate(typeof(Func<ThermoPhaseHandle, double[], int>)));

            return pairs;
        });
    }

    readonly Lazy<SpeciesCollection> _species;

    public SpeciesCollection Species => _species.Value;

    public ThermoPhase(string filename, string? phasename = null)
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

    public void SetPair(ThermoPair pair, double first, double second)
    {
        if (!PairSetters.Value.TryGetValue(pair, out var setter))
        {
            throw new InvalidOperationException($"Cannot set thermo pair {pair}!");
        }

        InteropUtil.CheckReturn(setter(_handle, new[] { first, second }));
    }
}